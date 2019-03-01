from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi
import numpy as np
from pytriqs.applications.analytical_continuation.som import Som

def cut_coefficients(g_l, n_remaining_coeffs):
    g_cut = GfLegendre(indices = [0], beta = g_l.beta, n_points = n_remaining_coeffs)
    g_cut.data[:,:,:] = g_l.data[:n_remaining_coeffs, :, :]
    g_cut.enforce_discontinuity(np.identity(g_cut.N1))
    return g_cut

n_tau = 500                 # Number of tau-slices for the input GF
n_l = 30                    # Number of Legendre polynomials for the input GF
n_w = 801                   # Number of energy slices for the solution
energy_window = (-6.0, 6.0)  # Energy window to search the solution in

# Parameters for Som.run()
run_params = {'energy_window' : energy_window}
# Verbosity level
run_params['verbosity'] = 2
# Number of particular solutions to accumulate
run_params['l'] = 2000
# Number of global updates
run_params['f'] = 100
# Number of local updates per global update
run_params['t'] = 50
# Accumulate histogram of the objective function values
run_params['make_histograms'] = False


def SOM(input_filename="dmft_results.h5", output_filename="som_results.h5"):

# Read G(\tau) from archive
# Could be G(i\omega_n) or G_l as well.
    with HDFArchive(input_filename,'r') as QMC:
        G_tau = QMC["G_tau"]
        G_iw = QMC["G_iw"]
        G_l = QMC["G_l"]

    # Paramagnetic case: average spin up and spin down GF
    g_tau = (G_tau['up'] + G_tau['down']) / 2.
    g_iw = (G_iw['up'] + G_iw['down']) / 2.
    g_l = (G_l['up'] + G_l['down']) / 2.


    # Prepare input data: reduce the number of \tau-slices from 10001 to n_tau
    # reduce the number of Legendre coefficients to n_l
    g_tau_rebinned = rebinning_tau(g_tau, n_tau)
    g_l_cut = cut_coefficients(g_l, n_l)


    # Set the weight function S to a constant (all points of G_tau are equally important)
    S_tau = g_tau_rebinned.copy()
    S_tau.data[:] = 1.0

    S_l = g_l_cut.copy()
    S_l.data[:] = 1.0

    # Construct a SOM object
    #cont = Som(g_tau_rebinned, S_tau, kind = "FermionGf")
    cont = Som(g_l_cut, S_l, kind = "FermionGf")

    # Run!
    cont.run(**run_params)

    # Create a real frequency GF obtained with SOM
    g_w = GfReFreq(window = energy_window, n_points = n_w, indices = [0])
    g_w << cont

    # G(\tau) reconstructed from the SOM solution
    g_rec_tau = g_tau_rebinned.copy()
    g_rec_tau << cont

    # On master node, save results to an archive
    if mpi.is_master_node():
        with HDFArchive(output_filename,'w') as Results:
            Results['g_rec_tau'] = g_rec_tau
            Results['g_w'] = g_w


if __name__ == "__main__":
    SOM()