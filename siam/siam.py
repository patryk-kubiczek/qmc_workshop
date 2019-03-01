from pytriqs.gf.local import *
from pytriqs.operators import *
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi
from pytriqs.applications.impurity_solvers.cthyb import Solver

# Parameters
U = 5.0         # Local Coulomb repulsion
e_f = -U / 2.   # impurity energy level (= -U/2 at half filling)
V = 0.5         # bath-impurity hopping
D = 1.0         # half-width of the bath spectrum
beta = 20       # inverse temperature

def SIAM(U, e_f, V, D, beta, filename="qmc_results.h5"):
    # Create hybridization function
    Delta = V**2 * Flat(D)

    # Construct the impurity solver with the inverse temperature
    # and the structure of the Green's functions
    S = Solver(beta = beta, gf_struct = {'up':[0], 'down':[0]}, n_l = 50)

    # Initialize the non-interacting Green's function S.G0_iw
    for name, g0 in S.G0_iw: g0 << inverse(iOmega_n - e_f - Delta)

    # Run the solver. The results will be in S.G_tau, S.G_iw and S.G_l
    S.solve(h_int = U * n('up',0) * n('down',0),     # Local Hamiltonian
            n_cycles  = 2000000,                     # Number of QMC cycles
            length_cycle = 50,                       # Length of one cycle
            n_warmup_cycles = 20000,                 # Warmup cycles
            measure_g_l = True,                      # Measure G_l (representation of G in terms of Legendre polynomials)
            use_norm_as_weight=True,                 # Necessary option for the measurement of the density matrix
            measure_density_matrix=True,             # Measure reduced impurity density matrix
            measure_pert_order=True)                 # Measure histogram of k

    # Save the results in an HDF5 file (only on the master node)
    if mpi.is_master_node():
        with HDFArchive(filename,'w') as Results:
            Results["G_tau"] = S.G_tau
            Results["G_iw"] = S.G_iw
            Results["G_l"] = S.G_l
            Results["rho"] = S.density_matrix
            Results["k_histogram"] = S.perturbation_order_total
            Results["average_sign"] = S.average_sign

if __name__ == "__main__":
    SIAM(U, e_f, V, D, beta)
        

