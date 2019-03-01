from pytriqs.gf.local import *
from pytriqs.operators import *
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi
from pytriqs.applications.impurity_solvers.cthyb import Solver

# Parameters
U = 4.0         # Local Coulomb repulsion
e_f = -U / 2. - 2   # impurity energy level (= -U/2 at half filling)
V = 0.5         # bath-impurity hopping
D = 1.0         # half-width of the bath spectrum
beta = 5       # inverse temperature

def SIAM(U, e_f, V, D, beta, filename="qmc_results.h5"):
    Delta = V**2 * Flat(D)
    N_MC = 1e5
    l_max = 10
    independent_samples = 16
    for l in range(l_max + 1):
        for i in range(independent_samples):
            S = Solver(beta=beta, gf_struct={'up': [0], 'down': [0]})
            # Initialize the non-interacting Green's function S.G0_iw
            for name, g0 in S.G0_iw: g0 << inverse(iOmega_n - e_f - Delta)
            # Run the solver. The results will be in S.G_tau, S.G_iw and S.G_l
            S.solve(h_int = U * n('up',0) * n('down',0),        # Local Hamiltonian
                    n_cycles  = int(N_MC / 2**l),                              # Number of QMC cycles
                    length_cycle = 2**l,                        # Length of one cycle
                    n_warmup_cycles = int(N_MC / 2**l / 100),          #  Warmup cycles
                    measure_g_tau = False,                      #  Don't measure G_tau
                    measure_g_l = False,                        #  Don't measure G_l
                    perform_post_proc = False,                  #  Don't measure G_iw
                    use_norm_as_weight=True,                    # Necessary option for the measurement of the density matrix
                    measure_density_matrix=True,                # Measure reduced impurity density matrix
                    random_seed = i * 8521 + l * 14187 + mpi.rank * 7472)     # Random seed, very important!
            # Save the results in an HDF5 file (only on the master node)
            if mpi.is_master_node():
                with HDFArchive(filename) as Results:
                    Results["rho_l{}_i{}".format(l, i)] = S.density_matrix


if __name__ == "__main__":
    SIAM(U, e_f, V, D, beta)