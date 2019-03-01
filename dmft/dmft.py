from pytriqs.gf.local import *
from pytriqs.operators import *
from pytriqs.archive import HDFArchive
import pytriqs.utility.mpi as mpi
from pytriqs.applications.impurity_solvers.cthyb import Solver

# Parameters
U = 5.0         # Local Coulomb repulsion
e_d = -U / 2.   # energy level (= -U/2 at half filling)
t = 1.          # Bethe lattice hopping
beta = 10       # inverse temperature

# DMFT parameters

def DMFT(U, e_d, t, beta, filename="dmft_results.h5"):
    # Construct the CT-HYB-QMC solver
    S = Solver(beta = beta, gf_struct = {'up':[0], 'down':[0]}, n_l = 50)

    # Initialize Delta
    Delta = GfImFreq(beta = beta, indices = [0])
    Delta << t ** 2 * SemiCircular(half_bandwidth = 2 * t)

    # Now do the DMFT loop
    n_iter = 8
    for iter in range(n_iter):

        # Compute new S.G0_iw
        for name, g0 in S.G0_iw:
            g0 << inverse(iOmega_n - e_d - Delta)
        # Run the solver
        S.solve(h_int = U * n('up',0) * n('down',0),    # Local Hamiltonian
                n_cycles = 200000,                      # Number of QMC cycles
                length_cycle = 50,                     # Length of a cycle
                n_warmup_cycles = 2000,                 # How many warmup cycles
                measure_g_l = True)
        # Compute new Delta with the self-consistency condition while imposing paramagnetism
        g_l = (S.G_l['up'] + S.G_l['down']) / 2.
        Delta.set_from_legendre(t**2 * g_l)

        # Intermediate saves
        if mpi.is_master_node():
            with HDFArchive(filename) as Results:
                Results["G_tau_iter{}".format(iter)] = S.G_tau
                Results["G_iw_iter{}".format(iter)] = S.G_iw
                Results["G_l_iter{}".format(iter)] = S.G_l
                Results["Sigma_iter{}".format(iter)] = S.Sigma_iw

    if mpi.is_master_node():
        with HDFArchive(filename) as Results:
            Results["G_tau"] = S.G_tau
            Results["G_iw"] = S.G_iw
            Results["G_l"] = S.G_l
            Results["Sigma"] = S.Sigma_iw

if __name__ == "__main__":
    DMFT(U, e_d, t, beta)
