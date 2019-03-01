import matplotlib
matplotlib.use('TKAgg')
from pytriqs.plot.mpl_interface import *
from pytriqs.gf.local import *
from pytriqs.statistics.histograms import *
from pytriqs.archive import HDFArchive
from math import pi
import numpy as np

def load_from_QMC(filename):
    with HDFArchive(filename,'r') as QMC:
        G_tau = QMC["G_tau"]
        G_iw = QMC["G_iw"]
        G_l = QMC["G_l"]
        rho = QMC["rho"]
        k_histogram = QMC["k_histogram"]
        average_sign = QMC["average_sign"]
    return G_tau, G_iw, G_l, rho, k_histogram, average_sign

def load_from_SOM(filename):
    with HDFArchive(filename, 'r') as SOM:
        g_w = SOM["g_w"]
        g_rec_tau = SOM["g_rec_tau"]
    return g_w, g_rec_tau

if __name__ == "__main__":

    G_tau, G_iw, G_l, rho, k_histogram, average_sign = load_from_QMC("qmc_results.h5")

    print "Average sign:", average_sign

    diag_rho = [el[0][0] for el in rho]
    print "Diagonal of the reduced density matrix (up, down, double, 0):",  diag_rho

    max_k = 10
    norm = np.sum(k_histogram.data)
    plt.bar(range(max_k), k_histogram.data[:max_k] / norm,  align='center')
    plt.ylabel("$P(k)$")
    plt.xlim([-1, max_k + 1])
    plt.xlabel("$k$")
    plt.savefig("k_histogram.png")
    plt.close()

    oplotr(G_tau, '-')
    plt.legend(loc="best")
    plt.savefig("g_tau.png")
    plt.close()

    oplot(G_iw, "o", x_window=(0, 10))
    plt.legend(loc="best")
    plt.savefig("g_iw.png")
    plt.close()

    oplotr(G_l, "o")
    plt.legend(loc="best")
    plt.savefig("g_l.png")
    plt.close()

    g_w, g_rec_tau = load_from_SOM("som_results.h5")

    oplotr(G_tau, '-')
    oplotr(g_rec_tau, '-')
    plt.legend(loc="best")
    plt.ylabel("$G(\\tau)$")
    plt.savefig("g_rec_tau.png")
    plt.close()

    oploti(-1 / pi * g_w)
    plt.legend(loc="best")
    plt.ylabel("$A(\\omega)$")
    plt.savefig("spectrum.png")
    plt.close()

