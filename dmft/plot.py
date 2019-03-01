import matplotlib
matplotlib.use('TKAgg')
from pytriqs.plot.mpl_interface import *
from pytriqs.gf.local import *
from pytriqs.archive import HDFArchive
from math import pi

def load_from_DMFT(filename, n_iter = 0):
    with HDFArchive(filename,'r') as QMC:
        G_tau = QMC["G_tau"]
        G_iw = QMC["G_iw"]
        G_l = QMC["G_l"]
        Sigma_iw = QMC["Sigma"]
        G_iw_list = [QMC["G_iw_iter{}".format(i)] for i in range(n_iter)]
    return G_tau, G_iw, G_l, Sigma_iw, G_iw_list

def load_from_SOM(filename):
    with HDFArchive(filename, 'r') as SOM:
        g_w = SOM["g_w"]
        g_rec_tau = SOM["g_rec_tau"]
    return g_w, g_rec_tau

if __name__ == "__main__":

    G_tau, G_iw, G_l, Sigma_iw, G_iw_list = load_from_DMFT("dmft_results.h5", n_iter=8)

    oplotr(G_tau, '-')
    plt.legend(loc="best")
    plt.savefig("g_tau.png")
    plt.close()

    oplot(G_iw, "o", x_window=(0, 10))
    plt.legend(loc="best")
    plt.savefig("g_iw.png")
    plt.close()

    oplot(Sigma_iw, "o", x_window=(0, 10))
    plt.legend(loc="best")
    plt.savefig("sigma_iw.png")
    plt.close()

    oplotr(G_l, "o")
    plt.legend(loc="best")
    plt.savefig("g_l.png")
    plt.close()

    for i, g_iw in enumerate(G_iw_list):
        oploti(g_iw['up'], ".-", x_window=(0, 10), label="Iter {}".format(i))
    plt.legend(loc="best")
    plt.ylabel("$G(i\\omega_n)$")
    plt.savefig("g_iw_iter.png")
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

