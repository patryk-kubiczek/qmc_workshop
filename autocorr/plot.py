import matplotlib
matplotlib.use('TKAgg')
from matplotlib import pyplot as plt
from pytriqs.archive import HDFArchive
from math import pi, sqrt
import numpy as np

def make_diagonal_rho(rho):
    return [np.array(el[0][0]) for el in rho]

def load_from_QMC(filename, l_max = 10, n_samples = 16):
    rho = []
    for l in range(l_max + 1):
        rho.append([])
        for i in range(n_samples):
            with HDFArchive(filename,'r') as QMC:
                rho[-1].append(make_diagonal_rho(QMC["rho_l{}_i{}".format(l, i)]))
    return rho


if __name__ == "__main__":

    rho = load_from_QMC("qmc_results.h5")

    l_max = len(rho) - 1
    rho_av = []
    rho_std = []

    for rho_l in rho:
        rho_av.append(np.average(rho_l, axis=0))
        rho_std.append(np.std(rho_l, axis=0))
        #rho_std.append(np.sqrt(np.average([(x - rho_av[-1])**2 for x in rho_l], axis=0)))



    plt.plot([l for l in range(l_max + 1)], [r[2] for r in rho_av], 'o')
    plt.xlabel("$\mathrm{log}_2(\mathrm{cycle \,length})$")
    plt.ylabel("$d$")
    plt.savefig("d.png")
    plt.close()

    plt.plot([l for l in range(l_max + 1)], [np.average(r) for r in rho_std], 'o')
    plt.xlabel("$\mathrm{log}_2(\mathrm{cycle \, length})$")
    plt.ylabel("$\\Delta \\rho$")
    plt.savefig("deltad.png")
    plt.close()





