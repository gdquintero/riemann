"""Figuras log-log de s_n(alpha) a partir de una corrida de riemann_traj (por defecto
traj_1e6/), con el estilo de ../plot_all_alphas05.py.

Escribe fig_traj_05.pdf (alpha = 0.50..0.59), fig_traj_09.pdf (alpha = 0.90..1.00) y
fig_traj_all.pdf (alpha = 0.50, 0.60, ..., 1.00).

Uso:  python fig_traj.py [directorio]
"""
import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
folder = os.path.join(here, sys.argv[1] if len(sys.argv) > 1 else "traj_1e6")

plt.rcParams.update({"font.size": 12, "figure.figsize": [10, 6]})
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

def figure(alphas, name):
    for a in alphas:
        f = os.path.join(folder, f"output_{a:03d}.txt")
        if not os.path.exists(f):
            continue
        D = np.loadtxt(f)
        plt.loglog(D[:, 0], D[:, 2], lw=0.6, label=rf"$\alpha = {a / 100:.2f}$")
    plt.grid(True, which="both", linestyle="--", lw=0.4, alpha=0.5)
    plt.tick_params(axis="both", which="both", direction="in")
    plt.xlabel(r"$n$", fontsize=16)
    plt.ylabel(r"$s_n$", fontsize=16)
    plt.legend(loc="lower left", fontsize=10, frameon=True, edgecolor="black", fancybox=False,
               ncol=2)
    plt.savefig(os.path.join(here, name), bbox_inches="tight")
    plt.close()
    print(name)

figure(range(50, 60), "fig_traj_05.pdf")
figure(range(90, 101), "fig_traj_09.pdf")
figure(range(50, 101, 10), "fig_traj_all.pdf")
