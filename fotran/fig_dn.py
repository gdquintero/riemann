"""Figura 3 del paper: d_n(1/2)^2 (riemann_dn, dn_f90.txt) y s_n(1/2) exacto
(riemann_traj, traj_sn05/output_050.txt) para 2 <= n <= 5120, en escala log-log.
Escribe ../latex/fig_dn.pdf e imprime los valores de la Tabla 2.
"""
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
D = np.loadtxt(os.path.join(here, "dn_f90.txt"))
S = np.loadtxt(os.path.join(here, "traj_sn05", "output_050.txt"))
n_d, d2 = D[:, 0], D[:, 1]
n_s, sn = S[:, 0], S[:, 2]

plt.rcParams.update({"font.size": 9})
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

fig, ax = plt.subplots(figsize=(3.3, 2.7))
ax.loglog(n_s, sn, color="tab:red", lw=0.4, label=r"$s_n(1/2)$")
ax.loglog(n_d, d2, color="tab:blue", lw=1.0, label=r"$d_n(1/2)^2$")
ax.grid(True, which="major", linestyle="--", lw=0.4, alpha=0.6)
ax.tick_params(axis="both", which="both", direction="in")
ax.set_xlabel(r"$n$")
ax.legend(loc="lower left", fontsize=8, frameon=True, edgecolor="black", fancybox=False)
fig.savefig(os.path.join(here, "..", "latex", "fig_dn.pdf"), bbox_inches="tight")
plt.close(fig)

def sci(x):
    m, e = f"{x:.3e}".split("e")
    return f"${m}\\mathrm{{e}}{{{int(e)}}}$"

for n in [40, 160, 640, 1280, 2560, 5120]:
    print(f"${n}$ & {sci(d2[n_d == n][0])} & {sci(sn[n_s == n][0])} \\\\")
m = n_s >= 100
print(f"s_n(1/2), 100 <= n <= 5120: min {sn[m].min():.3e}, max {sn[m].max():.3e}, median {np.median(sn[m]):.3e}")
