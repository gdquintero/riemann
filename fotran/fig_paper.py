"""Figuras de la Seccion 4 del paper a partir de traj_1e6/ (corrida de riemann_traj con
riemann_traj_1e6.nml).  Escribe en ../latex:

  fig_sn.pdf        s_n(alpha) en escala log-log, alpha = 0.50, 0.60, ..., 1.00;
  fig_exponent.pdf  cambio de log10 max s_n y log10 min s_n por decada de n, para todos
                    los alphas, junto con la recta 1 - 2 alpha.
"""
import glob
import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

here = os.path.dirname(os.path.abspath(__file__))
folder = os.path.join(here, "traj_1e6")
latex = os.path.join(here, "..", "latex")

plt.rcParams.update({"font.size": 10})
plt.rc("text", usetex=True)
plt.rc("font", family="serif")

data = {}
for f in glob.glob(os.path.join(folder, "output_*.txt")):
    D = np.loadtxt(f)
    data[round(100 * D[0, 1])] = (D[:, 0], D[:, 2])

# ---- s_n(alpha) for alpha = 0.50, 0.60, ..., 1.00 ----
fig, ax = plt.subplots(figsize=(6.3, 3.6))
for a in range(50, 101, 10):
    n, s = data[a]
    ax.loglog(n, s, lw=0.5, label=rf"$\alpha = {a / 100:.2f}$")
ax.grid(True, which="major", linestyle="--", lw=0.4, alpha=0.6)
ax.tick_params(axis="both", which="both", direction="in")
ax.set_xlabel(r"$n$")
ax.set_ylabel(r"$s_n(\alpha)$")
ax.legend(loc="lower left", fontsize=8, ncol=3, frameon=True, edgecolor="black", fancybox=False)
fig.savefig(os.path.join(latex, "fig_sn.pdf"), bbox_inches="tight")
plt.close(fig)

# ---- change of log10 max / min of s_n per decade (1e3..1e4 -> 1e5..1e6) ----
alphas = np.array(sorted(data)) / 100
slope_max, slope_min = [], []
for a in sorted(data):
    n, s = data[a]
    d1, d3 = (n >= 1e3) & (n < 1e4), n >= 1e5
    slope_max.append(np.log10(s[d3].max() / s[d1].max()) / 2)
    slope_min.append(np.log10(s[d3].min() / s[d1].min()) / 2)

fig, ax = plt.subplots(figsize=(5.0, 3.2))
ax.plot(alphas, 1 - 2 * alphas, "k--", lw=0.8, label=r"$1 - 2\alpha$")
ax.plot(alphas, slope_max, "o", ms=3, mfc="none", label=r"$\max s_n$")
ax.plot(alphas, slope_min, "s", ms=2.5, label=r"$\min s_n$")
ax.grid(True, linestyle="--", lw=0.4, alpha=0.6)
ax.tick_params(axis="both", direction="in")
ax.set_xlabel(r"$\alpha$")
ax.set_ylabel(r"change of $\log_{10}$ per decade of $n$")
ax.legend(loc="lower left", fontsize=8, frameon=True, edgecolor="black", fancybox=False)
fig.savefig(os.path.join(latex, "fig_exponent.pdf"), bbox_inches="tight")
plt.close(fig)

dev = np.abs(np.array(slope_max + slope_min) - np.tile(1 - 2 * alphas, 2))
print(f"fig_sn.pdf, fig_exponent.pdf escritas; max |pendiente - (1-2 alpha)| = {dev.max():.3f}")
