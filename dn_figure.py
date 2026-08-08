"""
Figura para la Sugerencia 1 a partir de dn_results.npz (generado por dn_run.py).
Un solo panel: d_n^2 (optimo) vs s_n (coef. de Mobius), escala log-log.
El comportamiento de d_n^2 * log n frente a sum_rho 1/|rho|^2 va en la Tabla 4.
Guarda latex/fig_dn.pdf  (disenada a tamano final: width=0.70\\textwidth)
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, NullFormatter

# --- estilo: serif + matematica en Computer Modern, como el documento LaTeX ---
plt.rcParams.update({
    "font.family": "serif",
    "font.serif": ["DejaVu Serif"],     # cmr10 no trae "ö" (Möbius)
    "mathtext.fontset": "cm",           # la matematica sale en Computer Modern
    "axes.formatter.use_mathtext": True,
    "axes.unicode_minus": False,
    "font.size": 10,
    "axes.labelsize": 10.5,
    "legend.fontsize": 9,
    "xtick.labelsize": 9,
    "ytick.labelsize": 9,
    "axes.linewidth": 0.7,
    "xtick.direction": "out",
    "ytick.direction": "out",
    "xtick.major.width": 0.7,
    "ytick.major.width": 0.7,
    "xtick.minor.width": 0.5,
    "ytick.minor.width": 0.5,
    "savefig.bbox": "tight",
    "savefig.pad_inches": 0.02,
})

SURFACE = "#ffffff"
INK     = "#0b0b0b"
INK2    = "#52514e"
GRID    = "#e6e5e2"
C_OPT   = "#2a78d6"   # slot 1 (blue)   -> d_n^2 optimo
C_MOB   = "#eb6834"   # slot 2 (orange) -> s_n de Mobius

d = np.load("dn_results.npz")
ns, dn2, sn = d["ns"], d["dn2"], d["sn"]
C, D, C_ref = float(d["C"]), float(d["D"]), float(d["C_ref"])

fig, ax = plt.subplots(figsize=(4.4, 3.2))
fig.patch.set_facecolor(SURFACE)

ax.plot(ns, sn, "-s", lw=0.9, alpha=0.85, ms=4.5, mfc=SURFACE, mew=1.1,
        color=C_MOB, zorder=2, label=r"$s_n$  (Möbius $\mu(k)$)")
ax.plot(ns, dn2, "-o", lw=1.4, ms=4.5, color=C_OPT, mec=SURFACE, mew=0.8,
        zorder=3, label=r"$d_n^2$  (optimal)")

ax.set_xscale("log")
ax.set_yscale("log")
ax.set_xlim(8, 7000)
ax.set_ylim(2.2e-3, 2.6e-1)
ax.set_xlabel(r"$n$", color=INK)
ax.set_ylabel(r"residual $\|\,\cdot\,\|_{1/2}^{2}$", color=INK)
ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10) * 0.1,
                                      numticks=20))
ax.yaxis.set_minor_formatter(NullFormatter())

# ejes recesivos: sin marco superior/derecho, rejilla hairline solida
for side in ("top", "right"):
    ax.spines[side].set_visible(False)
for side in ("left", "bottom"):
    ax.spines[side].set_color(INK2)
ax.tick_params(colors=INK2, which="both")
ax.set_axisbelow(True)
ax.grid(True, which="major", color=GRID, lw=0.6, ls="-")
ax.grid(False, which="minor")

leg = ax.legend(loc="lower left", frameon=False, handlelength=1.8,
                borderaxespad=0.2, labelspacing=0.35, handletextpad=0.6)
for txt in leg.get_texts():
    txt.set_color(INK)

# etiquetas directas: la identidad no depende solo del color
ax.annotate(r"$s_n$", xy=(ns[-1], sn[-1]), xytext=(6, 2),
            textcoords="offset points", color=C_MOB, fontsize=9.5)
ax.annotate(r"$d_n^2$", xy=(ns[-1], dn2[-1]), xytext=(6, -3),
            textcoords="offset points", color=C_OPT, fontsize=9.5)

fig.tight_layout()
fig.savefig("latex/fig_dn.pdf")
fig.savefig("latex/fig_dn.png", dpi=220)
print("Guardado: latex/fig_dn.pdf y latex/fig_dn.png")
print(f"C={C:.5f}  D={D:.5f}  C_ref={C_ref}  (plateau -> Tabla 4)")
