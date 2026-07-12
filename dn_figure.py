"""
Figura para la Sugerencia 1 a partir de dn_results.npz (generado por dn_run.py).
Panel (a): d_n^2 (optimo) vs s_n (coef. de Mobius), escala log-log.
Panel (b): d_n^2 * log n vs n, con la constante conjeturada sum_rho 1/|rho|^2.
Guarda latex/fig_dn.pdf
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load("dn_results.npz")
ns, dn2, sn = d["ns"], d["dn2"], d["sn"]
C, D, C_ref = float(d["C"]), float(d["D"]), float(d["C_ref"])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.2, 3.7))

# ---- Panel (a) ----
ax1.loglog(ns, sn, "s", ms=5, mfc="none", color="C3",
           label=r"$s_n$  (Möbius $\mu(k)$)")
ax1.loglog(ns, dn2, "o", ms=5, color="C0",
           label=r"$d_n^2$  (optimal)")
ax1.set_xlabel(r"$n$")
ax1.set_ylabel(r"residual $\|\cdot\|_{1/2}^2$")
ax1.set_title(r"(a) optimal distance vs. Möbius")
ax1.grid(True, which="both", ls=":", alpha=0.5)
ax1.legend(frameon=False, loc="lower left")

# ---- Panel (b) ----
ax2.semilogx(ns, dn2 * np.log(ns), "o-", ms=5, color="C0",
             label=r"$d_n^2\,\log n$")
ax2.axhline(C_ref, ls="--", color="k",
            label=r"$\sum_\rho|\rho|^{-2}\approx%.4f$" % C_ref)
xf = np.linspace(ns[0], ns[-1], 400)
ax2.plot(xf, C + D / np.log(xf), ls="-.", color="C2", alpha=0.8,
         label=r"fit $C/\log n + D/(\log n)^2$")
ax2.set_xlabel(r"$n$")
ax2.set_ylabel(r"$d_n^2\,\log n$")
ax2.set_title(r"(b) plateau toward $\sum_\rho|\rho|^{-2}$")
ax2.set_ylim(0.035, 0.055)
ax2.grid(True, which="both", ls=":", alpha=0.5)
ax2.legend(frameon=False, loc="upper right")

fig.tight_layout()
fig.savefig("latex/fig_dn.pdf")
fig.savefig("latex/fig_dn.png", dpi=160)
print("Guardado: latex/fig_dn.pdf y latex/fig_dn.png")
print(f"C={C:.5f}  D={D:.5f}  C_ref={C_ref}")
