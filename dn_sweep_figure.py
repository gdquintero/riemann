"""Figura del barrido en alpha (dn_sweep.npz)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

d = np.load("dn_sweep.npz")
alphas, ns, DN2, floor = d["alphas"], d["ns"], d["DN2"], d["floor"]
gamma = np.array([-np.polyfit(np.log(ns[ns >= 160]),
                              np.log(DN2[i][ns >= 160]), 1)[0]
                  for i in range(len(alphas))])

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(9.4, 3.8))

# (a) familia d_n(alpha)^2 vs n
sel = [0.40, 0.47, 0.55, 0.65, 0.80]
for a in sel:
    i = int(np.argmin(np.abs(alphas - a)))
    ax1.loglog(ns, DN2[i], "o-", ms=4, label=fr"$\alpha={alphas[i]:.2f}$")
ax1.set_xlabel(r"$n$")
ax1.set_ylabel(r"$d_n(\alpha)^2$")
ax1.set_title(r"(a) optimal error vs. $n$")
ax1.grid(True, which="both", ls=":", alpha=0.5)
ax1.legend(frameon=False, fontsize=8)

# (b) exponente aparente gamma(alpha)
ax2.plot(alphas, gamma, "o-", color="C0", label=r"apparent $\gamma(\alpha)$")
# ajuste lineal en la porcion baja (alpha<=0.6) y su cruce por 0
m = alphas <= 0.6
sl, it = np.polyfit(alphas[m], gamma[m], 1)
xc = -it / sl
xf = np.linspace(xc, 0.7, 50)
ax2.plot(xf, sl * xf + it, "--", color="C2", alpha=0.7,
         label=r"linear fit ($\gamma=0$ at $\alpha\approx$" + f"{xc:.2f})")
ax2.axvline(0.5, ls=":", color="k", label=r"$\alpha=1/2$")
ax2.set_xlabel(r"$\alpha$")
ax2.set_ylabel(r"apparent exponent $\gamma$")
ax2.set_title(r"(b) $d_n^2\sim C\,n^{-\gamma}$ over $160\leq n\leq2560$")
ax2.grid(True, ls=":", alpha=0.5)
ax2.legend(frameon=False, fontsize=8, loc="upper left")

fig.tight_layout()
fig.savefig("latex/fig_sweep.pdf")
fig.savefig("latex/fig_sweep.png", dpi=160)
print(f"gamma(0.5)={gamma[np.argmin(np.abs(alphas-0.5))]:.3f}  "
      f"linear-fit zero at alpha={xc:.3f}")
print("Guardado: latex/fig_sweep.pdf/.png")
