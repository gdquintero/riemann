"""Maximo y minimo de s_n(alpha) en cada decada de n, para todos los alphas de una corrida
de riemann_traj con muestreo logaritmico (por defecto traj_1e6/).

Escribe tabla_env.tex: una fila por alpha y, para cada decada [10^k, 10^(k+1)], las
columnas max s_n y min s_n.  Con paso p (en centesimas) solo se usan los alphas
multiplos de p/100 y el archivo es tabla_env_p.tex.

Uso:  python tabla_env.py [directorio] [paso]
"""
import glob
import os
import sys
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
folder = os.path.join(here, sys.argv[1] if len(sys.argv) > 1 else "traj_1e6")
step = int(sys.argv[2]) if len(sys.argv) > 2 else 1

data = []
for f in sorted(glob.glob(os.path.join(folder, "output_*.txt"))):
    D = np.loadtxt(f)
    if round(100 * D[0, 1]) % step == 0:
        data.append((D[0, 1], D[:, 0], D[:, 2]))
data.sort(key=lambda d: d[0])

n = data[0][1]
k0, k1 = int(round(np.log10(n.min()))), int(round(np.log10(n.max())))
decades = [(10**k, 10**(k + 1)) for k in range(k0, k1)]

def sci(x):
    m, e = f"{x:.3e}".split("e")
    return f"${m}\\mathrm{{e}}{{{int(e)}}}$"

def inside(n, lo, hi):
    # [lo, hi) except for the last decade, which includes nmax
    return (n >= lo) & ((n < hi) if hi < n.max() else (n <= hi))

lines = ["% Generado por tabla_env.py", r"\begin{table}[htb!]", r"\centering", r"\footnotesize",
         r"\begin{tabular}{|c|" + "cc|" * len(decades) + "}", r"\hline",
         r"\multirow{2}{*}{$\alpha$} & " + " & ".join(
             f"\\multicolumn{{2}}{{c|}}{{$10^{{{k0 + i}}} \\leq n "
             f"{'<' if i < len(decades) - 1 else chr(92) + 'leq'} 10^{{{k0 + i + 1}}}$}}"
             for i in range(len(decades))) + r" \\",
         r"\cline{2-" + str(2 * len(decades) + 1) + "}",
         " & " + " & ".join([r"$\max s_n$ & $\min s_n$"] * len(decades)) + r" \\",
         r"\hline\hline"]
for a, n, s in data:
    cells = [f"{a:.2f}"]
    for lo, hi in decades:
        m = inside(n, lo, hi)
        cells += [sci(s[m].max()), sci(s[m].min())]
    lines.append(" & ".join(cells) + r" \\")
lines += [r"\hline", r"\end{tabular}",
          r"\caption{Largest and smallest value of the M\"obius error $s_n(\alpha)$ over each "
          r"decade of $n$, computed with \texttt{riemann\_traj} at " + str(len(n)) +
          r" log-spaced values of $n$.}",
          r"\label{tab:envelope}", r"\end{table}"]

name = "tabla_env.tex" if step == 1 else f"tabla_env_{step}.tex"
with open(os.path.join(here, name), "w") as f:
    f.write("\n".join(lines) + "\n")
print(f"{name} escrita ({len(data)} alphas, {len(decades)} decadas)")
