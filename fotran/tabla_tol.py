"""Tabla 1 del paper (n y tiempo para s_n(alpha) < eps) con los resultados de riemann_tol.

Escribe tabla_tol.tex: la tabla Fortran con el mismo formato que la del paper y, debajo,
la original de Python (../output_tol_*.txt) para comparar.  n es el indice real de S_n;
los scripts de Python reportaban n + 1.
"""
import os
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
alphas = [0.50, 0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95]
eps = ["1e{-}2", "1e{-}3", "1e{-}4"]

def sci(t):
    m, e = f"{t:.3e}".split("e")
    return f"${m}\\mathrm{{e}}{{{int(e):+d}}}$"

def row(data, a, fortran):
    cells = [f"{a:.2f}"]
    for D in data:
        i = int(np.argmin(np.abs(D[:, 0] - a)))
        if fortran:
            n = f"{int(D[i, 2])}" if D[i, 4] == 1 else f"$>{int(D[i, 2])}$"
        else:
            n = f"{int(D[i, 2])}"
        cells += [n, sci(D[i, 3])]
    return " & ".join(cells) + r" \\"

def table(data, fortran, caption, label):
    head = [r"\begin{table}[htb!]", r"\centering", r"\resizebox{0.5\textwidth}{!}{",
            r"\begin{tabular}{|c|cc|cc|cc|}", r"\hline",
            r"\multirow{2}{*}{$\alpha$} & " + " & ".join(
                f"\\multicolumn{{2}}{{c|}}{{$\\varepsilon = \\mathrm{{{e}}}$}}" for e in eps) + r" \\",
            r"\cline{2-7}", r" & $n$ & Time (s) & $n$ & Time (s) & $n$ & Time (s) \\",
            r"\hline\hline"]
    body = [row(data, a, fortran) for a in alphas]
    tail = [r"\hline", r"\end{tabular}", "}", f"\\caption{{{caption}}}", f"\\label{{{label}}}",
            r"\end{table}"]
    return "\n".join(head + body + tail)

f90 = [np.loadtxt(os.path.join(here, f"output_tol_{t}.txt")) for t in ["0.01", "0.001", "0.0001"]]
py = [np.loadtxt(os.path.join(here, "..", f"output_tol_{t}.txt")) for t in ["0.01", "0.005", "0.0005"]]

with open(os.path.join(here, "tabla_tol.tex"), "w") as f:
    f.write("% Generado por tabla_tol.py\n")
    f.write(table(f90, True,
                  r"Fortran (\texttt{riemann\_tol}): smallest $n$ with $s_n(\alpha) < \varepsilon$ and "
                  r"elapsed time; $>n$ means not reached for $n \leq 10^4$.", "tab:tolerance-f90"))
    f.write("\n\n")
    f.write(table(py, False,
                  r"Python (\texttt{scipy.integrate.quad}), as in the paper: entries at the "
                  r"$3600$~s cap were not attained; $n$ reported as $n+1$.", "tab:tolerance-py"))
    f.write("\n")
print("tabla_tol.tex escrita")
