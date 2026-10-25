"""Tablas 2 y 3 del paper (s_n(alpha) en n = 1000, 10000, ..., 100000) con los resultados
de riemann_traj (traj/output_0XX.txt).

Escribe tabla_traj.tex: para cada grupo de alphas, la tabla Fortran con el mismo formato
que la del paper y, debajo, la original de Python tal como aparece en ../latex/riemann.tex.
"""
import os
import re
import numpy as np

here = os.path.dirname(os.path.abspath(__file__))
rows = [1000] + list(range(10000, 100001, 10000))
groups = [("tab:alpha05", list(range(50, 60)), "0.50 \\leq \\alpha \\leq 0.59"),
          ("tab:alpha09", list(range(90, 100)), "0.90 \\leq \\alpha \\leq 0.99")]

def paper_table(label):
    tex = open(os.path.join(here, "..", "latex", "riemann.tex")).read()
    i = tex.index(r"\label{%s}" % label)
    blk = tex[tex.rindex(r"\begin{table}", 0, i):i]
    vals = {}
    for line in blk.splitlines():
        m = re.match(r"\s*(\d+)\s*&(.*)\\\\", line)
        if m:
            vals[int(m.group(1))] = [float(a) * 10.0**int(b) for a, b in
                                     re.findall(r"\$([\d.]+)\\mathrm\{e\}\{([-+]\d+)\}\$", m.group(2))]
    return vals

def sci(x):
    m, e = f"{x:.3e}".split("e")
    return f"${m}\\mathrm{{e}}{{{int(e)}}}$"

def table(cols, cell, title, caption, label):
    head = [r"\begin{table}[ht!]", r"\centering",
            r"\resizebox{\textwidth}{!}{\begin{tabular}{|c|" + "c|" * len(cols) + "}",
            r"\hline",
            r" \multirow{2}{*}{$n$} & \multicolumn{" + str(len(cols)) + r"}{c|}{" + title + r"}\\",
            r" \cline{2-" + str(len(cols) + 1) + "}",
            " & " + " & ".join(f"$\\alpha = {a / 100:.2f}$" for a in cols) + r"\\",
            r" \hline\hline"]
    body = [f"{n:<6} & " + " & ".join(cell(n, a) for a in cols) + r" \\" for n in rows]
    tail = [r"\hline", r"\end{tabular}}", f"\\caption{{{caption}}}", f"\\label{{{label}}}",
            r"\end{table}"]
    return "\n".join(head + body + tail)

out = ["% Generado por tabla_traj.py"]
for label, cols, rng in groups:
    f90 = {a: dict(zip(*np.loadtxt(os.path.join(here, "traj", f"output_{a:03d}.txt"),
                                   usecols=(0, 2), unpack=True))) for a in cols}
    paper = paper_table(label)
    out.append(table(cols, lambda n, a: sci(f90[a][n]), "$s_n$ (Fortran)",
                     f"Values of the M\\\"obius error $s_n(\\alpha)$ computed with "
                     f"\\texttt{{riemann\\_traj}}, for ${rng}$.", label + "-f90"))
    out.append(table(cols, lambda n, a: sci(paper[n][cols.index(a)]), "$s_n$ (Python)",
                     f"Values of the M\\\"obius error $s_n(\\alpha)$ in the paper "
                     f"(\\texttt{{scipy.integrate.quad}}), for ${rng}$.", label + "-py"))

with open(os.path.join(here, "tabla_traj.tex"), "w") as f:
    f.write("\n\n".join(out) + "\n")
print("tabla_traj.tex escrita")
