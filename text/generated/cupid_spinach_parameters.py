# cupid_spinach_parameters.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Mon 26 Jun 2023 16:48:17 BST

from pathlib import Path
import nmrespy as ne

RESULT_DIR = Path("~/Documents/DPhil/results/cupid").expanduser()
experiments = {
    "Four Multiplets": RESULT_DIR / "four_multiplets/estimators/estimator_0",
    "Sucrose": RESULT_DIR / "sucrose/estimator",
}

table = """
\\begin{table}[h!]
\\centering
\\begin{tabular}{<TABULAR_CONFIG>}
\\hline
<NAMES>
\\hline
<DATA>
\\hline
\\end{tabular}
\\caption[
    Experiment parameters for \\ac{2DJ} simulations run using \\textsc{Spinach}.
]{
    Experiment parameters for \\ac{2DJ} simulations run using \\textsc{Spinach}.
}
\\label{tab:spinach-jres-params}
\\end{table}
"""

table = table.replace("<TABULAR_CONFIG>", "".join((len(experiments) + 1) * ["c"]))
table = table.replace("<NAMES>", "Parameter & " + " & ".join(k for k in experiments) + "\\\\")
jres_data = []
for name, path in experiments.items():
    estimator = ne.Estimator2DJ.from_pickle(path)
    data = [
        estimator.sfo[1],
        estimator.offset()[1],
        estimator.offset(unit="ppm")[1],
        estimator.sw()[0],
        estimator.sw()[1],
        estimator.sw(unit="ppm")[1],
        estimator.default_pts[0],
        estimator.default_pts[1],
    ]
    data = [f"{x:.4g}" if isinstance(x, (int, float)) else x for x in data]
    jres_data.append(data)

# transpose
jres_data = list(map(list, zip(*jres_data)))
names = [
    "$f_{\\text{bf}}^{(1)} (\\unit{\\mega\\hertz})$",
    "$\\fofftwo$ (\\unit{\\hertz})",
    "$\\fofftwo$ (\\unit{\\partspermillion})",
    "$\\fswone$ (\\unit{\\hertz})",
    "$\\fswtwo$ (\\unit{\\hertz})",
    "$\\fswtwo$ (\\unit{\\partspermillion})",
    "$\\None$",
    "$\\Ntwo$",
]
s = ""
for row, name in zip(jres_data, names):
    row.insert(0, name)
    s += " & ".join(row) + "\\\\\n"
s = s[:-1]
table = table.replace("<DATA>", s)

with open("text/generated/cupid_spinach_parameters.tex", "w") as fh:
    fh.write(table)
