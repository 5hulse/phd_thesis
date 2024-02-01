# invrec_spinach_parameters.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Wed 31 Jan 2024 15:25:08 EST

from pathlib import Path
import nmrespy as ne

RESULT_DIR = Path("~/Documents/DPhil/results/invrec").expanduser()
experiments = {
    "Five Multiplets": RESULT_DIR / "five_multiplets/estimators/estimator_0",
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
\\caption{
    The parameters used for the ``five multiplets'' inversion recovery
    simulations run using \\textsc{Spinach}.
}
\\label{tab:spinach-invrec-params}
\\end{table}
"""

table = table.replace("<TABULAR_CONFIG>", "".join((len(experiments) + 1) * ["c"]))
table = table.replace("<NAMES>", "Parameter & " + " & ".join(k for k in experiments) + "\\\\")
jres_data = []
for name, path in experiments.items():
    estimator = ne.EstimatorInvRec.from_pickle(path)
    data = [
        estimator.sfo[0],
        estimator.offset()[0],
        estimator.offset(unit="ppm")[0],
        estimator.sw()[0],
        estimator.sw(unit="ppm")[0],
        estimator.default_pts[0],
        len(estimator.increments),
        estimator.increments[-1],
    ]
    data = [f"{x:.5g}" if isinstance(x, (int, float)) else x for x in data]
    jres_data.append(data)

# transpose
jres_data = list(map(list, zip(*jres_data)))
names = [
    "$f_{\\text{bf}}^{(1)} (\\unit{\\mega\\hertz})$",
    "$\\foffone$ (\\unit{\\hertz})",
    "$\\foffone$ (\\unit{\\partspermillion})",
    "$\\fswone$ (\\unit{\\hertz})",
    "$\\fswone$ (\\unit{\\partspermillion})",
    "$\\None$",
]
s = ""
for row, name in zip(jres_data, names):
    row.insert(0, name)
    s += " & ".join(row) + "\\\\\n"
s = s[:-1]
table = table.replace("<DATA>", s)

with open("text/generated/invrec_spinach_parameters.tex", "w") as fh:
    fh.write(table)
