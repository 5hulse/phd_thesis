# diff_params.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Mon 17 Jul 2023 18:49:19 BST

import nmrespy as ne
from utils import RESULT_DIR

resdir = RESULT_DIR / "diffusion"
experiments = {
    "Andrographolide": resdir / "andrographolide/estimator",
    "Glucose/valine/threonine": resdir / "gluc_val_thre/estimator",
}


table = """
\\null\\vfill
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
    Noteworthy experiment parameters for the diffusion datasets acquired.
]{
    Noteworthy experiment parameters for the diffusion datasets acquired.
    NS: Number of scans,
    DS: Number of dummy scans,
    PLW1: Hard pulse power (\\unit{\\watt}),
    P1: Duration of $\\nicefrac{\\pi}{2}$ pulse,
    D1: Duration of relaxation delay.
}
\\label{tab:onedim_params}
\\end{table}
\\vfill\\null
"""

table = table.replace("<TABULAR_CONFIG>", "".join((len(experiments) + 1) * ["c"]))
table = table.replace("<NAMES>", " & " + " & ".join(k for k in experiments) + "\\\\")
data_table = []
for name, path in experiments.items():
    estimator = ne.Estimator2DJ.from_pickle(path)
    acqus = estimator.bruker_params["acqus"]
    data_table.append(
        [
            f"{float(x):.5g}" for x in [
                acqus["BF1"],
                acqus["O1"],
                acqus["O1"] / acqus["SFO1"],
                acqus["SW_h"],
                acqus["SW"],
                acqus["TD"],
                acqus["NS"],
                acqus["DS"],
                acqus["PLW"][1],
                acqus["P"][1],
                acqus["D"][1],
            ]
        ]
    )

# transpose
data_table = list(map(list, zip(*data_table)))
names = [
    "$f_{\\text{bf}}$ (\\unit{\\mega\\hertz})",
    "$\\foff$ (\\unit{\\hertz})",
    "$\\foff$ (\\unit{\\partspermillion})",
    "$\\fsw$ (\\unit{\\hertz})",
    "$\\fsw$ (\\unit{\\partspermillion})",
    "$N$",
    "NS",
    "DS",
    "PLW1 (\\unit{\\watt})",
    "P1 (\\unit{\\micro\\second})",
    "D1 (\\unit{\\second})",
]
data = ""
for row, name in zip(data_table, names):
    row.insert(0, name)
    data += " & ".join(row) + "\\\\\n"
table = table.replace("<DATA>", data)

with open("text/generated/diff_params.tex", "w") as fh:
    fh.write(table)
