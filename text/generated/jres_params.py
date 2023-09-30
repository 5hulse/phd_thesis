# jres_params.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 28 Sep 2023 20:08:12 BST

import nmrespy as ne
from utils import RESULT_DIR

resdir = RESULT_DIR / "cupid"
experiments = {
    "Quinine": resdir / "quinine/estimator",
    "Dexamethasone": resdir / "dexamethasone/estimator",
    "Camphor": resdir / "camphor/estimator",
    "Estradiol": resdir / "estradiol_low_snr/estimator",
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
    Noteworthy experiment parameters for the 2D J-Reolved and PSYCHE experiments used.
]{
    Noteworthy experiment parameters for the 2D J-Reolved and PSYCHE experiments used.
    NS: Number of scans,
    DS: Number of dummy scans,
    PLW1: Hard pulse power (\\unit{\\watt}),
    P1: Duration of $\\nicefrac{\\pi}{2}$ pulse,
    D1: Duration of relaxation delay.
}
\\label{tab:jres_params}
\\end{table}
"""

table = table.replace("<TABULAR_CONFIG>", "".join((len(experiments) + 1) * ["c"]))
table = table.replace("<NAMES>", " & " + " & ".join(k for k in experiments) + "\\\\")
data_table = []
for name, path in experiments.items():
    estimator = ne.Estimator2DJ.from_pickle(path)
    acqus = estimator.bruker_params["acqus"]
    acqu2s = estimator.bruker_params["acqu2s"]
    data_table.append(
        [
            f"{float(x):.5g}" for x in [
                acqus["BF1"],
                acqus["O1"],
                acqus["O1"] / acqus["SFO1"],
                acqu2s["SW_h"],
                acqus["SW_h"],
                acqus["SW"],
                acqu2s["TD"],
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
    "$\\fofftwo$ (\\unit{\\hertz})",
    "$\\fofftwo$ (\\unit{\\partspermillion})",
    "$\\fswone$ (\\unit{\\hertz})",
    "$\\fswtwo$ (\\unit{\\hertz})",
    "$\\fswtwo$ (\\unit{\\partspermillion})",
    "$\\None$",
    "$\\Ntwo$",
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

print(table)
with open("text/generated/jres_params.tex", "w") as fh:
    fh.write(table)
