# spinach_parameters.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 11 May 2023 20:22:13 BST

from pathlib import Path
import nmrespy as ne

RESULT_DIR = Path("~/Documents/DPhil/results").expanduser()
jres_experiments = {
    "Four Multiplets": RESULT_DIR / "cupid/four_multiplets/estimators/estimator_0",
    "Sucrose": RESULT_DIR / "cupid/sucrose/estimator",
}
invrec_experiments = {
    "Five Multiplets": RESULT_DIR / "invrec/five_multiplets/estimators/estimator_0",
}


jres_table = """
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

jres_table = jres_table.replace("<TABULAR_CONFIG>", "".join((len(jres_experiments) + 1) * ["c"]))
jres_table = jres_table.replace("<NAMES>", "Parameter & " + " & ".join(k for k in jres_experiments) + "\\\\")
jres_data = []
for name, path in jres_experiments.items():
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
jres_str = ""
for row, name in zip(jres_data, names):
    row.insert(0, name)
    jres_str += " & ".join(row) + "\\\\\n"
jres_str = jres_str[:-1]
jres_table = jres_table.replace("<DATA>", jres_str)



 # INVREC
invrec_table = """
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
    Experiment parameters for inversion recovery simulations run using \\textsc{Spinach}.
]{
    Experiment parameters for inversion recovery simulations run using \\textsc{Spinach}.
    $K$ specifies the number of increments run, and $\\tau_{\\text{max}}$
    specifies the largest delay time used. Delays were generated with linear spacings,
    with the first delay always being \qty{0}{\second}, such that the
    n\\textsuperscript{th} delay was $\\nicefrac{(n-1)\\tau_{\\text{max}}}{K -
    1}$.
}
\\label{tab:spinach-invrec-params}
\\end{table}
"""

invrec_table = invrec_table.replace("<TABULAR_CONFIG>", "".join((len(invrec_experiments) + 1) * ["c"]))
invrec_table = invrec_table.replace("<NAMES>", "Parameter &" + " & ".join(k for k in invrec_experiments) + "\\\\")
invrec_data = []
for name, path in invrec_experiments.items():
    estimator = ne.Estimator2DJ.from_pickle(path)
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
    invrec_data.append(data)

# transpose
invrec_data = list(map(list, zip(*invrec_data)))
names = [
    "$f_{\\text{bf}}^{(1)} (\\unit{\\mega\\hertz})$",
    "$\\foffone$ (\\unit{\\hertz})",
    "$\\foffone$ (\\unit{\\partspermillion})",
    "$\\fswone$ (\\unit{\\hertz})",
    "$\\fswone$ (\\unit{\\partspermillion})",
    "$\\None$",
    "$K$",
    "$\\tau_{\\text{max}}$ (\\unit{\\second})",
]
invrec_str = ""
for row, name in zip(invrec_data, names):
    row.insert(0, name)
    invrec_str += " & ".join(row) + "\\\\\n"
invrec_str = invrec_str[:-1]
invrec_table = invrec_table.replace("<DATA>", invrec_str)

tables = f"{jres_table}\n{invrec_table}"

with open("text/generated/spinach_parameters.tex", "w") as fh:
    fh.write(tables)
