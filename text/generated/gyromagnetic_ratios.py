# gyromagnetic_ratios.py
# Simon Hulse
# simonhulse@protonmail.com
# Last Edited: Wed 31 Jan 2024 16:15:04 EST

import re
import numpy as np

info = [
    ("1H", 2.792847351, 0.5, "99.9885"),
    ("2H", 0.857438231, 1.0, "0.0115"),
    ("6Li", 0.822043, 1.0, "7.59"),
    ("7Li", 3.256407, 1.5, "92.41"),
    ("12C", None, 0.0, "98.93"),
    ("13C", 0.702369, 0.5, "1.07"),
    ("14N", 0.403573, 1.0, "99.636"),
    ("15N", -0.2830569, 0.5, "0.364"),
    ("16O", None, 0.0, "99.756"),
    ("17O", -1.893543, 2.5, "0.038"),
    ("19F", 2.628321, 0.5, "100"),
    ("31P", 1.130925, 0.5, "100"),
]
planck = 6.62607015e-34 / (2 * np.pi)
bohr = 5.0507837461e-27

def gamma(mu, I):
    return (mu * bohr) / (I * planck)


table = """\\begin{table}
    \\begin{center}
        \\begin{tabular}{ c c c c c }
            \\toprule
            Nucleus & $I$ & $\\mu$ ($\\mu_{\\text{N}}$) & $\\gamma (\\si{\\radian\\per\\tesla\\per\\second})$ & Relative Abundance (\\%) \\\\
            \\midrule
<DATA>
            \\bottomrule
        \\end{tabular}
    \\end{center}
    \\caption[
        Statistics related to a number of nuclei which are regularly-encountered in \\acs{NMR}.
    ]{
        \\correction{
            Statistics related to a number of nuclei which are
            regularly-encountered in \\acs{NMR}.
            Also listed are common nuclei which do not possess spin.
            The gyromagnetic ratios were determined by obtaining the relevant
            nuclear magnetic dipole moments $\\mu$ in units of nuclear magneton
            $\\mu_{\\text{N}}=\\qty{5.05078e-27}{\\joule\\per\\tesla}$,
            and applying the equation $\\gamma = \\nicefrac{\\mu
            \\mu_{\\text{N}}}{I \hbar}$~\\cite{Stone2019,Tiesinga2021}).
        }
    }
    \\label{tab:nuclei}
\\end{table}"""

data = ""
for (nuc, mu, I, abun) in info:
    symbol_match = re.match(r"(\d+)([A-Za-z]+)", nuc)
    mass = symbol_match.group(1)
    elem = symbol_match.group(2)
    nuc_lab = f"\\ch{{^{{{mass}}}{elem}}}"


    if mu is None:
        mu_lab = "--"
        gamma_lab = "--"
        I_lab = "$0$"
    else:
        I_num, I_denom = I.as_integer_ratio()
        if I_num == I_denom:
            I_lab = "$1$"
        else:
            I_lab = f"$\\nicefrac{{{I_num}}}{{{I_denom}}}$"
        mu_lab = "\\num{{{0:.3f}}}".format(mu)
        gamma_lab = "\\num{{{0:.4e}}}".format(gamma(mu, I)).replace("+0", "")

    abun_lab = f"\\num{{{abun}}}"

    data += "            {0} & {1} & {2} & {3} & {4} \\\\\n".format(nuc_lab, I_lab, mu_lab, gamma_lab, abun_lab)
data = data[:-1]
table = table.replace("<DATA>", data)

with open("text/generated/gyromagnetic_ratios.tex", "w") as fh:
    fh.write(table)
