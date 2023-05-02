# shifts_and_couplings.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Tue 02 May 2023 22:12:09 BST

from pathlib import Path
import pickle
import re

RESULT_DIR = Path("~/Documents/DPhil/results/cupid").expanduser()

with open(RESULT_DIR / "sucrose/shifts.pkl", "rb") as fh:
    shifts = pickle.load(fh)
with open(RESULT_DIR / "sucrose/couplings.pkl", "rb") as fh:
    couplings = pickle.load(fh)

table = """
\\begin{longtable}[h!]{c c c c}
\\caption[
The isotropic chemical shifts and scalar couplings used in simulated datasets presented in this work.
]{
The isotropic chemical shifts ($\\delta$), corresponding rotating frame
frequencies ($\\omega_{\\text{rot}}$), and scalar couplings ($J$) used in simulated
datasets presented in this work.
}
\\label{tab:shifts_and_couplings}\\\\
\\hline
Spin & $\\delta$ (\\unit{\\partspermillion}) & $\\omega_{\\text{rot}}$ (\\unit{\\hertz}) & $J$ (\\unit{\\hertz}) \\\\
\\hline
\\endfirsthead
\\hline
Spin & $\\delta$ (\\unit{\\partspermillion}) & $\\omega_{\\text{rot}}$ (\\unit{\\hertz}) & $J$ (\\unit{\\hertz}) \\\\
\\hline
\\endhead
\\hline
\\endlastfoot
\\multicolumn{4}{r}{Continues on next page...}\\\\
\\hline
\\endfoot
<FOUR-MP>
<SUCROSE>
\\hline
\\end{longtable}
"""

per_spin_couplings = {}
for coupling in couplings:
    spin1 = coupling[0]
    if spin1 in per_spin_couplings:
        per_spin_couplings[spin1].append(coupling[1:])
    else:
        per_spin_couplings[spin1] = [coupling[1:]]
    spin2 = coupling[1]
    if spin2 in per_spin_couplings:
        per_spin_couplings[spin2].append(coupling[::2])
    else:
        per_spin_couplings[spin2] = [coupling[::2]]

sucrose_data = (
    "\\hline\n"
    "\\multicolumn{4}{c}{\\textbf{Sucrose}}\\\\\n"
    "\\hline\n"
)
for i, shift in enumerate(shifts, start=1):
    sucrose_data += f"\\textbf{{{chr(64 + i)}}} & {shift:.3f} & {shift*300:.1f} & "
    if i in per_spin_couplings:
        sucrose_data += ", ".join(
            [f"\\textbf{{{chr(64 + spin2)}}}: ${freq:.3f}$"
             for spin2, freq in per_spin_couplings[i]]
        )
    else:
        sucrose_data += "--"
    sucrose_data += " \\\\\n"
table = table.replace("<SUCROSE>", sucrose_data)

four_mp_data = ""

with open(RESULT_DIR / "four_multiplets/shifts_and_couplings.txt", "r") as fh:
    txt = fh.readlines()

shift_regex = re.compile(r"^\d: (-?\d\.\d+)")
coupling_regex = re.compile(r"^\d: (.*?)\n")
for i in range(1, 6):
    four_mp_data += f"\\hline\n\\multicolumn{{4}}{{c}}{{\\textbf{{Four Multiplets, Run {i}}}}}\\\\\n\\hline\n"
    for _ in range(5):
        txt.pop(0)
    for _ in range(2):
        txt.pop(4)
    for spin in range(4):
        shift_line = txt[spin]
        shift_ppm = float(re.search(shift_regex, shift_line).group(1))
        shift_hz = float(shift_ppm) * 500.
        coupling_line = txt[spin + 4]
        couplings = ", ".join(
            [
                f"\\textbf{{{chr(69 + i)}:}} {float(x):.3f}" for i, x in
                enumerate(re.search(coupling_regex, coupling_line).group(1).split(", "))
            ]
        )
        four_mp_data += f"\\textbf{{{chr(65 + spin)}}} & \\num{{{shift_ppm:.2e}}} & {shift_hz:.2f} & {couplings} \\\\\n"
    for _ in range(9):
        txt.pop(0)

table = table.replace("<FOUR-MP>", four_mp_data)

with open("text/generated/shifts_and_couplings.tex", "w") as fh:
    fh.write(table)
