# shifts_and_couplings.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Fri 28 Apr 2023 19:02:04 BST

from pathlib import Path
import pickle
import re

RESULT_DIR = Path("~/Documents/DPhil/results/cupid").expanduser()

with open(RESULT_DIR / "sucrose/shifts.pkl", "rb") as fh:
    shifts = pickle.load(fh)
with open(RESULT_DIR / "sucrose/couplings.pkl", "rb") as fh:
    couplings = pickle.load(fh)

table = """
\\begin{longtable}[h!]{c c c}
\\caption[
    The isotropic chemical shifts and scalar couplings used in simulated datasets.
]{
The isotropic chemical shifts and scalar couplings used in the simulated examples presented in Section \\ref{subsec:simulated_results}.
}
\\label{tab:shifts_and_couplings}\\\\
\\hline
\\textbf{Spin} & \\textbf{Shift} (\\unit{\\partspermillion}/\\unit{\\hertz}) & \\textbf{Couplings} (\\unit{\\hertz}) \\\\
\\hline
\\endfirsthead
\\hline
\\textbf{Spin} & \\textbf{Shift} (\\unit{\\partspermillion}/\\unit{\\hertz}) & \\textbf{Couplings} (\\unit{\\hertz}) \\\\
\\hline
\\endhead
\\hline
\\endlastfoot
\\multicolumn{3}{r}{Continues on next page...}\\\\
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
    "\\multicolumn{3}{c}{\\textbf{Sucrose}}\\\\\n"
    "\\hline\n"
)
for i, shift in enumerate(shifts, start=1):
    sucrose_data += f"\\textbf{{{chr(64 + i)}}} & {shift:.4g}/{shift*300:.4g} & "
    if i in per_spin_couplings:
        sucrose_data += ", ".join(
            [f"\\textbf{{{chr(64 + spin2)}}}: ${freq:.2f}$"
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
    four_mp_data += f"\\hline\n\\multicolumn{{3}}{{c}}{{\\textbf{{Four Multiplets, Run {i}}}}}\\\\\n\\hline\n"
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
                f"\\textbf{{{chr(69 + i)}:}} {float(x):.4g}" for i, x in
                enumerate(re.search(coupling_regex, coupling_line).group(1).split(", "))
            ]
        )
        four_mp_data += f"\\textbf{{{chr(65 + spin)}}} & {shift_ppm:.4g}/{shift_hz:.4g} & {couplings} \\\\\n"
    for _ in range(9):
        txt.pop(0)

table = table.replace("<FOUR-MP>", four_mp_data)

with open("text/generated/shifts_and_couplings.tex", "w") as fh:
    fh.write(table)
