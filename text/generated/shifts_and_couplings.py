# shifts_and_couplings.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Mon 26 Jun 2023 19:09:41 BST

from pathlib import Path
import pickle
import re

RESULT_DIR = Path("~/Documents/DPhil/results/cupid").expanduser()


def pickle_load(path):
    with open(path, "rb") as fh:
        return pickle.load(fh)


def get_coupling_info(couplings):
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
    return per_spin_couplings


def make_tabular(shifts, couplings, t1s=None, t2s=None):
    if t1s is None:
        t1s = len(shifts) * [None]
    if t2s is None:
        t2s = len(shifts) * [None]
    tabular = ""
    for i, (shift, t1, t2) in enumerate(zip(shifts, t1s, t2s), start=1):
        tabular += f"\\textbf{{{chr(64 + i)}}} & {shift:.3f} & {shift*300:.1f} & "
        if i in couplings:
            tabular += ", ".join(
                [f"\\textbf{{{chr(64 + spin2)}}}: ${freq:.3f}$"
                 for spin2, freq in couplings[i]]
            )
        else:
            tabular += "--"

        for t in (t1, t2):
            if t is None:
                tabular += "& --"
            else:
                tabular += f"& {t:.2f}"

        tabular += " \\\\\n"
    return tabular


table = """
\\begin{longtable}[h!]{c c c c c c}
\\caption[
The isotropic chemical shifts, scalar couplings and relaxation times
associated with spin systems used in \\textsc{Spinach} simulations.
]{
The isotropic chemical shifts ($\\delta$), corresponding rotating frame
frequencies ($\\omega_{\\text{rot}}$), scalar couplings ($J$), and
relaxation times ($T_1$, $T_2$, if applicable)
associated with spin systems used in \\textsc{Spinach} simulations.
}
\\label{tab:shifts_and_couplings}\\\\
\\hline
Spin & $\\delta$ (\\unit{\\partspermillion}) & $\\omega_{\\text{rot}}$ (\\unit{\\hertz}) & $J$ (\\unit{\\hertz}) & $T_1$ (\\unit{\\second}) & $T_2$ (\\unit{\\second})\\\\
\\hline
\\endfirsthead
\\hline
Spin & $\\delta$ (\\unit{\\partspermillion}) & $\\omega_{\\text{rot}}$ (\\unit{\\hertz}) & $J$ (\\unit{\\hertz}) & $T_1$ (\\unit{\\second}) & $T_2$ (\\unit{\\second})\\\\
\\hline
\\endhead
\\hline
\\endlastfoot
\\multicolumn{6}{r}{Continues on next page...}\\\\
\\hline
\\endfoot
<FOUR-MP>
<SUCROSE>
<FIVE-MP>
<STRYCHININE>
\\hline
\\end{longtable}
"""

# SUCROSE CUPID
sucrose_shifts = pickle_load(RESULT_DIR/ "cupid/sucrose/shifts.pkl")
sucrose_couplings = get_coupling_info(pickle_load(RESULT_DIR/ "cupid/sucrose/couplings.pkl"))
sucrose_tabular = (
    "\\hline\n"
    "\\multicolumn{6}{c}{\\textbf{Sucrose}}\\\\\n"
    "\\hline\n"
)
sucrose_tabular += make_tabular(sucrose_shifts, sucrose_couplings)
table = table.replace("<SUCROSE>", sucrose_tabular)

# STRYCHININE INVREC
strychinine_shifts = pickle_load(RESULT_DIR / "invrec/strychinine/shifts.pkl")
strychinine_couplings = get_coupling_info(pickle_load(RESULT_DIR / "invrec/strychinine/couplings.pkl"))
strychinine_t1s = pickle_load(RESULT_DIR / "invrec/strychinine/t1s.pkl")
strychinine_t2s = pickle_load(RESULT_DIR / "invrec/strychinine/t2s.pkl")
strychinine_tabular = (
    "\\hline\n"
    "\\multicolumn{6}{c}{\\textbf{Strychinine}}\\\\\n"
    "\\hline\n"
)
strychinine_tabular += make_tabular(
    strychinine_shifts,
    strychinine_couplings,
    strychinine_t1s,
    strychinine_t2s,
)
table = table.replace("<STRYCHININE>", strychinine_tabular)


def get_info_from_file(file, title):
    with open(file, "r") as fh:
        txt = fh.readlines()
    tabular = ""
    shift_t1_t2_regex = re.compile(r"^\d: (-?\d\.\d+)")
    coupling_regex = re.compile(r"^\d: (.*?)\n")

    i = 1
    while True:
        if not txt:
            break
        tabular += f"\\hline\n\\multicolumn{{6}}{{c}}{{\\textbf{{{title}, Run {i}}}}}\\\\\n\\hline\n"
        for _ in range(5):
            txt.pop(0)
        if i == 1:
            nspins = 0
            while True:
                match = re.search(shift_t1_t2_regex, txt[nspins])
                if match is None:
                    break
                else:
                    nspins += 1
        for _ in range(2):
            txt.pop(nspins)

        rows = []
        for s in range(nspins):
            shift_line = txt[s]
            shift_ppm = float(re.search(shift_t1_t2_regex, shift_line).group(1))
            shift_hz = float(shift_ppm) * 500.
            coupling_line = txt[s + nspins]
            couplings = ", ".join(
                [
                    f"\\textbf{{{chr(65 + nspins + i)}:}} {float(x):.3f}" for i, x in
                    enumerate(re.search(coupling_regex, coupling_line).group(1).split(", "))
                ]
            )
            rows.append(f"\\textbf{{{chr(65 + s)}}} & \\num{{{shift_ppm:.2e}}} & {shift_hz:.2f} & {couplings} & <T1> & <T2>\\\\\n")
        for _ in range(2 * nspins):
            txt.pop(0)
        print(txt)
        if txt[0] == "\n":
            for s in range(nspins):
                rows[s] = rows[s].replace("<T1>", "--")
                rows[s] = rows[s].replace("<T2>", "--")
        else:
            for _ in range(2):
                txt.pop(0)
            for s in range(nspins):
                t1_line = txt.pop(0)
                t1 = float(re.search(shift_t1_t2_regex, t1_line).group(1))
                rows[s] = rows[s].replace("<T1>", f"{t1:.3f}")
                rows[s] = rows[s].replace("<T2>", "--")

        tabular += "\n".join(rows)
        txt.pop(0)
        print(txt)
        i += 1

    return tabular

# FOUR MULTIPLETS CUPID
four_mp_tabular = get_info_from_file(
    RESULT_DIR / "cupid/four_multiplets/shifts_and_couplings.txt",
    "Four Multiplets",
)
table = table.replace("<FOUR-MP>", four_mp_tabular)

# FIVE MULTIPLETS INVREC
five_mp_tabular = get_info_from_file(
    RESULT_DIR / "invrec/five_multiplets/shifts_and_couplings.txt",
    "Five Multiplets",
)
table = table.replace("<FIVE-MP>", five_mp_tabular)

with open("text/generated/shifts_and_couplings.tex", "w") as fh:
    fh.write(table)
