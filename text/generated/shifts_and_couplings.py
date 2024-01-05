# shifts_and_couplings.py
# Simon Hulse
# simon.hulse@chem.ox.ac.uk
# Last Edited: Thu 04 Jan 2024 16:01:31 GMT

from pathlib import Path
import pickle
import re

RESULT_DIR = Path("~/Documents/DPhil/results").expanduser()


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


def make_tabular(shifts, couplings, sfo, t1s=None):
    if t1s is None:
        t1s = len(shifts) * [None]
    tabular = ""
    for i, (shift, t1) in enumerate(zip(shifts, t1s), start=1):
        tabular += f"({chr(64 + i)}) & {shift:.3f} & {shift*sfo:.1f} & "
        if i in couplings:
            tabular += ", ".join(
                [f"({chr(64 + spin2)}): ${freq:.3f}$"
                 for spin2, freq in couplings[i]]
            )
        else:
            tabular += "--"

        if t1 is None:
                tabular += "& --"
        else:
            tabular += f"& {t1:.2f}"

        tabular += " \\\\\n"
    return tabular


table = """
\\begin{longtable}[h!]{c c c c c}
\\caption[
The isotropic chemical shifts, scalar couplings and relaxation times
associated with spin systems used in \\textsc{Spinach} simulations.
]{
The isotropic chemical shifts ($\\delta$), corresponding rotating frame
frequencies ($\\omega_0$), and scalar couplings ($J$)
associated with spin systems used in \\textsc{Spinach} simulations.
For the ``Five multiplets'' spin systems, the associated $T_1$ times are provided too.
}
\\label{tab:shifts_and_couplings}\\\\
\\hline
Spin & $\\delta$ (\\unit{\\partspermillion}) & $\\omega_0$ (\\unit{\\hertz}) & $J$ (\\unit{\\hertz}) & $T_1$ (\\unit{\\second}) \\\\
\\hline
\\endfirsthead
\\hline
Spin & $\\delta$ (\\unit{\\partspermillion}) & $\\omega_0$ (\\unit{\\hertz}) & $J$ (\\unit{\\hertz}) & $T_1$ (\\unit{\\second}) \\\\
\\hline
\\endhead
\\hline
\\endlastfoot
\\multicolumn{5}{r}{Continues on next page...}\\\\
\\hline
\\endfoot
<FIVE-MP>
<SUCROSE>
<FOUR-MP>
<STRYCHININE>
\\hline
\\end{longtable}
"""

# SUCROSE CUPID
sucrose_shifts = pickle_load(RESULT_DIR/ "cupid/sucrose/shifts.pkl")
sucrose_couplings = get_coupling_info(pickle_load(RESULT_DIR/ "cupid/sucrose/couplings.pkl"))
sucrose_tabular = (
    "\\hline\n"
    "\\multicolumn{5}{c}{\\textbf{Sucrose}}\\\\\n"
    "\\hline\n"
)
sucrose_tabular += make_tabular(sucrose_shifts, sucrose_couplings, 300)
table = table.replace("\n<SUCROSE>", "")
# If suucrose used, uncomment this line
# table = table.replace("<SUCROSE>", sucrose_tabular)

# STRYCHININE CUPID
strychinine_shifts = pickle_load(RESULT_DIR / "invrec/strychinine/shifts.pkl")
strychinine_couplings = get_coupling_info(pickle_load(RESULT_DIR / "invrec/strychinine/couplings.pkl"))
strychinine_tabular = (
    "\\hline\n"
    "\\multicolumn{5}{c}{\\textbf{Strychinine}}\\\\\n"
    "\\hline\n"
)

# Re-order
order = [i[0] for i in reversed(sorted(enumerate(strychinine_shifts), key=lambda x: x[1]))]
strychinine_shifts = [strychinine_shifts[i] for i in order]
strychinine_couplings_new = {}
for i, j in enumerate(order):
    old_key = j + 1
    new_key = i + 1
    new_value = []
    for (old_spin, coupling) in strychinine_couplings[old_key]:
        new_spin = order.index(old_spin - 1) + 1
        new_value.append((new_spin, coupling))
    strychinine_couplings_new[new_key] = new_value

strychinine_tabular += make_tabular(
    strychinine_shifts,
    strychinine_couplings_new,
    500.,
    None,
    # strychinine_t1s,
    # strychinine_t2s,
)
table = table.replace("<STRYCHININE>", strychinine_tabular)


def get_info_from_file(file, title, shift=0.):
    with open(file, "r") as fh:
        txt = fh.readlines()
    tabular = ""
    shift_t1_t2_regex = re.compile(r"^\d: (-?\d\.\d+)")
    coupling_regex = re.compile(r"^\d: (.*?)\n")

    i = 1
    while True:
        if not txt:
            break
        tabular += f"\\hline\n\\multicolumn{{5}}{{c}}{{\\textbf{{{title}, Run {i}}}}}\\\\\n\\hline\n"
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
            shift_ppm = float(re.search(shift_t1_t2_regex, shift_line).group(1)) - shift
            shift_hz = float(shift_ppm) * 500.
            coupling_line = txt[s + nspins]
            couplings = ", ".join(
                [
                    f"({chr(65 + nspins + i)}): {float(x):.3f}" for i, x in
                    enumerate(re.search(coupling_regex, coupling_line).group(1).split(", "))
                ]
            )
            rows.append(f"({chr(65 + s)}) & \\num{{{shift_ppm:.2e}}} & {shift_hz:.2f} & {couplings} & <T1> \\\\\n")
            # rows.append(f"\\textbf{{{chr(65 + s)}}} & \\num{{{shift_ppm:.2e}}} & {shift_hz:.2f} & {couplings} & <T1> & <T2>\\\\\n")
        for _ in range(2 * nspins):
            txt.pop(0)
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
    shift=1.4,
)
table = table.replace("<FIVE-MP>", five_mp_tabular)

with open("text/generated/shifts_and_couplings.tex", "w") as fh:
    fh.write(table)
