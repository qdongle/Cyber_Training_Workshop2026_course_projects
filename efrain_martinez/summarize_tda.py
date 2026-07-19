from pathlib import Path
import csv
import re
import sys


materials = ["Cs2_Sn_Br6", "Cs2_Sn_I6", "Cs2_Ti_Br6", "Cs2_Ti_I6", "Cs2_Zr_Br6", "Cs2_Zr_I6"]
models = [("M0", "M0"), ("M1", "M1")]
columns = [
    "Material",
    "Model",
    "All_5_Converged",
    "State1_Energy_eV",
    "State1_Oscillator_Strength",
    "State1_Dominant_Transition",
    "State1_NTO_Weight",
    "State1_Hole_Component",
    "State1_Electron_Component",
    "Brightest_State",
    "Brightest_Energy_eV",
    "Brightest_Oscillator_Strength",
    "Brightest_Dominant_Transition",
    "Brightest_NTO_Weight",
    "Brightest_Hole_Component",
    "Brightest_Electron_Component",
]


def parse_states(text):
    states = {}
    before_nto = text.split("SELECTED NATURAL TRANSITION ORBITALS", 1)[0]
    pattern = re.compile(r"^State\s+([1-5])\s*$([\s\S]*?)(?=^State\s+[1-5]\s*$|\Z)", re.M)
    for match in pattern.finditer(before_nto):
        number = int(match.group(1))
        block = match.group(2)
        states[number] = {
            "converged": re.search(r"Converged:\s*(True|False)", block).group(1),
            "energy": re.findall(r"Excitation energy:\s*([-+0-9.eE]+)\s*eV", block)[-1],
            "osc": re.search(r"Oscillator strength:\s*([-+0-9.eE]+)", block).group(1),
            "transition": re.search(r"Dominant transitions:\s*\n\s*(.+)", block).group(1).strip(),
        }
    return states


def parse_ntos(text):
    ntos = {}
    nto_text = text.split("SELECTED NATURAL TRANSITION ORBITALS", 1)[1]
    pattern = re.compile(r"^State:\s*([1-5])\s*$([\s\S]*?)(?=^State:\s*[1-5]\s*$|\Z)", re.M)
    for match in pattern.finditer(nto_text):
        number = int(match.group(1))
        block = match.group(2)
        ntos[number] = {
            "weight": re.search(r"Dominant NTO pair weight:\s*([-+0-9.eE]+)", block).group(1),
            "hole": re.search(r"Dominant hole NTO composition:\s*\n\s*(.+)", block).group(1).strip(),
            "electron": re.search(r"Dominant electron NTO composition:\s*\n\s*(.+)", block).group(1).strip(),
        }
    return ntos


def parse_log(root, material, folder, label):
    text = (root / material / folder / "tda.log").read_text(errors="replace")
    states = parse_states(text)
    ntos = parse_ntos(text)
    brightest = max(range(1, 6), key=lambda number: float(states[number]["osc"]))
    return {
        "Material": material,
        "Model": label,
        "All_5_Converged": str(all(states[number]["converged"] == "True" for number in range(1, 6))),
        "State1_Energy_eV": states[1]["energy"],
        "State1_Oscillator_Strength": states[1]["osc"],
        "State1_Dominant_Transition": states[1]["transition"],
        "State1_NTO_Weight": ntos[1]["weight"],
        "State1_Hole_Component": ntos[1]["hole"],
        "State1_Electron_Component": ntos[1]["electron"],
        "Brightest_State": str(brightest),
        "Brightest_Energy_eV": states[brightest]["energy"],
        "Brightest_Oscillator_Strength": states[brightest]["osc"],
        "Brightest_Dominant_Transition": states[brightest]["transition"],
        "Brightest_NTO_Weight": ntos[brightest]["weight"],
        "Brightest_Hole_Component": ntos[brightest]["hole"],
        "Brightest_Electron_Component": ntos[brightest]["electron"],
    }


def write_text(path, rows):
    by_key = {(row["Material"], row["Model"]): row for row in rows}
    lines = [
        "VODP TDA-HF/CIS Summary",
        "",
        "Full 12-row table",
        "Material      Model  S1_eV     S1_f       S1_NTO    Bright  Bright_eV  Bright_f   Bright_NTO",
    ]
    for row in rows:
        lines.append(
            f"{row['Material']:<12} {row['Model']:<5} "
            f"{float(row['State1_Energy_eV']):>8.6f} "
            f"{float(row['State1_Oscillator_Strength']):>10.6f} "
            f"{float(row['State1_NTO_Weight']):>8.6f} "
            f"{row['Brightest_State']:>6} "
            f"{float(row['Brightest_Energy_eV']):>10.6f} "
            f"{float(row['Brightest_Oscillator_Strength']):>10.6f} "
            f"{float(row['Brightest_NTO_Weight']):>10.6f}"
        )

    lines += ["", "M0 versus M1 comparison", "Material      M0_S1_eV  M1_S1_eV  Shift_eV  M0_Bright_eV  M1_Bright_eV  Bright_Shift_eV"]
    for material in materials:
        m0 = by_key[(material, "M0")]
        m1 = by_key[(material, "M1")]
        lines.append(
            f"{material:<12} "
            f"{float(m0['State1_Energy_eV']):>9.6f} "
            f"{float(m1['State1_Energy_eV']):>9.6f} "
            f"{float(m1['State1_Energy_eV']) - float(m0['State1_Energy_eV']):>9.6f} "
            f"{float(m0['Brightest_Energy_eV']):>13.6f} "
            f"{float(m1['Brightest_Energy_eV']):>13.6f} "
            f"{float(m1['Brightest_Energy_eV']) - float(m0['Brightest_Energy_eV']):>16.6f}"
        )
    path.write_text("\n".join(lines) + "\n")


root = Path(sys.argv[1])
rows = [parse_log(root, material, folder, label) for material in materials for folder, label in models]

csv_path = root / "tda_summary.csv"
txt_path = root / "tda_summary.txt"
csv_path.unlink(missing_ok=True)
txt_path.unlink(missing_ok=True)

with csv_path.open("w", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=columns)
    writer.writeheader()
    writer.writerows(rows)

write_text(txt_path, rows)
print("Wrote tda_summary.csv and tda_summary.txt")
