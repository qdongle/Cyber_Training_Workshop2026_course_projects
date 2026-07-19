import csv
import sys
import numpy as np
from pathlib import Path
from pyscf import scf, qmmm, tdscf
from pyscf.scf import chkfile

models = {"M0": ("M0", False), "M1": ("M1", True)}
hartree_to_ev = 27.211386245988

def mo_label(index, homo, lumo):
    if index <= homo:
        gap = homo - index
        return "HOMO" if gap == 0 else f"HOMO-{gap}"
    gap = index - lumo
    return "LUMO" if gap == 0 else f"LUMO+{gap}"

if len(sys.argv) != 2:
    raise SystemExit("Usage: python run_tda.py MODEL_FOLDER")

model_folder = Path(sys.argv[1]).resolve()
material = model_folder.parent.name
model_name = model_folder.name
label, embedded = models[model_name]

tda_path = model_folder / "tda.log"
chk_path = model_folder / "scf.chk"
pc_path = model_folder / "point_charges.csv"
tda_path.unlink(missing_ok=True)

mol, data = chkfile.load_scf(str(chk_path))
mo_coeff = np.asarray(data["mo_coeff"])
mo_occ = np.asarray(data["mo_occ"])
mo_energy = np.asarray(data["mo_energy"])

pc_coords = np.zeros((0, 3))
pc_values = np.zeros(0)
if embedded:
    pc_rows = list(csv.DictReader(pc_path.open()))
    pc_coords = np.array([[float(row["x_A"]), float(row["y_A"]), float(row["z_A"])] for row in pc_rows])
    pc_values = np.array([float(row["charge_e"]) for row in pc_rows])

mol.verbose = 0
mol.max_memory = 16000
mf = scf.RHF(mol).density_fit()
if embedded:
    mf = qmmm.mm_charge(mf, pc_coords, pc_values, unit="Angstrom")
mf.mo_coeff = mo_coeff
mf.mo_occ = mo_occ
mf.mo_energy = mo_energy
mf.e_tot = float(data["e_tot"])
mf.converged = True
mf.verbose = 0

with tda_path.open("w") as log:
    mol.stdout = log
    mf.stdout = log
    td = tdscf.TDA(mf)
    td.stdout = log
    td.verbose = 0
    td.singlet = True
    td.nstates = 5
    td.conv_tol = 1.0e-5
    td.max_cycle = 100

    active_indices = np.where(np.asarray(td.get_frozen_mask(), dtype=bool))[0]
    active_occ = mo_occ[active_indices]
    occupied_indices = active_indices[active_occ > 0]
    virtual_indices = active_indices[active_occ == 0]
    homo = int(occupied_indices[-1])
    lumo = int(virtual_indices[0])

    energies, xy = td.kernel()
    energies = np.asarray(energies, dtype=float)
    strengths = np.asarray(td.oscillator_strength(gauge="length"), dtype=float)
    brightest_state = int(np.argmax(strengths))
    selected_states = [0] if brightest_state == 0 else [0, brightest_state]

    log.write("VODP TDA-HF/CIS RESULTS\n")
    log.write(f"Material: {material}\n")
    log.write(f"Model: {label}\n")
    log.write(f"SCF energy: {mf.e_tot:.12f} Hartree\n")
    log.write(f"Point charges: {len(pc_values)}\n\n")

    for state in range(5):
        converged = td.converged[state] if isinstance(td.converged, (list, tuple, np.ndarray)) else td.converged
        x = np.asarray(td.xy[state][0] if isinstance(td.xy[state], (list, tuple)) else td.xy[state])
        weights = np.abs(x) ** 2
        weights /= np.sum(weights)

        log.write(f"State {state + 1}\n")
        log.write(f"Converged: {bool(converged)}\n")
        log.write(f"Excitation energy: {energies[state]:.8f} Hartree\n")
        log.write(f"Excitation energy: {energies[state] * hartree_to_ev:.6f} eV\n")
        log.write(f"Oscillator strength: {strengths[state]:.6f}\n\nDominant transitions:\n")
        for item in np.argsort(weights.ravel())[::-1][:3]:
            i, a = np.unravel_index(item, x.shape)
            i_abs = int(occupied_indices[i])
            a_abs = int(virtual_indices[a])
            log.write(f"    {mo_label(i_abs, homo, lumo)} (MO {i_abs + 1}) -> {mo_label(a_abs, homo, lumo)} (MO {a_abs + 1})\n")
            log.write(f"        amplitude = {float(x[i, a]):.10f}\n")
            log.write(f"        relative weight = {float(weights[i, a] * 100.0):.4f} %\n")
        log.write("\n")

    log.write("SELECTED NATURAL TRANSITION ORBITALS\n\n")
    for state in selected_states:
        nto_weights, _ = td.get_nto(state=state + 1, threshold=0.10, verbose=0)
        x = np.asarray(td.xy[state][0] if isinstance(td.xy[state], (list, tuple)) else td.xy[state])
        x = x / np.linalg.norm(x)
        u, _, vh = np.linalg.svd(x, full_matrices=False)

        log.write(f"State: {state + 1}\n")
        log.write(f"Reason selected: {'lowest energy' if state == 0 else 'brightest'}\n")
        log.write(f"Excitation energy: {energies[state] * hartree_to_ev:.6f} eV\n")
        log.write(f"Oscillator strength: {strengths[state]:.6f}\n")
        log.write(f"Dominant NTO pair weight: {float(nto_weights[0]):.6f}\n")
        log.write(f"Remaining NTO-pair weight: {max(0.0, 1.0 - float(nto_weights[0])):.6f}\n\n")

        log.write("Dominant hole NTO composition:\n")
        hole_weights = np.abs(u[:, 0]) ** 2
        for item in np.argsort(hole_weights)[::-1][:3]:
            mo = int(occupied_indices[item])
            log.write(f"    {mo_label(mo, homo, lumo)} (MO {mo + 1}): coefficient = {float(u[item, 0]):.4f}, weight = {float(hole_weights[item] * 100.0):.2f} %\n")

        log.write("\nDominant electron NTO composition:\n")
        electron_weights = np.abs(vh[0, :]) ** 2
        for item in np.argsort(electron_weights)[::-1][:3]:
            mo = int(virtual_indices[item])
            log.write(f"    {mo_label(mo, homo, lumo)} (MO {mo + 1}): coefficient = {float(vh[0, item]):.4f}, weight = {float(electron_weights[item] * 100.0):.2f} %\n")
        log.write("\n")

print(f"{material} {label}: 5 singlet TDA states converged; E1={energies[0] * hartree_to_ev:.4f} eV, f1={strengths[0]:.6f}")
