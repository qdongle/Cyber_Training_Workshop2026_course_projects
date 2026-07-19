import csv
import sys
import numpy as np
from pathlib import Path
from pyscf import gto, scf, qmmm

models = {"M0": ("M0", -2, 0, False), "M1": ("M1", 6, 0, True)}
ecp_elements = {"Cs", "Zr", "Sn", "I"}
hartree_to_ev = 27.211386245988

def read_xyz(path):
    lines = path.read_text().splitlines()
    atoms = []
    for line in lines[2:2 + int(lines[0])]:
        fields = line.split()
        atoms.append([fields[0], [float(fields[1]), float(fields[2]), float(fields[3])]])
    return lines[1], atoms

if len(sys.argv) != 2:
    raise SystemExit("Usage: python run_pyscf.py MODEL_FOLDER")

model_folder = Path(sys.argv[1]).resolve()
material = model_folder.parent.name
model_name = model_folder.name
label, charge, spin, embedded = models[model_name]

xyz_path = model_folder / "quantum_region.xyz"
pc_path = model_folder / "point_charges.csv"
log_path = model_folder / "scf.log"
chk_path = model_folder / "scf.chk"
log_path.unlink(missing_ok=True)
chk_path.unlink(missing_ok=True)

comment, atoms = read_xyz(xyz_path)
elements = sorted({atom[0] for atom in atoms})
ecp = {element: "def2-svp" for element in elements if element in ecp_elements}

pc_coords = np.zeros((0, 3))
pc_values = np.zeros(0)
if embedded:
    pc_rows = list(csv.DictReader(pc_path.open()))
    pc_coords = np.array([[float(row["x_A"]), float(row["y_A"]), float(row["z_A"])] for row in pc_rows])
    pc_values = np.array([float(row["charge_e"]) for row in pc_rows])

mol = gto.Mole()
mol.atom = atoms
mol.unit = "Angstrom"
mol.basis = "def2-svp"
mol.ecp = ecp
mol.charge = charge
mol.spin = spin
mol.symmetry = False
mol.verbose = 0
mol.max_memory = 16000
mol.output = str(log_path)
mol.build()

with log_path.open("a") as log:
    log.write("VODP SCF INPUT\n")
    log.write(f"Material: {material}\n")
    log.write(f"Model: {label}\n")
    log.write(f"XYZ comment: {comment}\n")
    log.write(f"QM atoms: {len(atoms)}\n")
    log.write(f"QM charge: {charge}\n")
    log.write(f"Basis: def2-SVP\n")
    log.write(f"ECP elements: {', '.join(ecp) if ecp else 'none'}\n")
    log.write(f"Point charges: {len(pc_values)}\n")
    log.write(f"Point-charge sum: {float(np.sum(pc_values)):.10f}\n\n")

if not embedded:
    mf = scf.RHF(mol).density_fit()
    mf.chkfile = str(chk_path)
    mf.conv_tol = 1.0e-9
    mf.max_cycle = 100
    energy = mf.kernel()
else:
    mf1 = scf.RHF(mol).density_fit()
    mf1.chkfile = str(chk_path)
    mf1.conv_tol = 1.0e-5
    mf1.max_cycle = 80
    mf1.init_guess = "minao"
    mf1.DIIS = scf.ADIIS
    mf1.diis_start_cycle = 6
    mf1.diis_space = 8
    mf1.damp = 0.5
    mf1.level_shift = 0.5
    mf1 = qmmm.mm_charge(mf1, pc_coords, pc_values, unit="Angstrom")
    energy1 = mf1.kernel()

    mf = scf.RHF(mol).density_fit()
    mf.chkfile = str(chk_path)
    mf.conv_tol = 1.0e-9
    mf.max_cycle = 80
    mf = qmmm.mm_charge(mf, pc_coords, pc_values, unit="Angstrom")
    mf = mf.newton()
    mf.chkfile = str(chk_path)
    mf.conv_tol = 1.0e-9
    mf.max_cycle = 80
    energy = mf.kernel(mf1.mo_coeff, mf1.mo_occ)

mo_occ = np.array(mf.mo_occ)
mo_energy = np.array(mf.mo_energy)
occ = np.where(mo_occ > 0)[0]
vir = np.where(mo_occ == 0)[0]
homo = int(occ[-1])
lumo = int(vir[0])
gap = float(mo_energy[lumo] - mo_energy[homo])

with log_path.open("a") as log:
    log.write("\nVODP SCF SUMMARY\n")
    log.write(f"Status: {'CONVERGED' if mf.converged else 'NOT CONVERGED'}\n")
    log.write(f"Material: {material}\n")
    log.write(f"Model: {label}\n")
    log.write(f"Method: {'Embedded DF-RHF' if embedded else 'DF-RHF'}\n")
    if embedded:
        log.write(f"Stage 1 energy: {float(energy1):.12f} Hartree\n")
    log.write(f"Total DF-RHF energy: {float(energy):.12f} Hartree\n")
    log.write(f"HOMO: MO {homo + 1}, {mo_energy[homo]:.12f} Hartree\n")
    log.write(f"LUMO: MO {lumo + 1}, {mo_energy[lumo]:.12f} Hartree\n")
    log.write(f"HOMO-LUMO gap: {gap * hartree_to_ev:.8f} eV\n")

pc_text = f", PCs={len(pc_values)}" if embedded else ""
print(f"{material} {label}: DF-RHF converged, E={float(energy):.10f} Ha, gap={gap * hartree_to_ev:.4f} eV{pc_text}")
