import numpy as np
import pyspawn
import pyspawn.general
import pyspawn.process_geometry as pg
import pyspawn.dictionaries as dicts

# Read geometry
natoms, atoms, _, comment = pg.process_geometry("geometry.xyz")

# Masses in atomic units
mass_dict = dicts.get_atomic_masses()
masses = np.asarray([mass_dict[atom] * 1822.0 for atom in atoms for _ in range(3)])

# Gaussian widths
# These are generic starting values; verify/adjust later.
widths_dict = {
    "C": 22.7,
    "H": 4.7,
    "O": 16.0,
}
widths = np.asarray([widths_dict[atom] for atom in atoms for _ in range(3)])

# Wigner sampling
wigner_temp = 0
seed = 69884

# Propagators
clas_prop = "vv"
qm_prop = "fulldiag"
qm_ham = "adiabatic"
potential = "molcas_cas"

# Time settings
t0 = 0.0
ts = 0.5
tfinal = 750.0

numdims = natoms * 3
numstates = 4

# OpenMolcas options for TBHP SA(4)-CASSCF(6e,4o)/cc-pVDZ
molcas_options = {
    "method": "casscf",
    "pt2": False,
    "imaginary": 0.0,
    "ipea": 0.0,
    
    "basis": "cc-pVDZ",
    "atoms": atoms,
    "charge": 0,
    "spinmult": 1,
    "nactel": 6,
    "actorb": 4,
    "inactive": 22,

    "cassinglets": numstates,
    "castargetmult": 1,

    "cas_energy_states": [0, 1, 1],
    "cas_energy_mults": [1, 1, 1],

    "python3": "/usr/bin/python3",
    #"python3": "/cvmfs/soft.ccr.buffalo.edu/versions/2023.01/compat/usr/bin/python",
    "project": "TBHP",
}

traj_params = {
    "time": t0,
    "timestep": ts,
    "maxtime": tfinal,
    "istate": 1,
    "spawnthresh": (0.5 * np.pi) / ts / 20.0,
    "widths": widths,
    "atoms": molcas_options["atoms"],
    "masses": masses,
    "molcas_options": molcas_options,
}

sim_params = {
    "quantum_time": traj_params["time"],
    "timestep": traj_params["timestep"],
    "max_quantum_time": traj_params["maxtime"],
    "qm_amplitudes": np.ones(1, dtype=np.complex128),
    "qm_energy_shift": 0.0,
}

# Import propagation methods
exec("pyspawn.import_methods.into_simulation(pyspawn.qm_integrator." + qm_prop + ")")
exec("pyspawn.import_methods.into_simulation(pyspawn.qm_hamiltonian." + qm_ham + ")")
exec("pyspawn.import_methods.into_traj(pyspawn.potential." + potential + ")")
exec("pyspawn.import_methods.into_traj(pyspawn.classical_integrator." + clas_prop + ")")

pyspawn.general.check_files()

# Initial trajectory
traj1 = pyspawn.traj(numdims, numstates)
traj1.set_numstates(numstates)
traj1.set_numdims(numdims)
traj1.set_parameters(traj_params)

traj1.initial_wigner(seed)

# Simulation
sim = pyspawn.simulation()
sim.add_traj(traj1)
sim.set_parameters(sim_params)

# Optional SSAIMS
sim.enable_ssaims(
    epsilon=1e-10,
    ss_seed=527516,
    suspend_during_spawn=True,
    spawn_delay_steps=10,
    min_tbf_to_start=2,
    verbose=True,
)

sim.propagate()
