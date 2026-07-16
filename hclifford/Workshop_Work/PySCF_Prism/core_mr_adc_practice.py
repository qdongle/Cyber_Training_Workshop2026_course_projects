from pyscf import gto, scf, mcscf
import prism.interface
import prism.mr_adc
import numpy as np

#furan 
mol = gto.Mole()
mol.atom = """
O 1.07108 0.02372 0.29535
C 0.58811 0.11195 -0.95745
H 1.30631 0.21974 -1.73774
C -0.76328 0.04720 -0.95265
H -1.40422 0.09175 -1.80540
C -1.14074 -0.09336 0.42696
H -2.12555 -0.17687 0.83087
C 0.01847 -0.09994 1.12458
H 0.23850 -0.17769 2.16469
"""
mol.basis = "6-31g"
mol.verbose = 4
mol.charge = 0
mol.build()

mf = scf.RHF(mol).density_fit()
mf.kernel()

mc = mcscf.CASSCF(mf, 6, 6).density_fit()
emc = mc.mc1step()[0]

interface = prism.interface.PYSCF(mf, mc, backend = 'opt_einsum').density_fit()
mr_adc = prism.mr_adc.MRADC(interface)
mr_adc.method = "mr-adc(2)-x"
mr_adc.method_type = "cvs-ip"
mr_adc.ncvs = 3
mr_adc.nroots = 8

e_diff, intensity, x = mr_adc.kernel()

from prism.tools.spectrum import plot
plot(e_diff, intensity, broadening = 0.5, omega_min = 540, omega_max = 550, plot = True, x_label = "Energy, eV", y_label = "Intensity", title = "XPS spectrum", filename = "mr_adc")

n_states = 3
weights = np.ones(n_states)/n_states
mc = mcscf.CASSCF(mf, 6, 6).state_average_(weights).density_fit()
mc.fix_spin_(ss=0) 
emc = mc.mc1step()[0]

state1_ci = mc.ci[0]
state2_ci = mc.ci[1]
state3_ci = mc.ci[2]
mo = mc.mo_coeff.copy()

mc = mcscf.CASCI(mf, 6, 6).density_fit()
emc = mc.casci(mo_coeff=mo, ci0=state1_ci)[0]

interface = prism.interface.PYSCF(mf, mc, backend = 'opt_einsum').density_fit()
mr_adc = prism.mr_adc.MRADC(interface)
mr_adc.method = "mr-adc(2)-x"
mr_adc.method_type = "cvs-ip"
mr_adc.max_cycle = 100
mr_adc.ncvs = 3
mr_adc.nroots = 8

e_diff, intensity, x = mr_adc.kernel()

plot(e_diff, intensity, broadening = 0.5, omega_min = 540, omega_max = 550, plot = True, x_label = "Energy, eV", y_label = "Intensity", title = "XPS spectrum", filename = "mr_adc_S0")

mc = mcscf.CASCI(mf, 6, 6).density_fit()
emc = mc.casci(mo_coeff=mo, ci0=state2_ci)[0]

interface = prism.interface.PYSCF(mf, mc, backend = 'opt_einsum').density_fit()
mr_adc = prism.mr_adc.MRADC(interface)
mr_adc.method = "mr-adc(2)-x"
mr_adc.method_type = "cvs-ip"
mr_adc.max_cycle = 100
mr_adc.ncvs = 3
mr_adc.nroots = 8

e_diff, intensity, x = mr_adc.kernel()

plot(e_diff, intensity, broadening = 0.5, omega_min = 540, omega_max = 550, plot = True, x_label = "Energy, eV", y_label = "Intensity", title = "XPS spectrum", filename = "mr_adc_S1")


mc = mcscf.CASCI(mf, 6, 6).density_fit()
emc = mc.casci(mo_coeff=mo, ci0=state3_ci)[0]

interface = prism.interface.PYSCF(mf, mc, backend = 'opt_einsum').density_fit()
mr_adc = prism.mr_adc.MRADC(interface)
mr_adc.method = "mr-adc(2)-x"
mr_adc.method_type = "cvs-ip"
mr_adc.max_cycle = 100
mr_adc.ncvs = 3
mr_adc.nroots = 8

e_diff, intensity, x = mr_adc.kernel()

plot(e_diff, intensity, broadening = 0.5, omega_min = 540, omega_max = 550, plot = True, x_label = "Energy, eV", y_label = "Intensity", title = "XPS spectrum", filename = "mr_adc_S2")


