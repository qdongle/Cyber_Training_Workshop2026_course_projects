import numpy as np
from pyscf import gto, scf, cc
from pyscf import lib
from pyscf.tools import cubegen
from pyscf.tools import molden

mol = gto.Mole()
#Geometry for napthalene
mol.atom = '''
C 2.41229 0.75716 0.00000
C 2.44089 -0.65953 0.00000
C 1.27099 -1.37135 0.00000
C 0.01433 -0.70687 -0.00000
C -0.01412 0.70690 -0.00000
C 1.21477 1.42144 0.00000
C -1.27093 1.37111 -0.00000
C -2.44092 0.65942 -0.00000
C -2.41252 -0.75728 -0.00000
C -1.21473 -1.42112 -0.00000
H 3.33396 1.30672 0.00000
H 3.38393 -1.17155 0.00000
H 1.29099 -2.44499 0.00000
H 1.19157 2.49502 0.00000
H -1.29097 2.44474 -0.00000
H -3.38391 1.17153 -0.00000
H -3.33414 -1.30692 -0.00000
H -1.19148 -2.49469 -0.00000
'''
mol.basis = 'aug-cc-pvdz'
mol.verbose = 4
mol.charge = 0
mol.build()

mf = scf.RHF(mol).density_fit("aug-cc-pvdz-jkfit")
mf.kernel()

mycc = cc.CCSD(mf).set_frozen().density_fit("aug-cc-pvdz-ri")
mycc.verbose = 4
mycc.kernel()

eS = mycc.eomee_ccsd_singlet(nroots=10)[0]

for state in range(len(eS)):
    print ("State #%d, dE = %9.6f (a.u.), %6.3f (eV)" % (state+1, eS[state], eS[state] * 27.2114))
