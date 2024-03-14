import protocol
#protocol.loadPDB(file='WRresult_1pbu_multialign_md_B3_3.pdb',deleteUnknownAtoms=True)
protocol.loadPDB(file='PDB_NAME',deleteUnknownAtoms=True)

from ccrPotTools import create_CCRPot
from noePotTools import create_NOEPot

pot = create_NOEPot('all','NOE_CONSTRAINS')
#pot=create_CCRPot('ccr',file='constraints/1pbu_xplor_noe_test.tbl',tauc=2.5e-9)
#pot.setScale(2.3)
pot.setShowAllRestraints(True)

print(pot.calcEnergy())

pot.setShowAllRestraints(True)
print(pot.showViolations())

from ivm import IVM
dyn = IVM()

protocol.torsionTopology(dyn)

from monteCarlo import randomizeTorsions
randomizeTorsions(dyn)
