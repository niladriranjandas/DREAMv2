
"""
 get the pdb 
 get the distance bounds and calculate the violations
"""

from pymol import cmd
from math import sqrt,pow
from string import replace
import os.path

def getAcoViol(pdb_obj, aco_file, rmsd):
	'''
	ARGUMENTS
	'''
	rmsd=0
	diff=0
	sum1=0
	sum2=0
	count=0
	viol_count=0
	
	with open(aco_file) as f:
		for line in f:
			count = count + 1
			cols = line.split()

			resi_no = int(cols[0])
			resi_nm = cols[1]
			phi_psi = cols[2]
			lo_ang  = float(cols[3])
			hi_ang  = float(cols[4])

			if ( phi_psi == "PHI" ): #| phi_psi == "phi" | ph_psi == "Phi" ):
				n_i_1 = resi_no - 1
				atm_1 = '%s///%d/c'%(pdb_obj,n_i_1)
				atm_2 = '%s///%d/n'%(pdb_obj,resi_no)
				atm_3 = '%s///%d/ca'%(pdb_obj,resi_no)
				atm_4 = '%s///%d/c'%(pdb_obj,resi_no)				
			elif ( phi_psi == "PSI" ): #| phi_psi == "psi" | phi_psi == "Psi" ):
				n_i_1 = resi_no + 1
				atm_1 = '%s///%d/n'%(pdb_obj,resi_no)
				atm_2 = '%s///%d/ca'%(pdb_obj,resi_no)
				atm_3 = '%s///%d/c'%(pdb_obj,resi_no)
				atm_4 = '%s///%d/n'%(pdb_obj,n_i_1)				
			try:							
				dihed = cmd.get_dihedral(atm_1,atm_2,atm_3,atm_4)				
				if ( dihed < lo_ang ):
					viol_count = viol_count + 1
					diff = lo_ang - dihed
					sum1 = sum1 + diff
					sum2 = sum2 + pow(diff,2)
					print "---------\n"
					print '(resi_no:%d dihed:%f lo_ang:%f diff:%f)'%(resi_no,dihed,lo_ang,diff)
				elif (dihed > hi_ang ):
					viol_count = viol_count + 1
					diff = hi_ang - dihed
					sum1 = sum1 + diff
					sum2 = sum2 + pow(diff,2)
					print "---------\n"
					print '(resi_no:%d dihed:%f hi_ang:%f diff:%f)'%(resi_no,dihed,hi_ang,diff)
			except:
				#traceback.print_stack()
				print 'dihed error:%s,%s,%s,%s'%(atm_1,atm_2,atm_3,atm_4)
	
	rmsd = sqrt(sum2/count)
	mean = sum1/count
	std  = sqrt(sum2/count - pow(mean,2))
	print '(mean: %f, std: %f, rmsd: %f, violcount: %f, count: %f)'%(mean,std,rmsd,viol_count,count)
	return mean,std,rmsd,viol_count,count
# make it directly callable in PyMOL:

cmd.extend('getAcoViol',getAcoViol)			
