#!/usr/bin/env python
"""usage: pymol -cq getViolations.py -- "/home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/gromacs_pos_restraint/2m4k_trr_extract/script_1" 2m4k_1_clusters_0001 /home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/pymol_ensemble_rmsd/2m4k.dist.1.upl_cpy

          pymol -cq getViolations.py -- <path to pdbs> <pdb1,pdb2...> <full_path to upl file> <full path to lo file>
	  can be called without lo file but upl file has to be mentioned
"""

# Tell PyMOL to launch quiet (-q), fullscreen (-e) and without internal GUI (-i)
import __main__
__main__.pymol_argv = [ 'pymol', '-qei' ]

import pymol, sys

# Call the function below before using any PyMOL modules.
pymol.finish_launching()

from pymol import cmd
import os.path
cmd.reinitialize()

print len(sys.argv)
if (len(sys.argv) != 4):
	print "\n usage: pymol -cq runAcoViol.py -- <path to pdbs> <pdb1,pdb2..> <full path aco path>"
	cmd.exit(1)
	sys.exit(1)
    
pdb_base_dir = '%s' %(sys.argv[1])
pdbs = '%s' %(sys.argv[2])
list_pdbs = pdbs.split(",")
aco_file = '%s' %(sys.argv[3])


if os.path.exists(aco_file):
	file_count=0
	sum_rmsd=0
	sum_mean=0
	sum_std=0
	for pdb in list_pdbs:
		print pdb
		if os.path.exists('%s/%s.pdb'%(pdb_base_dir,pdb)):
			file_count += 1
			cmd.load('%s/%s.pdb'%(pdb_base_dir,pdb))
			cmd.do('run acoViolation.py')
        
			getAcoViolCmd = "[mean,std,rmsd,viol_count,count]=getAcoViol('%s','%s',0)" %(pdb,aco_file)
			cmd.do(getAcoViolCmd)
			sum_rmsd += rmsd
			sum_mean += mean
			sum_std  += std
		else:
			print 'Error:%s not found '%(pdb)
        
	if (file_count > 0):
		aco_rmsd = sum_rmsd/file_count
		aco_mean = sum_mean/file_count
		aco_std = sum_std/file_count
	else:
		aco_rmsd = 0
		aco_mean = 0
		aco_std = 0

	print "(rmsd:%f ,mean:%f ,std:%f)" %(aco_rmsd,aco_mean,aco_std)
else:
	print "Error aco file %s does not exist." %(aco_file)
	cmd.exit(1)
	sys.exit(1)
        

