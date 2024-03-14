#!/usr/bin/env python
"""usage: pymol -cq getViolations.py -- "/home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/gromacs_pos_restraint/2m4k_trr_extract/script_1" 2m4k_1_clusters_0001 /home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/pymol_ensemble_rmsd/2m4k.dist.1.upl_cpy
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

pdb_base_dir = '%s' %(sys.argv[1])
pdbs = '%s' %(sys.argv[2])
flag_up = 0
flag_lo = 0
if (len(sys.argv) > 3):
	up_file = '%s' %(sys.argv[3])
	if os.path.exists(up_file):
		flag_up = 1
	else:
		print "\n upper bound file must exist"
		cmd.quit()
		sys.exit(1)
else:
	print "\n up file must exist"
	cmd.quit()
	sys.exit(1)

if (len(sys.argv) > 4):
	lo_file = '%s' %(sys.argv[4])
	if os.path.exists(lo_file):
		flag_lo = 1
#result_folder = '%s' %(sys.argv[4])

sum_up=0
sum_lo=0
count_up=0
count_lo=0
list_pdbs = pdbs.split(',')
for pdb in list_pdbs:
    if os.path.exists('%s/%s.pdb'%(pdb_base_dir,pdb)):
        cmd.load('%s/%s.pdb'%(pdb_base_dir,pdb))
        cmd.do('run findGetViolations.py')

        if flag_up:
            getViolationCmd = "[tmp_up,stv_up,mean_up,num_viol] = getRMSD('%s','%s',1,0)" %(pdb,up_file)
	    cmd.do(getViolationCmd)	    
            sum_up += tmp_up
            count_up += 1
	    print "PDB:%s (%s,%s,%s)" %(pdb,tmp_up,stv_up,mean_up)

        if flag_lo:
            getViolationCmd = "[tmp_lo,stv_lo,mean_lo,num_viol] = getRMSD('%s','%s',2,0)" %(pdb,lo_file)
            cmd.do(getViolationCmd)
	    sum_lo += tmp_lo
            count_lo += 1
	    print "PDB:%s (%s,%s,%s)" %(pdb,tmp_lo,stv_lo,mean_lo)
    else:
        print 'Error:%s not found '%(pdb)

if (count_up > 0):
	rmsd_up_avg = sum_up/count_up
else:
	rmsd_up_avg = 0

if (count_lo > 0):
	rmsd_lo_avg = sum_lo/count_lo
else:
	rmsd_lo_avg = 0

print "(rmsd_up_avg:%s,rmsd_lo_avg:%s,rmsd_up_count:%s,rmsd_lo_count:%s)" %(rmsd_up_avg,rmsd_lo_avg,count_up,count_lo)

cmd.quit()
