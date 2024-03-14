#!/usr/bin/env python
"""usage: pymol -cq ensmbleSplitnSave.py -- <ensemble  pdb name without .pdb extn> <folder_name> <upl file> <lo_file>
   in case of more than one upl file append them on below the other in a single upl/lo file
"""

# Tell PyMOL to launch quiet (-q), fullscreen (-e) and without internal GUI (-i)
import __main__
__main__.pymol_argv = [ 'pymol', '-qei' ]

import pymol, sys, os.path

# Call the function below before using any PyMOL modules.
pymol.finish_launching()

from pymol import cmd
import os.path
cmd.reinitialize()

full_name = '%s/%s.pdb' %(sys.argv[2],sys.argv[1])
if os.path.exists(full_name):
	cmd.load(full_name)
	cmd.do('run splitnSave.py')
	splitsave_cmd = "pdbs = splitnSave('%s','%s')" %(sys.argv[1],sys.argv[2])
	cmd.do(splitsave_cmd)
	print pdbs	
	
	flag_up = 0
	flag_lo = 0
	if os.path.exists(sys.argv[3]): # upl file	
		flag_up = 1	
		count_up = 0
		sum_up = 0
		cmd.do('run findGetViolations.py')		
		for pdb_i in pdbs:		
			if flag_up:
		        	getViolationCmd = "[tmp_up,stv_up,mean_up] = getRMSD('%s','%s',1,0)" %(pdb_i,sys.argv[3],1,0)
	    			cmd.do(getViolationCmd)	    
            			sum_up += tmp_up
            			count_up += 1
	    			print "PDB:%s (%s,%s,%s)" %(pdb_i,tmp_up,stv_up,mean_up)

		        if flag_lo:
            			getViolationCmd = "[tmp_lo,stv_lo,mean_lo] = getRMSD('%s','%s',2,0)" %(pdb_i,sys.argv[4],2,0)
            			cmd.do(getViolationCmd)
	    			sum_lo += tmp_lo
            			count_lo += 1
	    			print "PDB:%s (%s,%s,%s)" %(pdb,tmp_lo,stv_lo,mean_lo)	
			else:
				print "warning lo file %s not found." 
		
		avg_rmsd_up = sum_up/count_up
		print "avg_rmsd_up:%s" %(avg_rmsd_up)

	cmd.quit()
else:
	print "Error: file %s doesn't exist." %(full_name)


