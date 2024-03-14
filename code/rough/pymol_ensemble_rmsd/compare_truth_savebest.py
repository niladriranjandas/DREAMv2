#!/usr/bin/env python
"""
call a
usage:pymol -cq compare_truth_savebest.py -- <gen_pdb_name without .pdb extn> <path to gen_pdb> <ground truth pdb without .pdb extn> <path to ground truth> <name of pse file to save without .pse>
      pymol -cq compare_truth_savebest.py -- grp_registration /home/niladri/Documents/Disco_etc_all_in_1/SPROS_try/localization_try/pdb_resi 2KK1 /home/niladri/Documents/Disco_etc_all_in_1/SPROS_try/proteins/2kk1 2kk1
"""
 
# Tell PyMOL to launch quiet (-q), fullscreen (-e) and without internal GUI (-i)
import __main__
__main__.pymol_argv = [ 'pymol', '-qei' ]
 
import pymol, sys

# Call the function below before using any PyMOL modules.
pymol.finish_launching()
 
from pymol import cmd
cmd.reinitialize()
gen_res = '%s/%s.pdb' %(sys.argv[2],sys.argv[1])
truth = '%s/%s.pdb' %(sys.argv[4],sys.argv[3])
cmd.load(gen_res)
cmd.load(truth)
cmd.do('run align_all_savebest.py')
align_all_savebest_cmd = 'align_all_savebest target=%s, org_mol=%s, save_pse=%s' %(sys.argv[1],sys.argv[3],sys.argv[5])
#cmd.do('align_all_savebest target=grp_registration, org_mol=2KK1')
cmd.do(align_all_savebest_cmd)
