#!/usr/bin/env python
"""usage: pymol -cq ensembleAvgRMSD.py -- <folder_name> <ensembla name without .pdb extn>
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

pdb_name = '%s/%s.pdb' %(sys.argv[1],sys.argv[2])
if os.path.exists(pdb_name):
    cmd.load(pdb_name)
    cmd.do('run averageRMS.py')
    #avg_rmsd_cmd = "avgRMS('%s','backbone',1,0,'yesfit_ALL','yes',0)" %(sys.argv[2])
    avg_rmsd_cmd = "avgRMS('%s','name ca',1,0,'yesfit_ALL','yes',0)" %(sys.argv[2])
    cmd.do(avg_rmsd_cmd)
else:
    print("File %snot found" %(pdb_name))
