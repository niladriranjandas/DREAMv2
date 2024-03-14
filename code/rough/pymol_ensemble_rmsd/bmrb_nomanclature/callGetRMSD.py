#!/usr/bin/env python
"""usage: pymol -cq callSplitnSave.py -- <ensemble  pdb name without .pdb extn> <folder_name>
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
