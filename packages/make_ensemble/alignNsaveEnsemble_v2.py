"""usage: pymol -cq alignNsaveEnsemble_v2.py -- 1.pdb,2.pdb,3.pdb,4.pdb,5.pdb try.pdb
"""

# Tell PyMOL to launch quiet (-q), fullscreen (-e) and without internal GUI (-i)
import __main__
__main__.pymol_argv = [ 'pymol', '-qei' ]

import pymol, sys, os.path

def alignNsave(oppPDBprefix):
    '''
    
    '''
    
    listofobjs = cmd.get_object_list()
    arg_str = ' '.join(listofobjs)
    cmd.extra_fit('name CA','*','super')
    cmd.join_states('ensemble','*',0)
    cmd.multisave(oppPDBprefix,'ensemble',0)
    
# Call the function below before using any PyMOL modules.
pymol.finish_launching()

cmd.reinitialize()

allpdbs = sys.argv[1]
oppdb   = sys.argv[2]

list_allpdbs = allpdbs.split(',')

for pdb in list_allpdbs:
     cmd.load(pdb)

alignNsave(oppdb)

cmd.quit()


