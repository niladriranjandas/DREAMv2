
import os
from pymol import cmd


def splitnSave(object,dest_folder):
	'''
	splits the input pdbs into the models and saves the pdb's separately. Also returns the names of pdb
	bundles as a list
 	     object: pdb name without pdb extension (after loaded in pymol
	dest_folder: folder where the individual pdb's have to be saved
	'''				
	object = str(object)
	names = ""
	if os.path.isdir(dest_folder):
		states = cmd.count_states(object)
		if (states > 1):
			cmd.split_states(object)
			names = cmd.get_names()	
			print names
			names.remove(object)

			for name_i in names:
				print name_i
				file_i = '%s/%s.pdb' %(dest_folder,name_i)
				cmd.save(file_i,name_i)
		else:
			print "Number of states in %s is 1" %(object)
	else:
		print "Error: destination folder doesn't exist %s" %(dest_folder)

	return names
# make it directly callable in PyMOL:
cmd.extend('splitnSave',splitnSave)

