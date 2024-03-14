
from modeller import *
from modeller.automodel import *    # Load the automodel class
log.verbose()
env = environ()

# directories for input atom files
env.io.atom_files_directory = ['.', '../atom_files']
class MyModel(automodel):
       def special_patches(self, aln):
               self.rename_segments(segment_ids='', renumber_residues=109)
       def select_atoms(self):
		return selection(self.residue_range('109','111'),
				self.residue_range('246','249'))


a = MyModel(env, alnfile = 'alignment.ali',
knowns = 'grp_1tfb', sequence = 'grp_1tfb_fill')
a.starting_model= 1
a.ending_model  = 5

a.make()
