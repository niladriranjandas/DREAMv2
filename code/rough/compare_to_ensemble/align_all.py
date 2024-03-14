#! /usr/bin/env python

from pymol import cmd

#def align_all(target=None,org_mol=None,mobile_selection='name ca',target_selection='name ca',cutoff=2, cycles=5, method='align'):
def align_all(target=None,org_mol=None,mobile_selection='backbone',target_selection='backbone',cutoff=2, cycles=5, method='align'):
 """
 Aligns target to all the structures
 online help: http://pldserver1.biochem.queensu.ca/~rlc/work/pymol/align_all.py
 usage: align_all target=grp_registration, org_mol=2KK1
 """
 cutoff = int(cutoff)
 cycles = int(cycles)

 cmd.split_states('%s'%(org_mol))

 object_list = cmd.get_names()
 object_list.remove(target)
 object_list.remove(org_mol)

 rmsd = {}
 rmsd_list = []

 for i in range(len(object_list)):
	if method == 'align':
		rms = cmd.align('%s & %s'%(object_list[i],mobile_selection),'%s &%s'%(target,target_selection),cutoff=cutoff,cycles=cycles)
	elif method == 'super':
		rms = cmd.super('%s & %s'%(object_list[i],mobile_selection),'%s &%s'%(target,target_selection),cutoff=cutoff,cycles=cycles)
	elif method == 'cealign':
		rms = cmd.cealign('%s & %s' %(target,target_selection),'%s & %s' %(object_list[i],mobile_selection))
	else:
		print "methods allowed for alignment: align super cealign"
		sys.exit(-1)
	rmsd[object_list[i]] = (rms[0],rms[1])
	rmsd_list.append((object_list[i],rms[0],rms[1]))

 rmsd_list.sort(lambda x,y: cmp(x[1],y[1]))
# loop over dictionary and print out matrix of final rms values
#print "Aligning against:",target
#for object_name in object_list:
#	print "%s: %6.3f using %d atoms" % (object_name,rmsd[object_name][0],rmsd[object_name][1])

 print "\nSorted from best match to worst:"
 for r in rmsd_list:
	print "%s: %6.3f using %d atoms" % r

 object_name = 'align_to_%s' %(rmsd_list[0][0])
 if method == 'align':
	rms_best = cmd.align('%s & %s'%(rmsd_list[0][0],mobile_selection),'%s &%s'%(target,target_selection),cutoff=cutoff,cycles=cycles,object=object_name)
 elif method == 'super':
	rms_best = cmd.align('%s & %s'%(rmsd_list[0][0],mobile_selection),'%s &%s'%(target,target_selection),cutoff=cutoff,cycles=cycles,object=object_name)
 elif method == 'cealign':
	rms_best = cmd.cealign('%s & %s' %(target,target_selection),'%s & %s' %(rmsd_list[0][0],mobile_selection))

 print rms_best
 
 cmd.hide("all")
 cmd.show("cartoon",target)
 cmd.show("cartoon",rmsd_list[0][0])
 filename = 'pymol_%s.pse' %(org_mol)
 cmd.save(filename)

cmd.extend('align_all',align_all)         
