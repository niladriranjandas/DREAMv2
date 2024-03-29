#
# Python module "average3d.py" for averaging 3D coordinates of multiple states 
# of a single PyMOL molecular object and storing residue-specific RMSDs for the
# bundle of structures in the b-factor field of the resulting averaged model.  
# Could be used in conjuction with "putty"-style PyMOL cartoons to illustrate
# residue-specific RMSDs. 
#

"""
AUTHOR:
	Cameron Mura <cmura@ucsd.edu>; 07/2005

SYNOPSIS:
	Python module "average3d.py" for averaging 3D coordinates of multiple states 
	of a single PyMOL molecular object and storing residue-specific RMSDs for the
	bundle of structures in the b-factor field of the resulting averaged model.  
	Could be used in conjuction with "putty"-style PyMOL cartoons to illustrate
	residue-specific RMSDs. 

USAGE:
	See this source file for now... in particular, any notes under individual 
	function definitions, such as 'avgRMS()'.

"""

from pymol import cmd
from math import sqrt,pow
from string import replace

def avgRMS(object='all',object_sel=None,first=1,last=0,newobj=None,fitverdict='no',\
		verb=0):
	'''
	ARGUMENTS:
	object (string [defaults to 'all']):
		Starting PyMOL molecular object consisting of >1 states to be averaged
		
	object_sel (string [defaults to "name CA"]):
		An optional subset of atoms on which to operate (e.g., "name CA and resi 20-50")

	first (int [defaults to 1]):
		number of the first state to include in averaging

	last (int [defaults to last state]):
		Number of the last state to include in averaging. A value of '0' (zero) 
		means use the last state in the object.
	
	newobj (string [defaults to name of object with string "_avg_AtoZ" appended, 
			where A and Z are the numbers of the first and last frames included 
			in the averaging]):
		Desired name of the new output pymol object (i.e., averaged structure)
	
	fitverdict (string [defaults to "no", any value != "no" is taken as a "yes"]):
		Use the pymol function intra_fit() to fit all the states of "object & object_sel" 
		before calculating the average structure and RMSDs??
		'no' = NO, !'no' = yes
	
	verb (int [defaults to 0]):
		Verbosity level of output. 0 = quiet, !0 = verbose (e.g., write position-specific 
		RMSDs to standard output).
	
	NOTES:	
	  usage: avgRMS('1T3V','backbone',1,0,'yesfit_ALL','yes',0)
	  The state specified as 'first' is taken as the reference state for the (optional)
	  intra_fit() alignment. If no selection is given, the string "and name CA" will be 
	  appended to the object and the averaging and rmsd calculations will be restricted 
	  to "object and name CA" (i.e., C-alpha atoms). To avoid this default behavior, 
	  specify a valid pymol atom selection for the object_sel argument (e.g., "name CA" 
	  or "name CA and resi 20-50" or whatever).
	'''
	
	object = str(object)
	object_sel = str(object_sel)
	first=int(first)
	last=int(last)
	num_states_tot = cmd.count_states(object)
	if last<1:	
		last=num_states_tot
	num_states2avg = last - first + 1
	if newobj==None or newobj=='':
		newobj = "%s_avg_%dto%d"%(object,first,last)
	cmd.create(newobj,object,first,1)	

	rmsd_file_name = '%s_rmsd_file.txt'%(object)
	#rmsd_file_name = '%s_rmsd_file_ca.txt'%(object)
	rmsd_file = open(rmsd_file_name,'w')	
		
	obj_sel_string = (lambda o,s: 	(o and s and "%s and %s"%(o,s)) or \
									(o and not s and "%s and name CA"%o))(object,object_sel)
	newobj_sel_string = (lambda o,s: 	(o and s and "%s and %s"%(o,s)) or \
										(o and not s and "%s and name CA"%o))(newobj,object_sel)
	print '%s'%('-'*80)
	print 'Averaging %d states [%d, %d] of object "%s" (%d total states) to get new single-state object "%s"...' \
			% (num_states2avg,first,last,obj_sel_string,num_states_tot,newobj)
	print '%s'%('-'*80)
	
	tmpobject='%s_tmp4avg_%dto%d'%(object,first,last)
	i=0
	for eachstate in range(first,last+1):
		i+=1	# because pymol states are indexed starting from 1 (not 0)
		cmd.create(tmpobject,obj_sel_string,eachstate,i)
		
	if fitverdict != 'no':
		print "-> proceeding WITH FITTING to reference state..."
		tmpobject_intrafit_rmsds = cmd.intra_fit(tmpobject,1)
		if verb:
			print '   intrafit_rmsds = ',
			print tmpobject_intrafit_rmsds
	else:
		print "-> proceeding WITHOUT FITTING to reference state..."

	# create an atom index map between original (complete) object and subset to be averaged:
	atindex_map = [None]
	for at in cmd.get_model(newobj_sel_string,1).atom:
		atindex_map.append(at.index)
	if verb:
		print "-> atom index_map = ",
		print atindex_map

	newobj_chempy = cmd.get_model(tmpobject,1)
	for at in newobj_chempy.atom:
		this_at_idx = at.index
		sum_x = sum_y = sum_z = 0.0
		avg_x = avg_y = avg_z = 0.0
		state_coords=[]
		rmsd_sum = rmsd = 0.0
		for s in range(1,num_states2avg+1):
			this_x = cmd.get_model(tmpobject,s).atom[this_at_idx-1].coord[0]
			this_y = cmd.get_model(tmpobject,s).atom[this_at_idx-1].coord[1]
			this_z = cmd.get_model(tmpobject,s).atom[this_at_idx-1].coord[2]
			sum_x += this_x
			sum_y += this_y
			sum_z += this_z
			state_coords.append([this_x,this_y,this_z])
		avg_x = sum_x / num_states2avg
		avg_y = sum_y / num_states2avg
		avg_z = sum_z / num_states2avg
		cmd.alter_state(1,'%s and id %d'%(newobj,atindex_map[this_at_idx]),'x=%f'%avg_x)
		cmd.alter_state(1,'%s and id %d'%(newobj,atindex_map[this_at_idx]),'y=%f'%avg_y)
		cmd.alter_state(1,'%s and id %d'%(newobj,atindex_map[this_at_idx]),'z=%f'%avg_z)
		
		# position-specific RMSDs with respect to average structure:
		#for somept in state_coords:
		#	rmsd_sum += pow(calcDist(somept,[avg_x,avg_y,avg_z]),2.0)
		#rmsd = sqrt(rmsd_sum / num_states2avg)
		## DEBUG:
		## DEBUG print 'newobj = %s / id = %d / b = %0.3f'%(newobj,atindex_map[this_at_idx],rmsd)
		#cmd.alter('%s and id %d'%(newobj,atindex_map[this_at_idx]),'b=%0.3f'%rmsd)
		#if verb:
		#	print '"%s" index = %d'%(obj_sel_string,this_at_idx)
		#	print 'rmsd = %0.3f'%rmsd
		#avg_file.write('%d\t%0.3f\n'%(atindex_map[this_at_idx],rmsd))
        
        # added by me to save the average object
	#cmd.save('%s.pdb'%(newobj),newobj,0)
        #all_avg_rmsd = 0
        print '\n----rmsd---\n'
        avgobj = cmd.get_object_list(newobj)
        cmd.split_states(object)
	#print '\n%s'%zz
	myobjects = cmd.get_object_list()
	myobjects.remove(object)
	myobjects.remove(newobj)
	myobjects.remove(tmpobject)	
	#print myobjects
        for en in myobjects:		
		cycles=5
		cutoff=2
                atom_selection = 'backbone'                		
               #rms = cmd.align(curr_mol,newobj,cutoff=cutoff,cycles=cycles)
		print '\n %s to %s'%(en,avgobj[0])
                rms = cmd.align('%s & %s'%(en,object_sel),'%s &%s'%(avgobj[0],object_sel),cutoff=cutoff,cycles=cycles)
	        #all_avg_rms += rms
		print>>rmsd_file, rms
                print rms                
 
	rmsd_file.close()

def calcDist(vec1,vec2):
	sum=0.0
	for (coor1,coor2) in zip(vec1,vec2): sum += pow((coor1-coor2),2.0)
	return sqrt(sum)


# make it directly callable in PyMOL:
cmd.extend('avgRMS',avgRMS)

# scratch space:
# cmd.rebuild()	
# cmd.frame(first)
