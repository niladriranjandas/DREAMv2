
"""
 get the pdb 
 get the distance bounds and calculate the violations
"""

from pymol import cmd
from math import sqrt,pow
from string import replace
import os.path

def getRMSD(pdb_obj, bound_file, up_lo_flag, rmsd):
	'''
	ARGUMENTS:
	'''
	
	tol=0.2
	rmsd = 0
	mu_x_1 = 0
	mu_x_2 = 0
	var = 0
	sdv = 0
	sum_1 = 0
	sum_2 = 0
	count = 0
	violCount = 0
	max_viol = 0
	with open(bound_file) as f:
		for line in f:			
			cols = line.split()			

			resi_i = cols[0]
			resi_i_n = cols[1]
			resi_i_a = cols[2]

			resi_j = cols[3]
			resi_j_n = cols[4]
			resi_j_a = cols[5]
	
			dist_ij = float(cols[6])
			
			len_resi_i_a = len(resi_i_a)
			#if (len_resi_i_a > 3):
			#	tmp_1 = resi_i_a[:len_resi_i_a - 1]
			#	tmp_2 = resi_i_a[len_resi_i_a-1]
			#	resi_i_a = '%s%s'%(tmp_2,tmp_1)
				#tmp_resi_i_a = '%s%s'%(tmp_2,tmp_1)
				#print "***********%s"%(tmp_resi_i_a)
			#len_resi_j_a = len(resi_j_a)
			#if (len_resi_j_a > 3):
			#	tmp_1 = resi_j_a[:len_resi_j_a - 1]
			#	tmp_2 = resi_j_a[len_resi_j_a-1]
			#	resi_j_a = '%s%s'%(tmp_2,tmp_1)
				#tmp_resi_j_a = '%s%s'%(tmp_2,tmp_1)
				#print "----------%s"%(tmp_resi_j_a)
			i_str = '/%s///%s/%s'%(pdb_obj,resi_i,resi_i_a)
			j_str = '/%s///%s/%s'%(pdb_obj,resi_j,resi_j_a)

			#print i_str
			#print j_str

			try:
				dist_calc = cmd.get_distance(i_str,j_str)				
				count = count + 1
				if ( up_lo_flag == 1 ):
					#count = count + 1
					if (dist_calc > dist_ij ):
						diff = dist_calc - dist_ij
						if (diff > tol ):
							sum_1 = sum_1 + diff
							sum_2 = sum_2 + pow(diff,2)
							#count = count + 1
							violCount = violCount + 1
							if (diff > max_viol):
								max_viol = diff
					
				elif ( up_lo_flag == 2 ):
					#count = count + 1
					if (dist_calc < dist_ij ):
						diff = dist_ij - dist_calc
						if (diff > tol ):
							sum_1 = sum_1 + diff
							sum_2 = sum_2 + pow(diff,2)
							#count = count + 1
							violCount = violCount + 1
							if (diff > max_viol):
								max_viol = diff
			except:
				print 'error:(%s,%s)'%(i_str,j_str)

	if ( count > 0 ):
		#mu_x_1 = sum_1/count
		#mu_x_2 = sum_2/count
		mu_x_1 = sum_1/violCount
		mu_x_2 = sum_2/violCount
		print "mu_x1: %s, mu_x2: %s , count=%s" %(mu_x_1,mu_x_2,count)
		var = mu_x_2 - pow(mu_x_1,2)
		sdv = sqrt(var)
		rmsd = sqrt(sum_2/violCount)	
		#rmsd = sqrt(sum_2/count)	
	print "(RMSD,stdev,mean,violCount,maxViol):getRMSD:(%f,%f,%f,%f,%f)"%(rmsd,sdv,mu_x_1,violCount,max_viol)
	return rmsd,sdv,mu_x_1,violCount

# make it directly callable in PyMOL:
cmd.extend('getRMSD',getRMSD)


	

