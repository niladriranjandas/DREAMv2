
"""
 get the pdb 
 get the distance bounds and calculate the violations
"""

from pymol import cmd
from math import sqrt,pow
from string import replace
import numpy
import re
import os.path
import subprocess 
import shlex

from scipy.stats import norm
import matplotlib.mlab as mlab
#import matplotlib.pyplot as plt

#import pylab as plb
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
#from scipy import asarray as ar,exp
#from numpy import linspace
#import math

def gaussian(x, mean, amplitude, standard_deviation):
    return amplitude * numpy.exp( - ((x - mean) / standard_deviation) ** 2)

def gaus(x,a,sigma):
    return a*numpy.exp(-(x)**2/(2*sigma**2))

def getDist(coord1,coord2):
	'''
	ARGUMENTS:
	'''
	if ( coord1 is None or coord2 is None ):
		return -1
	if ( numpy.count_nonzero(coord1)==0 or numpy.count_nonzero(coord2)==0):
		return -1

	x=coord1[0] - coord2[0]
	y=coord1[1] - coord2[1]
	z=coord1[2] - coord2[2]

	tmp=pow(x,2)+pow(y,2)+pow(z,2)
	return math.sqrt(tmp)
	

def getPseudoCoord(pdb_obj,resi_name,resi_no,atm_name):
	'''
	ARGUMENTS:
	'''
	#print '/%s///%s/%s'%(pdb_obj,resi_no,atm_name)
	if "Q" not in atm_name:
		crd=[]
		getcrdcmd='/%s///%s/%s'%(pdb_obj,resi_no,atm_name)
		crd=cmd.get_coords(getcrdcmd,1)
		if crd is None:
			print '!!!! NOT FOUND %s %s'%(resi_no,atm_name)
			return [1]
		else:
			return crd[0]

	match_atms=""
	for line in open('pseudo.csv'):
		cols=line.split()
		if ( cols[0] == resi_name.lower() and cols[1] == atm_name ):
			match_atms = cols[2]

	atmnames=match_atms
	if ( atmnames != "" ):
		atm_nm = atmnames.split(",")
		count_at=0
		sum_coord =numpy.zeros(3)
		for atn_i in atm_nm:
			argsp=""
			outputp=""
			argsp = shlex.split('./getPDBFromBMRB_edit.sh %s %s'%(resi_name,atn_i))
			outputp,errorp = subprocess.Popen(argsp,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
			match_atm_i = outputp.rstrip()
			getcrdcmd='/%s///%s/%s'%(pdb_obj,resi_no,match_atm_i)
			crd=[]
			crd = cmd.get_coords(getcrdcmd,1)
			if (crd is None):
				print getcrdcmd
				print "could not find a match"
			else:
				coord_i=crd[0]
				sum_coord += crd[0]
				count_at+=1

		if (count_at!=0):
			return numpy.divide(sum_coord,count_at)
		else:
			return numpy.zeros(3)			
	else:
		return None

def getRMSD(pdb_obj, bound_file, up_lo_flag, rmsd):
	'''
	ARGUMENTS:
	'''
	tol = 0.5	
	#tol = 0

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
	viol_arr = []
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

			i_str = '/%s///%s/%s'%(pdb_obj,resi_i,resi_i_a)
			j_str = '/%s///%s/%s'%(pdb_obj,resi_j,resi_j_a)

			args=""
			output1=""
			args = shlex.split('./getPDBFromBMRB_edit.sh %s %s'%(resi_i_n,resi_i_a))
			output1,error1 = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()
				
			args=""
			output2=""
			args = shlex.split('./getPDBFromBMRB_edit.sh %s %s'%(resi_j_n,resi_j_a))
			output2,error2 = subprocess.Popen(args,stdout = subprocess.PIPE, stderr= subprocess.PIPE).communicate()

			dist_calc = -1
			if ( (output1.rstrip() == "XX") or (output2.rstrip() == "XX") ):
				if ( output1.rstrip() == "XX" ):
					coord_i = getPseudoCoord(pdb_obj,resi_i_n,resi_i,resi_i_a)
				else:
					coord_i = getPseudoCoord(pdb_obj,resi_i_n,resi_i,output1.rstrip())

				if ( output2.rstrip() == "XX" ):
					coord_j = getPseudoCoord(pdb_obj,resi_j_n,resi_j,resi_j_a)
				else:
					coord_j = getPseudoCoord(pdb_obj,resi_j_n,resi_j,output2.rstrip())
				
				if ( (len(coord_i) == 1 ) or (len(coord_j) == 1) ):
					dist_calc = -1
				else: 		
					dist_calc=getDist(coord_i,coord_j)
			else:
				i_str = '/%s///%s/%s'%(pdb_obj,resi_i,output1.rstrip())
				j_str = '/%s///%s/%s'%(pdb_obj,resi_j,output2.rstrip())

				try:
					dist_calc = cmd.get_distance(i_str,j_str)				
				except:
					print '%s\n%s\n'%(i_str,j_str)
			if ( dist_calc != -1 ):
				count = count + 1
				if ( up_lo_flag == 1 ):
					#count = count + 1
					#if (dist_calc > dist_ij ):
					if (dist_calc > (dist_ij + tol) ):						
						diff = dist_calc - dist_ij
						print '\n----dist_viol %s \t %s \t %s'%(i_str,j_str,diff)
						sum_1 = sum_1 + diff
						sum_2 = sum_2 + pow(diff,2)
						#count = count + 1
						violCount = violCount + 1
						viol_arr.append(diff)
						if (diff > max_viol):
							max_viol = diff
				elif ( up_lo_flag == 2 ):
					#count = count + 1
					if (dist_calc < dist_ij ):
						diff = dist_ij - dist_calc
						sum_1 = sum_1 + diff
						sum_2 = sum_2 + pow(diff,2)
						#count = count + 1
						violCount = violCount + 1
						viol_arr.append(diff)
						if (diff > max_viol):
							max_viol = diff
			else:
				print "dist could not be calculated"
	if ( count > 0 ):
		mu_x_1 = sum_1/count
		mu_x_2 = sum_2/count
		#mu_x_1 = sum_1/violCount
		#mu_x_2 = sum_2/violCount
		print "mu_x1: %s, mu_x2: %s , count=%s" %(mu_x_1,mu_x_2,count)
		var = mu_x_2 - pow(mu_x_1,2)
		sdv = sqrt(var)
		#rmsd = sqrt(sum_2/violCount)
		rmsd = sqrt(sum_2/count)	
	print "(RMSD,stdev,mean,violCount,maxViol):getRMSD:(%f,%f,%f,%f,%f)"%(rmsd,sdv,mu_x_1,violCount,max_viol)
	print len(viol_arr)
	tmp = range(1,len(viol_arr)+1)
	print len(tmp)
	#thefile = open('test_array.txt','w')
	#for item in viol_arr:
	#	thefile.write("%s\t" %item)
	#plt.hist(viol_arr)
	#plt.show()
	
	###tosave = plt.figure()
	# best fit of data
	#mu, std = norm.fit(viol_arr)
	#plt.hist(viol_arr, normed=True)
	#xmin, xmax = plt.xlim()
	#x = numpy.linspace(xmin, xmax, len(viol_arr))
	#p = norm.pdf(x,mu,std)
	#plt.plot(x,p,'k',linewidth=2)
	#plt.show()
	#print "mu=%f, std=%f" %(mu,std)
	(mu, sigma) = norm.fit(viol_arr)
	print " mu=%f, sigma=%f" %(mu,sigma)
	# the histogram of the data	
	n, bins, patches = plt.hist(viol_arr,normed=1,color='cyan')
	y = mlab.normpdf( bins, mu, sigma)
	#l = plt.plot(bins, y, 'r--', linewidth=2, label='norm-fit')
	#plt.show()

	#popt,pcov = curve_fit(gaus,ar(viol_arr),tmp,p0=None)
	#popt,pcov = curve_fit(gaus,linspace(1,len(viol_arr)+1,len(viol_arr)),ar(sorted(viol_arr)),p0=None)
	#plt.plot(gaus(ar(sorted(viol_arr)),*popt))
	#print popt
	#print pcov
	#plt.show()

	###bin_heights, bin_borders, _ = plt.hist(viol_arr, normed=True)
	###bin_centers = bin_borders[:-1] + numpy.diff(bin_borders) / 2
	###popt, _ = curve_fit(gaus, bin_centers, bin_heights, p0=None)

	###x_interval_for_fit = numpy.linspace(bin_borders[0], bin_borders[-1], len(tmp))
	###plt.plot(x_interval_for_fit, gaus(x_interval_for_fit, *popt), label='fit standard normal')

	# Empirical average and variance are computed
	#avg = numpy.mean(viol_arr)
	#var = numpy.var(viol_arr)
	# From that, we know the shape of the fitted Gaussian.
	#print "avg: %f	var: %f" %(avg,var)
	#pdf_x = numpy.linspace(numpy.min(viol_arr),numpy.max(viol_arr),100)
	#pdf_y = 1.0/numpy.sqrt(2*numpy.pi*var)*numpy.exp(-0.5*(pdf_x-avg)**2/var)

	# Then we plot :
	#plt.figure()
	#plt.hist(viol_arr,normed=True)
	#plt.plot(pdf_x,pdf_y,'k--', label='theoretical')
	#plt.legend(("Fit","Data"),"best")
	###plt.legend()
	###plt.suptitle('NOESY upper bound voilation', fontsize=20)
	###plt.ylabel('normalized')
	###plt.xlabel('violations (in amstrong)')
	###plt.show()
	###tosave.savefig('destination_path.eps', format='eps', dpi=1000)

	#print viol_arr
	#thefile.close()
	return rmsd,sdv,mu_x_1,violCount

# make it directly callable in PyMOL:
cmd.extend('getRMSD',getRMSD)


	

