import numpy as np


def returnStats(iplist):
    '''

    '''
    #std_all = np.std(np.array(iplist),axis=0)
    #quantile = 3

    #listnew = [i for i in iplist if i <= quantile * std_all]
    #if len(listnew):
    #   ret_std = np.std(np.array(listnew),axis=0)
    #   ret_mean = np.mean(np.array(listnew),axis=0)
    #else:
    #	ret_std = 0
    #	ret_mean = 0
    
    #if len(iplist) > 0 and ret_mean == 0 :
    ret_mean = np.mean(np.array(iplist),axis=0)
    ret_std  = np.std(np.array(iplist),axis=0)

    return ret_std, ret_mean

def parseViolFile(xplorviolfile, opprefix):
    '''

    '''

    strong=2.7
    medium=3.5
    weak=5.5

    #tol=0.5  # for rest
    tol=0.1   # for xplor and spros_em_wref

    pseudo_strong_file = opprefix+'_pseudo_strong.voils'
    pseudo_medium_file = opprefix+'_pseudo_meidum.voils'
    pseudo_weak_file = opprefix+'_pseudo_weak.voils'

    nonpseudo_strong_file = opprefix+'_nonpseudo_strong.voils'
    nonpseudo_medium_file = opprefix+'_nonpseudo_meidum.voils'
    nonpseudo_weak_file = opprefix+'_nonpseudo_weak.voils'
   
    pseudo_viol_strong = []
    pseudo_viol_medium = []
    pseudo_viol_weak   = []
    nonpseudo_viol_strong = []
    nonpseudo_viol_medium = []
    nonpseudo_viol_weak   = []    
    

    with open(xplorviolfile) as fp, open(pseudo_strong_file,'w') as p_s, open(pseudo_medium_file,'w') as p_m, open(pseudo_weak_file,'w') as p_w, open(nonpseudo_strong_file,'w') as np_s, open(nonpseudo_medium_file,'w') as np_m, open(nonpseudo_weak_file,'w') as np_w:
         lines = fp.readlines()
         for i in range(0, len(lines)):
             pseudo_flag = 0
             line = lines[i]
             if line.strip() :                 
                 first_word = line.split()[0]
                 if (first_word.isnumeric()) or (first_word == '*'):
                        if(i<len(lines)-1):
                           if (lines[i+1].strip()):
                             nextline_firstword = lines[i+1].split()[0]
                             if ( nextline_firstword == '(' ):
                                  pseudo_flag = 1

                        parts = line.split(')')
                        resi_atom_i = parts[0].split('(')[1]
                        resi_i = resi_atom_i.split()[0]
                        resi_atom_j = parts[1].split('(')[1]
                        resi_j = resi_atom_j.split()[0]
                    
                        dist_info = parts[2].split()
                        calc_dist = float(dist_info[0])
                        given_dist_range = dist_info[1].split('..')
                        given_dist_lo = float(given_dist_range[0])
                        given_dist_hi = float(given_dist_range[1])
                        diff = float(dist_info[2])
                        
                        ub_diff = calc_dist - given_dist_hi
                        if (ub_diff > tol) :
                            if (given_dist_hi > 0) and  (given_dist_hi <=  strong) :
                                if (pseudo_flag == 1) :
                                     pseudo_viol_strong.append(diff)
                                     p_s.write(line)
                                else :
                                     nonpseudo_viol_strong.append(diff)
                                     np_s.write(line)

                            if (given_dist_hi > strong) and (given_dist_hi <= medium) :
                                if (pseudo_flag == 1) :
                                     pseudo_viol_medium.append(diff)
                                     p_m.write(line)
                                else :
                                     nonpseudo_viol_medium.append(diff)
                                     np_m.write(line) 

                            if (given_dist_hi > medium) and (given_dist_hi <= weak) :
                                if (pseudo_flag == 1) :
                                     pseudo_viol_weak.append(diff)
                                     p_w.write(line)
                                else :
                                     nonpseudo_viol_weak.append(diff)
                                     np_w.write(line) 
    
    if len(pseudo_viol_strong) :
        ps_ret_std, ps_ret_mean = returnStats(pseudo_viol_strong)
    else:
    	ps_ret_std = 0
    	ps_ret_mean = 0

    if len(pseudo_viol_medium) :
        pm_ret_std, pm_ret_mean = returnStats(pseudo_viol_medium)
    else:
    	pm_ret_std = 0
    	pm_ret_mean = 0

    if len(pseudo_viol_weak) :
        pw_ret_std, pw_ret_mean = returnStats(pseudo_viol_weak)
    else:
    	pw_ret_std = 0
    	pw_ret_mean = 0

    if len(nonpseudo_viol_strong) :
        nps_ret_std, nps_ret_mean = returnStats(nonpseudo_viol_strong)
    else:
    	nps_ret_std = 0
    	nps_ret_mean = 0        

    if len(nonpseudo_viol_medium) :
        npm_ret_std, npm_ret_mean = returnStats(nonpseudo_viol_medium)
    else:
    	npm_ret_std = 0
    	npm_ret_mean = 0        

    if len(nonpseudo_viol_weak) :
        npw_ret_std, npw_ret_mean = returnStats(nonpseudo_viol_weak)
    else:
    	npw_ret_std = 0
    	npw_ret_mean = 0

    #print('\n Pseudo: strong: (%d,%f,%f) medium: (%d,%f,%f) weak: (%d,%f,%f)'%(len(pseudo_viol_strong),ps_ret_mean,ps_ret_std,len(pseudo_viol_medium),pm_ret_mean,pm_ret_std,len(pseudo_viol_weak),pw_ret_mean,pw_ret_std))
    #print('\n Non-Pseudo: strong: (%d,%f,%f) medium: (%d,%f,%f) weak: (%d,%f,%f)'%(len(nonpseudo_viol_strong),nps_ret_mean,nps_ret_std,len(nonpseudo_viol_medium),npm_ret_mean,npm_ret_std,len(nonpseudo_viol_weak),npw_ret_mean,npw_ret_std))

    return len(pseudo_viol_strong),ps_ret_mean,ps_ret_std,len(pseudo_viol_medium),pm_ret_mean,pm_ret_std,len(pseudo_viol_weak),pw_ret_mean,pw_ret_std,len(nonpseudo_viol_strong),nps_ret_mean,nps_ret_std,len(nonpseudo_viol_medium),npm_ret_mean,npm_ret_std,len(nonpseudo_viol_weak),npw_ret_mean,npw_ret_std    

def parseViolFile_v2(xplorviolfile, opprefix):
    '''

    '''

    strong=2.7
    medium=3.5
    weak=5.5

    tol=0.5

    pseudo_strong_file = opprefix+'_pseudo_strong.voils'
    pseudo_medium_file = opprefix+'_pseudo_meidum.voils'
    pseudo_weak_file = opprefix+'_pseudo_weak.voils'

    nonpseudo_strong_file = opprefix+'_nonpseudo_strong.voils'
    nonpseudo_medium_file = opprefix+'_nonpseudo_meidum.voils'
    nonpseudo_weak_file = opprefix+'_nonpseudo_weak.voils'
   
    pseudo_viol_strong = []
    pseudo_viol_medium = []
    pseudo_viol_weak   = []
    nonpseudo_viol_strong = []
    nonpseudo_viol_medium = []
    nonpseudo_viol_weak   = []    
    
    p_strong_c = 0
    p_medium_c = 0
    p_weak_c   = 0
    np_strong_c = 0
    np_medium_c = 0
    np_weak_c   = 0

    with open(xplorviolfile) as fp, open(pseudo_strong_file,'w') as p_s, open(pseudo_medium_file,'w') as p_m, open(pseudo_weak_file,'w') as p_w, open(nonpseudo_strong_file,'w') as np_s, open(nonpseudo_medium_file,'w') as np_m, open(nonpseudo_weak_file,'w') as np_w:
         lines = fp.readlines()
         for i in range(0, len(lines)):
             pseudo_flag = 0
             line = lines[i]
             if line.strip() :                 
                 first_word = line.split()[0]
                 if (first_word.isnumeric()) or (first_word == '*'):
                        if(i<len(lines)-1):
                           if (lines[i+1].strip()):
                             nextline_firstword = lines[i+1].split()[0]
                             if ( nextline_firstword == '(' ):
                                  pseudo_flag = 1

                        parts = line.split(')')
                        resi_atom_i = parts[0].split('(')[1]
                        resi_i = resi_atom_i.split()[0]
                        resi_atom_j = parts[1].split('(')[1]
                        resi_j = resi_atom_j.split()[0]
                    
                        dist_info = parts[2].split()
                        calc_dist = float(dist_info[0])
                        given_dist_range = dist_info[1].split('..')
                        given_dist_lo = float(given_dist_range[0])
                        given_dist_hi = float(given_dist_range[1])
                        diff = float(dist_info[2])

                        if (pseudo_flag == 1):
                            if (given_dist_hi <= strong):                            	
                        	    p_strong_c = p_strong_c+1
                            if (given_dist_hi > strong) and (given_dist_hi <= medium) :
                        	    p_medium_c = p_medium_c+1
                            if (given_dist_hi > medium) and (given_dist_hi <= weak) :
                        	    p_weak_c = p_weak_c+1
                        else:
                            if (given_dist_hi <= strong):                            	
                        	    np_strong_c = np_strong_c+1
                            if (given_dist_hi > strong) and (given_dist_hi <= medium) :
                        	    np_medium_c = np_medium_c+1
                            if (given_dist_hi > medium) and (given_dist_hi <= weak) :
                        	    np_weak_c = np_weak_c+1

                        ub_diff = calc_dist - given_dist_hi
                        if (ub_diff > tol) :
                            if given_dist_hi <=  strong :
                                if (pseudo_flag == 1) :
                                     pseudo_viol_strong.append(diff)
                                     p_s.write(line)
                                else :
                                     nonpseudo_viol_strong.append(diff)
                                     np_s.write(line)

                            if (given_dist_hi > strong) and (given_dist_hi <= medium) :
                                if (pseudo_flag == 1) :
                                     pseudo_viol_medium.append(diff)
                                     p_m.write(line)
                                else :
                                     nonpseudo_viol_medium.append(diff)
                                     np_m.write(line) 

                            if (given_dist_hi > medium) and (given_dist_hi <= weak) :
                                if (pseudo_flag == 1) :
                                     pseudo_viol_weak.append(diff)
                                     p_w.write(line)
                                else :
                                     nonpseudo_viol_weak.append(diff)
                                     np_w.write(line) 
    
    if len(pseudo_viol_strong) :
        ps_ret_std, ps_ret_mean = returnStats(pseudo_viol_strong)
    else:
    	ps_ret_std = 0
    	ps_ret_mean = 0

    if len(pseudo_viol_medium) :
        pm_ret_std, pm_ret_mean = returnStats(pseudo_viol_medium)
    else:
    	pm_ret_std = 0
    	pm_ret_mean = 0

    if len(pseudo_viol_weak) :
        pw_ret_std, pw_ret_mean = returnStats(pseudo_viol_weak)
    else:
    	pw_ret_std = 0
    	pw_ret_mean = 0

    if len(nonpseudo_viol_strong) :
        nps_ret_std, nps_ret_mean = returnStats(nonpseudo_viol_strong)
    else:
    	nps_ret_std = 0
    	nps_ret_mean = 0        

    if len(nonpseudo_viol_medium) :
        npm_ret_std, npm_ret_mean = returnStats(nonpseudo_viol_medium)
    else:
    	npm_ret_std = 0
    	npm_ret_mean = 0        

    if len(nonpseudo_viol_weak) :
        npw_ret_std, npw_ret_mean = returnStats(nonpseudo_viol_weak)
    else:
    	npw_ret_std = 0
    	npw_ret_mean = 0

    #print('\n Pseudo: strong: (%d,%f,%f) medium: (%d,%f,%f) weak: (%d,%f,%f)'%(len(pseudo_viol_strong),ps_ret_mean,ps_ret_std,len(pseudo_viol_medium),pm_ret_mean,pm_ret_std,len(pseudo_viol_weak),pw_ret_mean,pw_ret_std))
    #print('\n Non-Pseudo: strong: (%d,%f,%f) medium: (%d,%f,%f) weak: (%d,%f,%f)'%(len(nonpseudo_viol_strong),nps_ret_mean,nps_ret_std,len(nonpseudo_viol_medium),npm_ret_mean,npm_ret_std,len(nonpseudo_viol_weak),npw_ret_mean,npw_ret_std))

    if p_strong_c == 0:
    	part_p_strong = 0
    else:
        part_p_strong = float(len(pseudo_viol_strong)/p_strong_c)
    if p_medium_c == 0:
    	part_p_medium = 0
    else:
    	part_p_medium = len(pseudo_viol_medium)/p_medium_c 
    if p_weak_c == 0:
    	part_p_weak = 0
    else:
        part_p_weak = len(pseudo_viol_weak)/p_weak_c
    if np_strong_c == 0:
    	part_np_strong = 0
    else:
        part_np_strong = len(nonpseudo_viol_strong)/np_strong_c
    if np_medium_c == 0:
    	part_np_medium = 0
    else:
    	part_np_medium = len(nonpseudo_viol_medium)/np_medium_c 
    if np_weak_c == 0:
    	part_np_weak = 0
    else:
        part_np_weak = len(nonpseudo_viol_weak)/np_weak_c


    return part_p_strong,ps_ret_mean,ps_ret_std,part_p_medium,pm_ret_mean,pm_ret_std,part_p_weak,pw_ret_mean,pw_ret_std,part_np_strong,nps_ret_mean,nps_ret_std,part_np_medium,npm_ret_mean,npm_ret_std,part_np_weak,npw_ret_mean,npw_ret_std    
    
