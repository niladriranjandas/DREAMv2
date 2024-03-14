import os
import sys
import math
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mp
#import plotly.graph_objects as go
import plotly.express as px
import pandas as pd


def defineTol():
    global strong
    strong = 2.7
    global medium
    medium = 3.5
    global weak
    weak = 5.5
    global tol
    tol = 0.5; #0.1

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

def parseViolFile(xplorviolfile, opprefix, destfolder="."):
    '''

    '''

    strong_file = destfolder+'/'+opprefix+'_strong.voils'
    medium_file = destfolder+'/'+opprefix+'_meidum.voils'
    weak_file = destfolder+'/'+opprefix+'_weak.voils'
   
    viol_strong = []
    viol_medium = []
    viol_weak   = []

    with open(xplorviolfile) as fp, open(strong_file,'w') as s, open(medium_file,'w') as m, open(pseudo_weak_file,'w') as w:
         lines = fp.readlines()
         for i in range(0, len(lines)):
             pseudo_flag = 0
             line = lines[i]
             if line.strip() :                 
                 first_word = line.split()[0]
                 if (first_word.isnumeric()) or (first_word == '*'):
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
                            if given_dist_hi <=  strong :
                                viol_strong.append(diff)
                                s.write(line)

                            if (given_dist_hi > strong) and (given_dist_hi <= medium) :
                                viol_medium.append(diff)
                                m.write(line) 

                            if (given_dist_hi > medium) and (given_dist_hi <= weak) :
                                voil_weak.append(diff)
                                w.write(line) 
    
    if len(viol_strong) :
        s_ret_std, s_ret_mean = returnStats(viol_strong)
    else:
    	s_ret_std = 0
    	s_ret_mean = 0

    if len(viol_medium) :
        m_ret_std, m_ret_mean = returnStats(viol_medium)
    else:
    	m_ret_std = 0
    	m_ret_mean = 0

    if len(viol_weak) :
        w_ret_std, w_ret_mean = returnStats(viol_weak)
    else:
    	w_ret_std = 0
    	w_ret_mean = 0

    return len(viol_strong),s_ret_mean,s_ret_std, len(viol_medium),m_ret_mean,m_ret_std, len(viol_weak),w_ret_mean,w_ret_std

def parseViolFile_v2(xplorviolfile, opprefix, destfolder="."):
    '''

    '''

    strong_file = destfolder+'/'+opprefix+'_strong.voils'
    medium_file = destfolder+'/'+opprefix+'_meidum.voils'
    weak_file = destfolder+'/'+opprefix+'_weak.voils'
   
    viol_strong = []
    viol_medium = []
    viol_weak   = []
    
    strong_c = 0
    medium_c = 0
    weak_c   = 0

    with open(xplorviolfile) as fp, open(strong_file,'w') as s, open(medium_file,'w') as m, open(weak_file,'w') as w:
         lines = fp.readlines()
         for i in range(0, len(lines)):
             line = lines[i]
             if line.strip() :                 
                 first_word = line.split()[0]
                 if (first_word.isnumeric()) or (first_word == '*'):

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

                        if (given_dist_hi <= strong):                            	
                      	    strong_c = strong_c+1
                        if (given_dist_hi > strong) and (given_dist_hi <= medium) :
                       	    medium_c = medium_c+1
                        if (given_dist_hi > medium) and (given_dist_hi <= weak) :
                       	    weak_c = weak_c+1

                        ub_diff = calc_dist - given_dist_hi
                        if (ub_diff > tol) :
                            if given_dist_hi <=  strong :
                                viol_strong.append(diff)
                                s.write(line)

                            if (given_dist_hi > strong) and (given_dist_hi <= medium) :
                                 viol_medium.append(diff)
                                 m.write(line)

                            if (given_dist_hi > medium) and (given_dist_hi <= weak) :
                                 viol_weak.append(diff)
                                 w.write(line)
    
    if len(viol_strong) :
        s_ret_std, s_ret_mean = returnStats(viol_strong)
    else:
    	s_ret_std = 0
    	s_ret_mean = 0

    if len(viol_medium) :
        m_ret_std, m_ret_mean = returnStats(viol_medium)
    else:
    	m_ret_std = 0
    	m_ret_mean = 0

    if len(viol_weak) :
        w_ret_std, w_ret_mean = returnStats(viol_weak)
    else:
    	w_ret_std = 0
    	w_ret_mean = 0


    if strong_c == 0:
    	part_strong = 0
    else:
        part_strong = float(len(viol_strong)/strong_c)
    if medium_c == 0:
    	part_medium = 0
    else:
    	part_medium = len(viol_medium)/medium_c 
    if weak_c == 0:
    	part_weak = 0
    else:
        part_weak = len(viol_weak)/weak_c

    return part_strong,s_ret_mean,s_ret_std, part_medium,m_ret_mean,m_ret_std, part_weak,w_ret_mean,w_ret_std     
    

def parselines(line_str):
    '''

    '''
    parts = line_str.split(')')
    resi_atom_i = parts[0].split('(')[1]
    resi_i = int(resi_atom_i.split()[0])
    resi_atom_j = parts[1].split('(')[1]
    resi_j = int(resi_atom_j.split()[0])
                    
    dist_info = parts[2].split()
    calc_dist = float(dist_info[0])
    given_dist_range = dist_info[1].split('..')
    given_dist_lo = float(given_dist_range[0])
    given_dist_hi = float(given_dist_range[1])
    diff = float(dist_info[2])

    viol_str = '%s - %s' %(resi_atom_i.lstrip(),resi_atom_j.lstrip())

    return viol_str, diff

def plotViols(strongfile, mediumfile, weakfile, part_strong, part_medium, part_weak, savefilename, strong_tol, medium_tol, weak_tol, tol):
    '''

    '''    
 
    data_normalizer = mp.colors.Normalize()
    color_map = mp.colors.LinearSegmentedColormap(
          "my_map",
          {
              "red": [(0, 1.0, 1.0),
                      (1.0, .5, .5)],
              "green": [(0, 0.5, 0.5),
                        (1.0, 0, 0)],
              "blue": [(0, 0.50, 0.5),
                        (1.0, 0, 0)]
          }
    )

    strong_str = 'strong - violation fraction:%s'%('{0:.2g}'.format(round(part_strong,2)))
    viols = [strong_str]
    viols_vals = [0]
    strong_line_pos = len(viols)

    with open(strongfile) as s:
        s_lines = s.readlines()
        for s_line in s_lines:
            if s_line.strip():
                viol_str, diff = parselines(s_line)
                viols.append(viol_str)
                viols_vals.append(diff)

    medium_str = 'medium - violation fraction:%s'%('{0:.2g}'.format(round(part_medium,2)))
    viols.append(medium_str)
    viols_vals.append(0)
    medium_line_pos = len(viols)

    with open(mediumfile) as m:
        m_lines = m.readlines()
        for m_line in m_lines:
            if m_line.strip():
                viol_str, diff = parselines(m_line)
                viols.append(viol_str)
                viols_vals.append(diff)

    weak_str = 'weak - violation fraction:%s'%('{0:.2g}'.format(round(part_weak,2)))
    viols.append(weak_str)
    viols_vals.append(0)       
    weak_line_pos = len(viols) 

    with open(weakfile) as w:
        w_lines = w.readlines()
        for w_line in w_lines:
            if w_line.strip():
                viol_str, diff = parselines(w_line)
                viols.append(viol_str)
                viols_vals.append(diff)


    y_pos = np.arange(len(viols))
    # Colorize the graph based on likeability:
    likeability_scores = np.array(viols_vals)

    fig = plt.figure(tight_layout=True) # need tight_layout to make everything fit
    ax = plt.subplot(111)
    ax.barh(y_pos, viols_vals, color=color_map(data_normalizer(likeability_scores)))
    ax.set_yticks(y_pos)
    #ax.set_yticklabels(viols)
    ax.set_yticklabels([])    
    for i, yi in enumerate(y_pos):
        ax.text(-0.0, yi, viols[i], horizontalalignment='left', verticalalignment='center')

    ax.axhline(y=strong_line_pos - 1.5, color='b', linestyle='--')
    ax.axhline(y=medium_line_pos - 1.5, color='b', linestyle='--')
    ax.axhline(y=weak_line_pos - 1.5, color='b', linestyle='--')

    #plt.tight_layout()
    #figure = plt.gcf()  # get current figure
    #fig.set_size_inches(2*32, 18) # set figure's size manually to your full screen (32x18)
    
    plt.savefig(savefilename, bbox_inches='tight', dpi=300) # bbox_inches removes extra white spaces

    #plt.show()

def plotViols_plotly(strongfile, mediumfile, weakfile, part_strong, part_medium, part_weak, savefilename, strong_tol, medium_tol, weak_tol, tol):
    '''

    '''    
 
    data_normalizer = mp.colors.Normalize()
    color_map = mp.colors.LinearSegmentedColormap(
          "my_map",
          {
              "red": [(0, 1.0, 1.0),
                      (1.0, .5, .5)],
              "green": [(0, 0.5, 0.5),
                        (1.0, 0, 0)],
              "blue": [(0, 0.50, 0.5),
                        (1.0, 0, 0)]
          }
    )

    strong_str = 'strong - violation fraction:%s'%('{0:.2g}'.format(round(part_strong,2)))
    viols = [strong_str]
    viols_vals = [0]
    strong_line_pos = len(viols)

    with open(strongfile) as s:
        s_lines = s.readlines()
        for s_line in s_lines:
            if s_line.strip():
                viol_str, diff = parselines(s_line)
                viols.append(viol_str)
                viols_vals.append(diff)

    medium_str = 'medium - violation fraction:%s'%('{0:.2g}'.format(round(part_medium,2)))
    viols.append(medium_str)
    viols_vals.append(0)
    medium_line_pos = len(viols)

    with open(mediumfile) as m:
        m_lines = m.readlines()
        for m_line in m_lines:
            if m_line.strip():
                viol_str, diff = parselines(m_line)
                viols.append(viol_str)
                viols_vals.append(diff)

    weak_str = 'weak - violation fraction:%s'%('{0:.2g}'.format(round(part_weak,2)))
    viols.append(weak_str)
    viols_vals.append(0)       
    weak_line_pos = len(viols) 

    with open(weakfile) as w:
        w_lines = w.readlines()
        for w_line in w_lines:
            if w_line.strip():
                viol_str, diff = parselines(w_line)
                viols.append(viol_str)
                viols_vals.append(diff)


    y_pos = np.arange(len(viols))
    # Colorize the graph based on likeability:
    likeability_scores = np.array(viols_vals)

    #fig = go.Figure(go.Bar(
    #        x=viols_vals,
    #        y=viols,
    #        orientation='h'))

    #fig.show()
    
    df = pd.DataFrame(
         {'viols':viols,
           'vals':viols_vals
         }
    )
    
    figheight = max(1024,math.ceil(3096/180 * len(viols))) #max(500,math.ceil(3096/180 * len(viols)))
    print(figheight)
    fig = px.bar(df, y='viols', x='vals',color='vals', orientation="h",color_continuous_scale='Bluered_r', height=figheight)
    fig.write_image(savefilename)
    #fig.show()


def writeViols(strongfile, mediumfile, weakfile, part_strong, part_medium, part_weak, savefilename, strong_tol, medium_tol, weak_tol, tol):
    '''

    '''    
 
    data_normalizer = mp.colors.Normalize()
    color_map = mp.colors.LinearSegmentedColormap(
          "my_map",
          {
              "red": [(0, 1.0, 1.0),
                      (1.0, .5, .5)],
              "green": [(0, 0.5, 0.5),
                        (1.0, 0, 0)],
              "blue": [(0, 0.50, 0.5),
                        (1.0, 0, 0)]
          }
    )

    strong_str = 'strong - violation fraction:%s'%('{0:.2g}'.format(round(part_strong,2)))
    viols = [strong_str]
    viols_vals = [0]
    strong_line_pos = len(viols)

    with open(strongfile) as s:
        s_lines = s.readlines()
        for s_line in s_lines:
            if s_line.strip():
                viol_str, diff = parselines(s_line)
                viols.append(viol_str)
                viols_vals.append(diff)

    medium_str = 'medium - violation fraction:%s'%('{0:.2g}'.format(round(part_medium,2)))
    viols.append(medium_str)
    viols_vals.append(0)
    medium_line_pos = len(viols)

    with open(mediumfile) as m:
        m_lines = m.readlines()
        for m_line in m_lines:
            if m_line.strip():
                viol_str, diff = parselines(m_line)
                viols.append(viol_str)
                viols_vals.append(diff)

    weak_str = 'weak - violation fraction:%s'%('{0:.2g}'.format(round(part_weak,2)))
    viols.append(weak_str)
    viols_vals.append(0)       
    weak_line_pos = len(viols) 

    with open(weakfile) as w:
        w_lines = w.readlines()
        for w_line in w_lines:
            if w_line.strip():
                viol_str, diff = parselines(w_line)
                viols.append(viol_str)
                viols_vals.append(diff)


    # simply write to a file
    data = np.column_stack([viols, viols_vals])
    np.savetxt(savefilename, data, fmt=['%s\t','%s'])

def main():
    xplorviolfile = sys.argv[1]
    opprefix = sys.argv[2]
    destfolder = os.path.dirname(xplorviolfile)

    defineTol()

    part_strong,s_ret_mean,s_ret_std, part_medium,m_ret_mean,m_ret_std, part_weak,w_ret_mean,w_ret_std = parseViolFile_v2(xplorviolfile, opprefix, destfolder)
    strong_file = destfolder+'/'+opprefix+'_strong.voils'
    medium_file = destfolder+'/'+opprefix+'_meidum.voils'
    weak_file = destfolder+'/'+opprefix+'_weak.voils'
    
    #savefilename = destfolder+'/'+opprefix+".png"
    #plotViols(strong_file, medium_file, weak_file, part_strong, part_medium, part_weak, savefilename, strong, medium, weak, tol)
    savefilename = destfolder+'/'+opprefix+".png"
    plotViols_plotly(strong_file, medium_file, weak_file, part_strong, part_medium, part_weak, savefilename, strong, medium, weak, tol)
    #savefilename = destfolder+'/'+opprefix+".viols.txt"
    #writeViols(strong_file, medium_file, weak_file, part_strong, part_medium, part_weak, savefilename, strong, medium, weak, tol)

if __name__ == "__main__":
    main()    