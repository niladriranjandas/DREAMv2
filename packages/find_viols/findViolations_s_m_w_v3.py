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
    #   ret_std = 0
    #   ret_mean = 0
    
    #if len(iplist) > 0 and ret_mean == 0 :
    ret_mean = np.mean(np.array(iplist),axis=0)
    ret_std  = np.std(np.array(iplist),axis=0)

    return ret_std, ret_mean

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

    viol_str = '%s,%s,%f' %(resi_atom_i.lstrip(),resi_atom_j.lstrip(),diff)

    return resi_atom_i.lstrip(), resi_atom_j.lstrip(), diff

def makeSortedDataframe(violfile):
    '''
    
    '''
    
    df = pd.DataFrame(columns = ['name_i', 'name_j', 'viol'])
    
    with open(violfile) as v:
        v_lines = v.readlines()
        for v_line in v_lines:
            if v_line.strip():
                resi_atom_i, resi_atom_j, violval = parselines(v_line)
                df = df.append({'name_i' : resi_atom_i, 'name_j' : resi_atom_j, 'viol' : violval}, ignore_index = True)
                

    sorted_df = df.sort_values(by=['viol'], ascending=True)
    return sorted_df, df    

def makeDataFrameForPlot(df_strong, part_strong, df_medium, part_medium, df_weak, part_weak):
    '''
    
    '''
    
    viols = []
    viols_vals = []
    
    df_strong['combined'] = df_strong['name_i'] + '-' + df_strong['name_j']
    df_medium['combined'] = df_medium['name_i'] + '-' + df_medium['name_j']
    df_weak['combined']   = df_weak['name_i']   + '-' + df_weak['name_j']
    
    strong_str = 'strong - violation fraction:%s'%('{0:.2g}'.format(round(part_strong,2)))
    viols.append(strong_str)
    viols_vals.append(0)
    viols.extend(df_strong['combined'].to_list())
    viols_vals.extend(df_strong['viol'].to_list())
    
    medium_str = 'medium - violation fraction:%s'%('{0:.2g}'.format(round(part_medium,2)))
    viols.append(medium_str)
    viols_vals.append(0)
    viols.extend(df_medium['combined'].to_list())
    viols_vals.extend(df_medium['viol'].to_list())
    
    weak_str = 'weak - violation fraction:%s'%('{0:.2g}'.format(round(part_weak,2)))
    viols.append(weak_str)
    viols_vals.append(0)
    viols.extend(df_weak['combined'].to_list())
    viols_vals.extend(df_weak['viol'].to_list())     
    
    df = pd.DataFrame({'viols':viols,'vals':viols_vals})
        
    return df


def plotViols_plotly(df, savefilename):
    '''
    
    '''
    ylabel = 'cases'
    xlabel =  'violations(\u212B deviation)'
    
    figheight = max(1024,math.ceil(3096/180 * df.shape[0])) #max(500,math.ceil(3096/180 * len(viols)))
    print(figheight)
    fig = px.bar(df, y='viols', x='vals',color='vals', orientation="h",color_continuous_scale='Bluered_r', height=figheight, labels=dict(vals=xlabel, viols=ylabel))
    fig.write_image(savefilename)

def plotViols_plotly_v2(df, savefilename):
    '''
    
    '''
    ylabel = 'cases'
    xlabel =  'violations(\u212B deviation)'    
    
    figheight = max(1024,math.ceil(3096/180 * df.shape[0])) #max(500,math.ceil(3096/180 * len(viols)))
    print(figheight)
    fig = px.bar(df, y='viols', x='vals',orientation="h",height=figheight, labels=dict(vals=xlabel, viols=ylabel))
    fig.update_xaxes(title_font=dict(size=32), tickfont=dict(size=18)) #, family='Courier', color='crimson'))
    fig.update_yaxes(title_font=dict(size=32), tickfont=dict(size=18)) #, family='Courier', color='crimson'))
    fig.write_image(savefilename)


def main():
    xplorviolfile = sys.argv[1]
    opprefix = sys.argv[2]
    destfolder = os.path.dirname(xplorviolfile)

    defineTol()

    #part_strong,s_ret_mean,s_ret_std, part_medium,m_ret_mean,m_ret_std, part_weak,w_ret_mean,w_ret_std = parseViolFile_v2(xplorviolfile, opprefix)
    part_strong,s_ret_mean,s_ret_std, part_medium,m_ret_mean,m_ret_std, part_weak,w_ret_mean,w_ret_std = parseViolFile_v2(xplorviolfile, opprefix, destfolder)
    strong_file = destfolder+'/'+opprefix+'_strong.voils'
    medium_file = destfolder+'/'+opprefix+'_meidum.voils'
    weak_file = destfolder+'/'+opprefix+'_weak.voils'

    df_strong_sorted , df_strong = makeSortedDataframe(strong_file)
    df_medium_sorted , df_medium = makeSortedDataframe(medium_file)
    df_weak_sorted   , df_weak   = makeSortedDataframe(weak_file)

    df_plot = makeDataFrameForPlot(df_strong_sorted, part_strong, df_medium_sorted, part_medium, df_weak_sorted, part_weak)    

    #savefilename = destfolder+'/'+opprefix+".png"
    #plotViols_plotly(df_plot, savefilename)

    savefilename = destfolder+'/'+opprefix+".png"
    plotViols_plotly_v2(df_plot, savefilename)

if __name__ == "__main__":
    main()   
