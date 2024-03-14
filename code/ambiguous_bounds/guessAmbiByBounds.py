import sys
import pandas as pd

###
#./masterPreCheckUpl.sh -suf 2m4k_now -upl 2m4k/2m4k_5perc_1_added.upl -seq 2m4k/2m4k.seq -hb '' -aco 2m4k/2m4k_concat_dihed.aco

def defGloablVars(TOLSTRONG_ = 10, TOLMEDIUM_ = 10, TOLLOW_ = 5):
    
    global STRONG
    global MEDIUM
    global LOW   

    global TOLSTRONG
    global TOLMEDIUM
    global TOLLOW   

    STRONG = 2.7
    MEDIUM = 3.5
    LOW    = 6.0

    TOLSTRONG = TOLSTRONG_  # for strongCheck()
    TOLMEDIUM = TOLMEDIUM_  # for mediumCheck()
    TOLLOW    = TOLLOW_     # for lowCheck()

def parseUpls(uplfile):
    '''
    parse the upl file in cyana format
    and give dicts for strong medium low
    '''

    dict_strong = {}
    dict_medium = {}
    dict_low    = {}

    count = 0
    with open(uplfile) as up:
        lines = up.readlines()
        for line in lines:
            #count = count+1
            parts = line.split()

            resi_i      = int(parts[0])
            resi_i_nm   = parts[1]
            resi_i_atm  = parts[2]

            resi_j      = int(parts[3])
            resi_j_nm   = parts[4]
            resi_j_atm  = parts[5]

            upl_dist_ij = float(parts[6])

            if upl_dist_ij > 0:
                 count = count +1   # this is to maintain cosistency with raw_upl after preprocess.m
                 if upl_dist_ij >0 and upl_dist_ij <= STRONG:
                      dict_strong[count] = [resi_i,resi_i_nm,resi_i_atm, resi_j,resi_j_nm,resi_j_atm,upl_dist_ij]
                 elif upl_dist_ij >STRONG and upl_dist_ij <= MEDIUM:
                      dict_medium[count] = [resi_i,resi_i_nm,resi_i_atm, resi_j,resi_j_nm,resi_j_atm,upl_dist_ij]
                 elif upl_dist_ij >MEDIUM and upl_dist_ij <= LOW:
                      dict_low[count] = [resi_i,resi_i_nm,resi_i_atm, resi_j,resi_j_nm,resi_j_atm,upl_dist_ij]

    
    return dict_strong, dict_medium, dict_low

def strongCheck(dict_strong, dict_medium):
    '''
    check for strong with strong
    check for strong with medium
    '''

    dict_ss = {}
    dict_sm = {}

    for key in dict_strong:
        list_key  = dict_strong[key]
        key_i     = list_key[0]
        key_i_nm  = list_key[1]
        key_i_atm = list_key[2]
        key_j     = list_key[3]
        key_j_nm  = list_key[4]
        key_j_atm = list_key[5]

        count_ss = 0        
        for otherkey in set(dict_strong) - {key}:
            other_list_key  = dict_strong[otherkey]
            other_key_i     = other_list_key[0]
            other_key_i_nm  = other_list_key[1]
            other_key_i_atm = other_list_key[2]
            other_key_j     = other_list_key[3]
            other_key_j_nm  = other_list_key[4]
            other_key_j_atm = other_list_key[5]

            dist_i_1 = abs(key_i - other_key_i)
            dist_i_2 = abs(key_i - other_key_j)
            dist_i   = min(dist_i_1, dist_i_2)

            dist_j_1 = abs(key_j - other_key_i)
            dist_j_2 = abs(key_j - other_key_j)
            dist_j   = min(dist_j_1, dist_j_2)

            if max(dist_i, dist_j) <= TOLSTRONG:
                  count_ss = count_ss+1

        dict_ss[key] = count_ss

        count_sm = 0        
        for otherkey in dict_medium:
            other_list_key  = dict_medium[otherkey]
            other_key_i     = other_list_key[0]
            other_key_i_nm  = other_list_key[1]
            other_key_i_atm = other_list_key[2]
            other_key_j     = other_list_key[3]
            other_key_j_nm  = other_list_key[4]
            other_key_j_atm = other_list_key[5]

            dist_i_1 = abs(key_i - other_key_i)
            dist_i_2 = abs(key_i - other_key_j)
            dist_i   = min(dist_i_1, dist_i_2)

            dist_j_1 = abs(key_j - other_key_i)
            dist_j_2 = abs(key_j - other_key_j)
            dist_j   = min(dist_j_1, dist_j_2)

            if max(dist_i, dist_j) <= TOLSTRONG:
                  count_sm = count_sm+1

        dict_sm[key] = count_sm

    return dict_ss, dict_sm

def mediumCheck(dict_strong, dict_medium, dict_low):
    '''
    check for medium with strong
    check for medium with medium
    '''

    dict_mm = {}
    dict_ms = {}
    dict_ml = {}

    for key in dict_medium:
        list_key  = dict_medium[key]
        key_i     = list_key[0]
        key_i_nm  = list_key[1]
        key_i_atm = list_key[2]
        key_j     = list_key[3]
        key_j_nm  = list_key[4]
        key_j_atm = list_key[5]

        count_ss = 0        
        for otherkey in set(dict_medium) - {key}:
            other_list_key  = dict_medium[otherkey]
            other_key_i     = other_list_key[0]
            other_key_i_nm  = other_list_key[1]
            other_key_i_atm = other_list_key[2]
            other_key_j     = other_list_key[3]
            other_key_j_nm  = other_list_key[4]
            other_key_j_atm = other_list_key[5]

            dist_i_1 = abs(key_i - other_key_i)
            dist_i_2 = abs(key_i - other_key_j)
            dist_i   = min(dist_i_1, dist_i_2)

            dist_j_1 = abs(key_j - other_key_i)
            dist_j_2 = abs(key_j - other_key_j)
            dist_j   = min(dist_j_1, dist_j_2)

            if max(dist_i, dist_j) <= TOLMEDIUM:
                  count_ss = count_ss+1

        dict_mm[key] = count_ss

        count_ms = 0        
        for otherkey in dict_strong:
            other_list_key  = dict_strong[otherkey]
            other_key_i     = other_list_key[0]
            other_key_i_nm  = other_list_key[1]
            other_key_i_atm = other_list_key[2]
            other_key_j     = other_list_key[3]
            other_key_j_nm  = other_list_key[4]
            other_key_j_atm = other_list_key[5]

            dist_i_1 = abs(key_i - other_key_i)
            dist_i_2 = abs(key_i - other_key_j)
            dist_i   = min(dist_i_1, dist_i_2)

            dist_j_1 = abs(key_j - other_key_i)
            dist_j_2 = abs(key_j - other_key_j)
            dist_j   = min(dist_j_1, dist_j_2)

            if max(dist_i, dist_j) <= TOLMEDIUM:
                  count_ms = count_ms+1

        dict_ms[key] = count_ms

        count_ml = 0        
        for otherkey in dict_low:
            other_list_key  = dict_low[otherkey]
            other_key_i     = other_list_key[0]
            other_key_i_nm  = other_list_key[1]
            other_key_i_atm = other_list_key[2]
            other_key_j     = other_list_key[3]
            other_key_j_nm  = other_list_key[4]
            other_key_j_atm = other_list_key[5]

            dist_i_1 = abs(key_i - other_key_i)
            dist_i_2 = abs(key_i - other_key_j)
            dist_i   = min(dist_i_1, dist_i_2)

            dist_j_1 = abs(key_j - other_key_i)
            dist_j_2 = abs(key_j - other_key_j)
            dist_j   = min(dist_j_1, dist_j_2)

            if max(dist_i, dist_j) <= TOLMEDIUM:
                  count_ml = count_ml+1

        dict_ml[key] = count_ml


    return dict_mm, dict_ms, dict_ml

def lowCheck(dict_medium, dict_low):
    '''
    check for strong with strong
    check for strong with medium
    '''

    dict_ll = {}
    dict_lm = {}

    for key in dict_low:
        list_key  = dict_low[key]
        key_i     = list_key[0]
        key_i_nm  = list_key[1]
        key_i_atm = list_key[2]
        key_j     = list_key[3]
        key_j_nm  = list_key[4]
        key_j_atm = list_key[5]

        count_ll = 0        
        for otherkey in set(dict_low) - {key}:
            other_list_key  = dict_low[otherkey]
            other_key_i     = other_list_key[0]
            other_key_i_nm  = other_list_key[1]
            other_key_i_atm = other_list_key[2]
            other_key_j     = other_list_key[3]
            other_key_j_nm  = other_list_key[4]
            other_key_j_atm = other_list_key[5]

            dist_i_1 = abs(key_i - other_key_i)
            dist_i_2 = abs(key_i - other_key_j)
            dist_i   = min(dist_i_1, dist_i_2)

            dist_j_1 = abs(key_j - other_key_i)
            dist_j_2 = abs(key_j - other_key_j)
            dist_j   = min(dist_j_1, dist_j_2)

            if max(dist_i, dist_j) <= TOLLOW:
                  count_ll = count_ll+1

        dict_ll[key] = count_ll

        count_lm = 0        
        for otherkey in dict_medium:
            other_list_key  = dict_medium[otherkey]
            other_key_i     = other_list_key[0]
            other_key_i_nm  = other_list_key[1]
            other_key_i_atm = other_list_key[2]
            other_key_j     = other_list_key[3]
            other_key_j_nm  = other_list_key[4]
            other_key_j_atm = other_list_key[5]

            dist_i_1 = abs(key_i - other_key_i)
            dist_i_2 = abs(key_i - other_key_j)
            dist_i   = min(dist_i_1, dist_i_2)

            dist_j_1 = abs(key_j - other_key_i)
            dist_j_2 = abs(key_j - other_key_j)
            dist_j   = min(dist_j_1, dist_j_2)

            if max(dist_i, dist_j) <= TOLLOW:
                  count_lm = count_lm+1

        dict_lm[key] = count_lm

    return dict_ll, dict_lm


def flagTheBounds(upl_dict, oppfilename):
    '''

    '''

    if upl_dict:
        df_upl = pd.DataFrame(list(upl_dict.items()))
        df_upl.rename(columns={df_upl.columns[0]: "uplindex"}, inplace=True)
        df_upl.rename(columns={df_upl.columns[1]: "count"}, inplace=True)

        df_upl_mean = df_upl.iloc[:,1].mean()
        df_upl_std  = df_upl.iloc[:,1].std()

        flag_tol = df_upl_mean - df_upl_std

        df_upl['flag'] = df_upl.iloc[:,1]<flag_tol

        df_upl.to_csv(oppfilename)
  
        return df_upl
    else:
        return False

def driverCode(uplfile, prefix='parsed_upls'):
    '''
    parse upl
    get strong-srong
        strong-medium

        medium-strong
        medium-medium
        medium-low

        low-medium
        low-low

    flag the entries as true/false
    combine with or for each category
    '''

    ### parse and do the counting    
    dict_strong, dict_medium, dict_low = parseUpls(uplfile)    
    df_strong = pd.DataFrame(list(dict_strong.items()))
    df_medium = pd.DataFrame(list(dict_medium.items()))
    df_low    = pd.DataFrame(list(dict_low.items()))

    dict_ss, dict_sm          = strongCheck(dict_strong, dict_medium)
    dict_mm, dict_ms, dict_ml = mediumCheck(dict_strong, dict_medium, dict_low)
    dict_ll, dict_lm          = lowCheck(dict_medium, dict_low)

    name_ss = '%s_ss.%d.csv'%(prefix,TOLMEDIUM)
    name_sm = '%s_sm.%d.csv'%(prefix,TOLMEDIUM)

    name_mm = '%s_mm.%d.csv'%(prefix,TOLMEDIUM)
    name_ms = '%s_ms.%d.csv'%(prefix,TOLMEDIUM)
    name_ml = '%s_ml.%d.csv'%(prefix,TOLMEDIUM)

    name_ll = '%s_ll.%d.csv'%(prefix,TOLLOW)
    name_lm = '%s_lm.%d.csv'%(prefix,TOLLOW)

    ### flag the entries as true/false 
    if dict_ss:
       df_ss_flag = flagTheBounds(dict_ss, name_ss)
    else:
        df_ss_flag = ''
    if dict_sm:
        df_sm_flag = flagTheBounds(dict_sm, name_sm)
    else:
        df_sm_flag = ''

    df_mm_flag = flagTheBounds(dict_mm, name_mm)
    df_ms_flag = flagTheBounds(dict_ms, name_ms)
    df_ml_flag = flagTheBounds(dict_ml, name_ml)
    
    df_ll_flag = flagTheBounds(dict_ll, name_ll)
    df_lm_flag = flagTheBounds(dict_lm, name_lm)

    if not df_strong.empty:
        df_strong.rename(columns={df_strong.columns[0]: "uplindex"}, inplace=True)
        df_strong.rename(columns={df_strong.columns[1]: "uplinfo"}, inplace=True)
        df_strong['flag'] = df_ss_flag['flag'] | df_sm_flag['flag']

    df_medium.rename(columns={df_medium.columns[0]: "uplindex"}, inplace=True)
    df_medium.rename(columns={df_medium.columns[1]: "uplinfo"}, inplace=True)
    df_medium['flag'] = df_mm_flag['flag'] | df_ms_flag['flag'] | df_ml_flag['flag']

    df_low.rename(columns={df_low.columns[0]: "uplindex"}, inplace=True)
    df_low.rename(columns={df_low.columns[1]: "uplinfo"}, inplace=True)
    df_low['flag'] = df_lm_flag['flag'] | df_ll_flag['flag']
 
    return df_strong, df_medium, df_low 

def main():

    uplfile = sys.argv[1]
    suffix  = sys.argv[2]
   
    if len(sys.argv) > 3:
       TOLSTRONG_ = float(sys.argv[3])
       TOLMEDIUM_ = float(sys.argv[4])
       TOLLOW_    = float(sys.argv[5])    
       defGloablVars(TOLSTRONG_, TOLMEDIUM_, TOLLOW_)
    else:
       defGloablVars() 

    df_strong, df_medium, df_low = driverCode(uplfile, suffix)
    if not df_strong.empty:
       df_strong.to_csv('%s_strong_flag.csv'%(suffix))
    df_medium.to_csv('%s_medium_flag.csv'%(suffix))
    df_low.to_csv('%s_low_flag.csv'%(suffix))

if __name__ == "__main__":
    main()
