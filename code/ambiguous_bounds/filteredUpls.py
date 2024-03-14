import pandas as pd
import numpy as np
import argparse

def filterOutUplEntries(markeduplorg, stronguplflagged, mediumuplflagged, weakflagged):
    '''
                           org_upl-----   0         1            
    confusion_matrix =                 0  true-ve   false+ve
                                       1  false+ve  true+ve
                                       |
                                       |
                               predicted        

                       true-ve  : false predicted as false (non-ambiguous recognised)
                       false+ve : true  predicted as false (ambiguous wrongly flagged as not ambiguous)
                       false+ve : false predicted as true  (non-ambiguous predicted as ambiguous)
                       true+ve  : true  predicted as true  (ambiguous predicted as ambiguous)

    e.g.: import filteredUpls as fp
          confusion_matrix, confusion_matrix_part, true_ambicount, my_ambicount = fp.filterOutUplEntries('../2m4k_5perc_1_added.upl.raw_upl.marked.csv','','/data2/nmr/our_algo_final/protein/2m4k_ambi5r1_noe/ambiguous_bounds/2m4k_ambi5r1_noe_medium_overall.csv','/data2/nmr/our_algo_final/protein/2m4k_ambi5r1_noe/ambiguous_bounds/2m4k_ambi5r1_noe_weak_overall.csv')
          fp.writeNonAmbiguousBounds('../2m4k_5perc_1_added.upl.raw_upl.marked.csv','','/data2/nmr/our_algo_final/protein/2m4k_ambi5r1_noe/ambiguous_bounds/2m4k_ambi5r1_noe_medium_overall.csv','/data2/nmr/our_algo_final/protein/2m4k_ambi5r1_noe/ambiguous_bounds/2m4k_ambi5r1_noe_weak_overall.csv','2m4k_ambi5r1_noe_final.upl')

          or

          python filteredUpls.py '../2m4k_5perc_1_added.upl.raw_upl.marked.csv' '' '/data2/nmr/our_algo_final/protein/2m4k_ambi5r1_noe/ambiguous_bounds/2m4k_ambi5r1_noe_medium_overall.csv' '/data2/nmr/our_algo_final/protein/2m4k_ambi5r1_noe/ambiguous_bounds/2m4k_ambi5r1_noe_weak_overall.csv' '2m4k_ambi5r1_noe_final.upl'
    '''
    
    confusion_matrix = np.zeros((2,2), dtype=np.int64)
    confusion_matrix_part = np.zeros((2,2), dtype=np.float64)
    true_ambicount = 0
    my_ambicount   = 0

    df_org = pd.read_table(markeduplorg, header=None)
    #df_org = pd.read_csv(markeduplorg, header=None)
    df_org.index = np.arange(1, len(df_org)+1)

    du_flagged_upls = pd.DataFrame(columns=['uplindex','resii','resij','density_together','distij','uplinfo','flag','ambiguous'])
    if stronguplflagged:
       du_strong = pd.read_csv(stronguplflagged)
       du_flagged_upls = du_flagged_upls.append(du_strong)

    if mediumuplflagged:       
       du_medium = pd.read_csv(mediumuplflagged)
       du_flagged_upls = du_flagged_upls.append(du_medium)

    du_weak   = pd.read_csv(weakflagged)
    du_flagged_upls = du_flagged_upls.append(du_weak)

    du_flagged_upls.set_index('uplindex', inplace=True)

    #print(df_org)
    #print(du_flagged_upls)

    # ----- iterate through df_org and check with du_flagged_upls ----- #
    total_count_org = 0
    for index, row in df_org.iterrows():
        #print(row)
        flag_org = row[7]  # note that this '7' is the default column name starting from 0

        if index in du_flagged_upls.index:
              total_count_org = total_count_org+1

              row_predicted  = du_flagged_upls.loc[index]
              predicted_flag = row_predicted['ambiguous']

              # ----- populate the confusion matrix ----- #
              if flag_org == False and predicted_flag == False:
                    confusion_matrix[0,0] = confusion_matrix[0,0]+1

              if flag_org == True and predicted_flag == False:
                    confusion_matrix[0,1] = confusion_matrix[0,1]+1

              if flag_org == False and predicted_flag == True:
                    confusion_matrix[1,0] = confusion_matrix[1,0]+1

              if flag_org == True and predicted_flag == True:
                    confusion_matrix[1,1] = confusion_matrix[1,1]+1

              # ----- overall count ---- #
              if flag_org == True:
                    true_ambicount = true_ambicount +1

              if predicted_flag == True:
                    my_ambicount = my_ambicount+1
    

    true_nonambicount = total_count_org - true_ambicount
    confusion_matrix_part[0,0] = confusion_matrix[0,0]/true_nonambicount
    confusion_matrix_part[0,1] = confusion_matrix[0,1]/true_ambicount
    confusion_matrix_part[1,0] = confusion_matrix[1,0]/true_nonambicount
    confusion_matrix_part[1,1] = confusion_matrix[1,1]/true_ambicount

    return confusion_matrix, confusion_matrix_part, true_ambicount, my_ambicount


def writeFalseNegatives(markeduplorg, stronguplflagged, mediumuplflagged, weakflagged): #, oppambis):
    '''

    '''

    df_org = pd.read_table(markeduplorg, header=None)
    #df_org = pd.read_csv(markeduplorg, header=None)
    df_org.index = np.arange(1, len(df_org)+1)

    du_flagged_upls = pd.DataFrame(columns=['uplindex','resii','resij','density_together','distij','uplinfo','flag','ambiguous'])
    if stronguplflagged:
       du_strong = pd.read_csv(stronguplflagged)
       du_flagged_upls = du_flagged_upls.append(du_strong)

    if mediumuplflagged:       
       du_medium = pd.read_csv(mediumuplflagged)
       du_flagged_upls = du_flagged_upls.append(du_medium)

    du_weak   = pd.read_csv(weakflagged)
    du_flagged_upls = du_flagged_upls.append(du_weak)

    du_flagged_upls.set_index('uplindex', inplace=True)

    #print(df_org)
    #print(du_flagged_upls)

    # ----- iterate through df_org and check with du_flagged_upls ----- #
    total_count_org = 0
    for index, row in df_org.iterrows():
        #print(row)
        flag_org = row[7]  # note that this '7' is the default column name starting from 0

        if index in du_flagged_upls.index:
              total_count_org = total_count_org+1

              row_predicted  = du_flagged_upls.loc[index]
              predicted_flag = row_predicted['ambiguous']

              # ----- populate the ambis which got through i.e. (0,1) of confusion matrix ----- #

              if flag_org == True and predicted_flag == False:
                 #print(row)
                 print('%d\t%s\t%s\t%d\t%s\t%s\t%3.2f'%(int(row[0]),row[1].strip(),row[2].strip(),int(row[3]),row[4].strip(),row[5].strip(),float(row[6])))
                 #if count == 1:
                 #         op.write('%d\t%s\t%s\t%d\t%s\t%s\t%3.2f'%(int(uplinfo[0]),uplinfo[1].strip(),uplinfo[2].strip(),int(uplinfo[3]),uplinfo[4].strip(),uplinfo[5].strip(),float(uplinfo[6])))
                 #else:
                 #         op.write('\n%d\t%s\t%s\t%d\t%s\t%s\t%3.2f'%(int(uplinfo[0]),uplinfo[1].strip(),uplinfo[2].strip(),int(uplinfo[3]),uplinfo[4].strip(),uplinfo[5].strip(),float(uplinfo[6])))


    

def writeNonAmbiguousBounds(markeduplorg, stronguplflagged, mediumuplflagged, weakflagged, oppupl):
    '''

    '''
    du_flagged_upls = pd.DataFrame(columns=['uplindex','resii','resij','density_together','distij','uplinfo','flag','ambiguous'])
    if stronguplflagged:
       du_strong = pd.read_csv(stronguplflagged)
       du_flagged_upls = du_flagged_upls.append(du_strong)

    if mediumuplflagged:
       du_medium = pd.read_csv(mediumuplflagged)
       du_flagged_upls = du_flagged_upls.append(du_medium)

    du_weak   = pd.read_csv(weakflagged)
    du_flagged_upls = du_flagged_upls.append(du_weak)

    du_flagged_upls.set_index('uplindex', inplace=True)
    du_flagged_upls.sort_index(ascending=True, inplace=True)

    df_org = pd.read_table(markeduplorg, header=None)
    df_org.index = np.arange(1, len(df_org)+1)

    # ----- iterate through df_org and check with du_flagged_upls ----- #
    count = 0
    with open(oppupl, 'w') as op:
        for index, row in df_org.iterrows():
             #flag_org = row[7]   # note that this '7' is the default column name starting from 0
             count = count+1
             if index in du_flagged_upls.index:
                 row_predicted  = du_flagged_upls.loc[index]
                 predicted_flag = row_predicted['ambiguous']

                 if predicted_flag == False:
                     uplinfo = row_predicted['uplinfo'].lstrip('[').rstrip(']').split(',')
                     if count == 1:
                          op.write('%d\t%s\t%s\t%d\t%s\t%s\t%3.2f'%(int(uplinfo[0]),uplinfo[1].strip(),uplinfo[2].strip(),int(uplinfo[3]),uplinfo[4].strip(),uplinfo[5].strip(),float(uplinfo[6])))
                     else:
                          op.write('\n%d\t%s\t%s\t%d\t%s\t%s\t%3.2f'%(int(uplinfo[0]),uplinfo[1].strip(),uplinfo[2].strip(),int(uplinfo[3]),uplinfo[4].strip(),uplinfo[5].strip(),float(uplinfo[6])))
             else:  # this means d_ij > 6 write it as it is
                 if count == 1:
                     op.write('%d\t%s\t%s\t%d\t%s\t%s\t%3.2f'%(row[0],row[1],row[2],row[3],row[4],row[5],row[6]))
                 else:
                     op.write('\n%d\t%s\t%s\t%d\t%s\t%s\t%3.2f'%(row[0],row[1],row[2],row[3],row[4],row[5],row[6]))


def main():
    my_parser = argparse.ArgumentParser(fromfile_prefix_chars='@')

    my_parser.add_argument('orguplmarked',
                       help='density files separated by comma')

    my_parser.add_argument('stronguplflagged',
                       help='strong upl list flagged as true/false')

    my_parser.add_argument('mediumuplflagged',
                       help='strong upl list flagged as true/false')

    my_parser.add_argument('weakuplflagged',
                       help='strong upl list flagged as true/false')    

    my_parser.add_argument('opp',
                       help='ooutput file')

    my_parser.add_argument('-v',
                       '--verbose',
                       action='store_true',
                       help='an optional argument')    

    args = my_parser.parse_args()

    writeNonAmbiguousBounds(args.orguplmarked, args.stronguplflagged, args.mediumuplflagged, args.weakuplflagged, args.opp)
    

if __name__ == "__main__":
    main()
