import pandas as pd 
import argparse

def findDuplidates(uplfile):
    '''

    '''

    #df = pd.read_table(uplfile, header=None)
    df = pd.read_csv(uplfile, header=None)
    df.index = df.index+1
    
    df_grp = df.groupby([0,1,2,3,4,5]).count()        
    df_grp.reset_index().to_csv('%s.groupby.csv'%(uplfile), header=True, index=False)

    df_grp_csv = pd.read_csv('%s.groupby.csv'%(uplfile))

    df_dups = df_grp_csv.loc[df_grp_csv['6'] > 1]
    
    for i in range(0,len(df_dups)):
        row = df_dups.iloc[i]

        resi_i_no = row[0]
        resi_i_nm = row[1]
        resi_i_at = row[2]
        resi_j_no = row[3]
        resi_j_nm = row[4]
        resi_j_at = row[5]

        df_tmp = df.loc[(df[0] == resi_i_no) & (df[1] == resi_i_nm) & (df[2] == resi_i_at) & (df[3] == resi_j_no) & (df[4] == resi_j_nm) & (df[5] == resi_j_at)]
        print(df_tmp)
        print('\n')

    # find entries where both i and j are same
    df_same_i_j = df.loc[ (df[0] == df[3]) & (df[1] == df[4]) & (df[2] == df[5]) ]
    print(df_same_i_j)

def main():
    my_parser = argparse.ArgumentParser(fromfile_prefix_chars='@')

    my_parser.add_argument('upl',
                       help='uplfile separated by tabs')

    my_parser.add_argument('-v',
                       '--verbose',
                       action='store_true',
                       help='an optional argument')    

    args = my_parser.parse_args()

    findDuplidates(args.upl)
    

if __name__ == "__main__":
    main()

        



