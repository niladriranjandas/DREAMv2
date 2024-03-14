import pandas as pd
import argparse

# e.g.: python parseAllFiles.py ../2m4k_ambi5r1_noe_medium.csv  2m4k_ambi5r1_noe_1_medium_flag.csv,2m4k_ambi5r1_noe_2_medium_flag.csv,2m4k_ambi5r1_noe_3_medium_flag.csv '2m4k_ambi5r1_noe_medium.csv'


def parseFiles(densityfiles, uplnearbyfiles, oppfile):
    '''

    '''

    files_density   = densityfiles.split(',')
    files_uplnearby = uplnearbyfiles.split(',')
    
    # the density files
    count = 0
    density_dataframe = pd.DataFrame(columns=['uplindex','resii','resij','density_together'])

    for files in files_density:
        data = pd.read_csv(files)
        count = count+1

        if count == 1:
            density_dataframe['uplindex'] = data['uplindex']
            density_dataframe['resii']    = data['resii']
            density_dataframe['resij']    = data['resij']
            density_dataframe['distij']   = data['distij']

            density_dataframe['density_together'] = data['density_together']
        else:
            density_tmp = density_dataframe.merge(data, how="inner", on=['uplindex'])
            density_dataframe['density_together'] = density_tmp['density_together_x'] | density_tmp['density_together_y']
            
           #density_dataframe['densityi_%d'%(count)]  = data['densityi']
           #density_dataframe['densityj_%d'%(count)]  = data['densityj']
           #density_dataframe['densityij_%d'%(count)] = data['densityij']
           #density_dataframe['density_together_%d'%(count)] = data['density_together']


       # the uplnearby files
    count = 0
    uplnearby_dataframe = pd.DataFrame(columns=['uplindex','uplinfo','flag'])

    for files in files_uplnearby:
        data = pd.read_csv(files)
        count = count+1

        if count == 1:
            #print(data)
            uplnearby_dataframe['uplindex'] = data['uplindex']
            uplnearby_dataframe['uplinfo']  = data['uplinfo']

            uplnearby_dataframe['flag'] = data['flag']
        else:
            uplnearby_tmp = uplnearby_dataframe.merge(data, how="inner", on=['uplindex'])
            uplnearby_dataframe['flag'] = uplnearby_tmp['flag_x'] | uplnearby_tmp['flag_y']

    
    all_together = density_dataframe.merge(uplnearby_dataframe, how='inner', on=['uplindex'])
    all_together['ambiguous'] = all_together['density_together'] & all_together['flag']

    all_together.to_csv(oppfile, index = False)

    return density_dataframe, uplnearby_dataframe, all_together

def main():
    my_parser = argparse.ArgumentParser(fromfile_prefix_chars='@')

    my_parser.add_argument('density',
                       help='density files separated by comma')

    my_parser.add_argument('uplnearby',
                       help='upls nearby files separated by comma')

    my_parser.add_argument('opp',
                       help='ooutput file')

    my_parser.add_argument('-v',
                       '--verbose',
                       action='store_true',
                       help='an optional argument')    

    args = my_parser.parse_args()

    density_dataframe, uplnearby_dataframe, all_together = parseFiles(args.density, args.uplnearby, args.opp)

if __name__ == "__main__":
    main()
  

