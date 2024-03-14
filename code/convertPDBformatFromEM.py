from __future__ import print_function
import sys
import pandas as pd

def mapCyanaToXplor(residue,cyanaatom, mapfile_dataframe):
    '''

    '''
    tmp = mapfile_dataframe.loc[ (mapfile_dataframe['residue'] == residue) & (mapfile_dataframe['gro'] == cyanaatom) ]
    tmp_ = tmp['xplor'].to_string().split()
    return tmp_[1]


def replaceLine(lineToConvert, mapfile_dataframe):
    '''

    '''
    items = lineToConvert.split()
    atom_name = items[2]
    resi_name = items[3]

    target_atom = mapCyanaToXplor(resi_name, atom_name, mapfile_dataframe)

    source='%s'%(atom_name)
    if (target_atom == 'HN'):
    	source='%s'%(atom_name)
    else:
        source='%s'%(atom_name)
    target='%s'%(target_atom)

    return lineToConvert.replace(source,target)

def readLinesAndConvert(pdbfilename, mappedFile, oppdb):
    '''    

    '''
    data_map = pd.read_csv(mappedFile)

    with open(pdbfilename) as fp, open(oppdb,'w') as wp :
         lines = fp.readlines()

         for eachline in lines:
             if eachline.startswith('ATOM'):
                  newline = replaceLine(eachline, data_map) 
                  #print('%s'%(newline),end =" ")
                  wp.write(newline)
             else:
                  #print('%s'%(eachline),end =" ")
                  wp.write(eachline)


def main():
    readLinesAndConvert(sys.argv[1], sys.argv[2], sys.argv[3])

if __name__ == "__main__":
    main()                           
