#!/usr/bin/env python

"""
   blah blah blah
"""
import numpy as np
import re
import argparse


def get_arguments():
    """ Get input and output filename from the commandline """
    parser = argparse.ArgumentParser(description='Merge itp file with'
        ' secondary structure exhibited by each residue in the output of'
        ' gmx do_dssp')
    parser.add_argument('posreitp',
                        help='posre.itp file')
    parser.add_argument('includeatomindices',
                        help='atom indices to be included from the gro file')
    parser.add_argument('outputfile',
                        help='output file for posre.itp')
    return parser.parse_args()

def makeListFromGro(filename):
    makelist = []
    with open(filename) as fp:
          for line in fp:
               noleadspc_line = line.lstrip()

               if re.match(r"^\d+.*$",noleadspc_line):
                    split_line = noleadspc_line.split()
                    makelist.append(int(split_line[0]))    
    return makelist

def makeListFromIncludeFile(filename):
    makelist = []
    with open(filename) as fp:
         for line in fp:
                   makelist.append(int(line))
    return makelist

def writeToFile(filename, listname):
    with open(filename,"w") as fid:
         for ele in listname:
                str = "%s\t1\t1000\t1000\t1000\n" %(ele)
                fid.write(str)
          

def main():
    """ Parse and plot the DSSP data """
    args = get_arguments()
    grolist = makeListFromGro(args.posreitp)
    grolist_nodup = list(set(grolist))
    includelist = makeListFromIncludeFile(args.includeatomindices)
    includelist_nodup = list(set(includelist))
    sortedlist  = sorted(grolist_nodup + includelist_nodup)
    writeToFile(args.outputfile, sortedlist)

if __name__ == '__main__':
    main()
    
