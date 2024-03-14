#!/usr/bin/env python

"""
abc
e.g. python readSpreadSheet.py 1xxe_rough/strongbounds_gt1_strand_pdb_not_1xxe_v2.xlsx 1xxe_resi124_139_medium -o 1xxe_rough/dist_1xxe_124_139
"""

import numpy as np
import argparse
import xlrd

def get_arguments():
    """ Get input and output filename from the commandline """
    parser = argparse.ArgumentParser(description='Calculate and plot the'
        ' secondary structure exhibited by each residue in the output of'
        ' gmx do_dssp')
    parser.add_argument('filename',
                        help='.ods/xlsx file name')
    parser.add_argument('sheetname',
                        help='sheetname in the xls/ods file')
    parser.add_argument('-o', '--output',
                        help='name and extension for output image')
    return parser.parse_args()

def populateListFromSheet(filename, sheetname):
    """abc"""
    wb = xlrd.open_workbook(filename)
    sheet = wb.sheet_by_name(sheetname)
   
    newlist = []
    for i in range(sheet.nrows):
         #print(sheet.row_values(i))
         newlist.append(sheet.row_values(i))
    return newlist
   
def main():
     """abc"""
     args = get_arguments()
     list_sheet = populateListFromSheet(args.filename, args.sheetname)

     count=0
     if args.output:
               with open(args.output,"w") as file_id:
                    for ele in list_sheet:
                           count = count+1
                           dist_string="distance dist%d, i. %d and n. %s, i. %d and n. %s\n" %(count,int(ele[0]),ele[1].encode('ascii','ignore'),int(ele[2]),ele[3].encode('ascii','ignore'))
                           file_id.write(dist_string)
                           file_id.write("hide labels, dist%d\n" %(count))
     else:
         for ele in list_sheet:
               count = count+1
               dist_string="distance dist%d, i. %d and n. %s, i. %d and n. %s" %(count,int(ele[0]),ele[1].encode('ascii','ignore'),int(ele[2]),ele[3].encode('ascii','ignore'))
               print(dist_string)
               print("hide labels, dist%d\n" %(count))          

if __name__ == '__main__':
    main()     
