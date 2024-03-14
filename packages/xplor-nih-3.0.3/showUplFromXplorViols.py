import sys
import re

# example run: python showUplFromXplorViols.py egg_viols.txt 18 34

STRONG = 2.7
MEDIUM = 3.5
LO = 6

def displayUplEntries(violfile, up1, up2, tol=0):
    '''

    '''    

    count = 0
    with open(violfile) as up:
         lines = up.readlines()
         num_lines = len(lines)
         for line in lines:
                #print(line)
                if "number of restraints" in line:   # this tends to be the last line in viol file generated
                    break

                count = count +1
                trim_spc = line.strip()
                if trim_spc:
                   if trim_spc[0].isdigit() :                
                        parts = line.split(')')

                        resi_atom_i = parts[0].split('(')[1]
                        resi_i = resi_atom_i.split()
                        resi_i_residue_no = resi_i[0]

                        resi_atom_j = parts[1].split('(')[1]
                        resi_j = resi_atom_j.split()
                        resi_j_residue_no = resi_j[0]

                        rest_parts = parts[2].strip().split()
                        calc_dist = rest_parts[0]
                        givn_dist = rest_parts[1]
                        givn_dist_parts = givn_dist.split('..')  
                        exp_dist = float(givn_dist_parts[1])                      
                        viol_diff = float(rest_parts[2])

                        if viol_diff > tol:
                            if (resi_i_residue_no >= up1 and resi_i_residue_no <= up2) or (resi_j_residue_no >= up1 and resi_j_residue_no <= up2):
                                if (count == num_lines):
                                     if exp_dist>0 and exp_dist <= STRONG:
                                          print(' %s\t%s\t%s\t%s\t%f\tSTRONG'%(resi_i,resi_j,calc_dist,givn_dist,viol_diff))
                                     elif exp_dist>STRONG and exp_dist <= MEDIUM:
                                          print(' %s\t%s\t%s\t%s\t%f\tMEDIUM'%(resi_i,resi_j,calc_dist,givn_dist,viol_diff))
                                     elif exp_dist>MEDIUM and exp_dist <= LO:   
                                          print(' %s\t%s\t%s\t%s\t%f\tLO'%(resi_i,resi_j,calc_dist,givn_dist,viol_diff))                                     
                                     
                                     #print(' %s\t%s\t%s\t%s\t%f'%(resi_i,resi_j,calc_dist,givn_dist,viol_diff))
                                     #print(' %s'%(line))
                                else:
                                     nextline = lines[count].strip()
                                     if nextline:
                                       if nextline[0] == '(':
                                          if exp_dist>0 and exp_dist <= STRONG:
                                                print('*%s\t%s\t%s\t%s\t%f\tSTRONG'%(resi_i,resi_j,calc_dist,givn_dist,viol_diff))
                                          elif exp_dist>STRONG and exp_dist <= MEDIUM:
                                                print('*%s\t%s\t%s\t%s\t%f\tMEDIUM'%(resi_i,resi_j,calc_dist,givn_dist,viol_diff))
                                          elif exp_dist>MEDIUM and exp_dist <= LO:   
                                                print('*%s\t%s\t%s\t%s\t%f\tLO'%(resi_i,resi_j,calc_dist,givn_dist,viol_diff))                                     

                                          #print('*%s\t%s\t%s\t%s\t%f'%(resi_i,resi_j,calc_dist,givn_dist,viol_diff))
                                       else:
                                          if exp_dist>0 and exp_dist <= STRONG:
                                                print(' %s\t%s\t%s\t%s\t%f\tSTRONG'%(resi_i,resi_j,calc_dist,givn_dist,viol_diff))
                                          elif exp_dist>STRONG and exp_dist <= MEDIUM:
                                                print(' %s\t%s\t%s\t%s\t%f\tMEDIUM'%(resi_i,resi_j,calc_dist,givn_dist,viol_diff))
                                          elif exp_dist>MEDIUM and exp_dist <= LO:   
                                                print(' %s\t%s\t%s\t%s\t%f\tLO'%(resi_i,resi_j,calc_dist,givn_dist,viol_diff))                                     

    up.close()


def main():
    violfile = sys.argv[1]
    up1 = sys.argv[2]
    up2 = sys.argv[3]

    displayUplEntries(violfile, up1, up2)

if __name__ == "__main__":
    main()
