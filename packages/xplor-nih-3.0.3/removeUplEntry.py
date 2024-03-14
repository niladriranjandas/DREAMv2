import sys
import re

def removeUplEntries(uplfile, pdbfile, oppuplfile):
    '''

    '''

    pdb = open(pdbfile, 'r')
    pdb_lines = pdb.read()

    with open(uplfile) as up, open(oppuplfile,'w') as oppupl:
         lines = up.readlines()
         for line in lines:
                #print(line)
                if 'assign' in line :
                        #print(line)
                        #print("---------------------------------")
                        parts = line.split(')')

                        resi_atom_i = parts[0].split('(')[1]
                        #print(resi_atom_i)
                        resi_i = resi_atom_i.split()
                        resi_i_residue_no = resi_i[1]
                        resi_i_atom_name = resi_i[4]

                        resi_atom_j = parts[1].split('(')[1]
                        #print(resi_atom_j)
                        resi_j = resi_atom_j.split()
                        resi_j_residue_no = resi_j[1]
                        resi_j_atom_name = resi_j[4]
                    
                        flag1 = ['1']
                        flag2 = ['1']

                        if not resi_i_atom_name.endswith('*'):
                               srch_str1 = "\s%s\s.* %s "%(resi_i_atom_name,resi_i_residue_no)
                               flag1 = re.findall(srch_str1, pdb_lines)
                        if not resi_j_atom_name.endswith('*'):
                               srch_str2 = "\s%s\s.* %s "%(resi_j_atom_name,resi_j_residue_no) 
                               flag2 = re.findall(srch_str2, pdb_lines)

                        if resi_i_atom_name.endswith('*'):
                            resi_i_atom_name = resi_i_atom_name[:-1]
                            srch_str1 = "\s%s[0-9]*\s.* %s "%(resi_i_atom_name,resi_i_residue_no)
                            flag1 = re.findall(srch_str1, pdb_lines)
                        if resi_j_atom_name.endswith('*'):
                            resi_j_atom_name = resi_j_atom_name[:-1]
                            srch_str2 = "\s%s[0-9]*\s.* %s "%(resi_j_atom_name,resi_j_residue_no)
                            flag2 = re.findall(srch_str2, pdb_lines)                               
                       
                        #print('%s - %s'%(srch_str1,len(flag1)))
                        #print('%s - %s'%(srch_str2,len(flag2)))
                        if len(flag1) == 0 or len(flag2) == 0 : 
                              continue
                        else :
                              oppupl.write('%s'%(line))

    pdb.close()


def main():
    uplfile = sys.argv[1]
    pdbfile = sys.argv[2]
    oppuplfile = sys.argv[3]
    
    removeUplEntries(uplfile, pdbfile, oppuplfile)

if __name__ == "__main__":
    main()