function atom_seq_cI = changeToAtomSeq(cI_resi,Comp_residue)
%% function to get the atom seq from the resi numbers
%   Input:      cI_resi: test.cI_resi = breakByDenseGraph(test.resi_adj_mat>0)
%          Comp_residue: Comp.residue
%  Output:  atom_seq_cI: the resi seq in cell of test.cI_resi changed to
%                        atoms seq
%%
   if nargin ~=2
       error('Module: changeToAtomSeq: incorrect number of inputs.');
   end
   
%%
   num_grp = length(cI_resi);   
   atom_seq_cI = cell(1,num_grp);
   
   for i=1:num_grp
     tmp = cI_resi{i};
     tmp_grp = [];
     for j=1:length(tmp)
        tmp_grp = [tmp_grp,find(Comp_residue==tmp(j))']; 
     end
    atom_seq_cI(i) = {tmp_grp};
   end
   
end