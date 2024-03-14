function resi_index = writeToPDB(file_name,atom_num, X,Comp)
%% This function calls writes only the residues for which the atoms are available.
%  Calls pdb_writer
%   Input:  file_name: output file name
%            atom_num: index of the atoms
%                   X: coordinates of the atoms (3xn)
%                Comp:  Comp structre
% Output:   resi_list: residues which has been written to the file.
%% check input


%%
   %[resi_list,pos_index] = getResiList(atom_num,Comp.residue); %this doesn't work
   
   resi_index = unique(Comp.residue(atom_num)); pos_index = atom_num;
   % prepare the Comp.residue, Comp.info, Comp.atom_names, Comp.residue_type for pdb_writer
   
   Comptmp.info         = Comp.info(:,pos_index);
   Comptmp.residue      = Comp.residue(pos_index);
   Comptmp.atom_names   = Comp.atom_names(pos_index);
   Comptmp.residue_type = Comp.residue_type(pos_index);
   
   % write into file
   tmp_index = zeros(1,length(pos_index));
   for i=1:length(pos_index)
      tmp_index(i) = find(atom_num == pos_index(i)); 
   end
   
   pdb_writer(file_name,X(:,tmp_index),Comptmp);

end

function [resi_list,pos_index] = getResiList(atom_num,Comp_residue)
%%

   resi_grp = Comp_residue(atom_num);
   
   resi_uniq = unique(resi_grp);
   pos_index = [];
   count = 0;
   
   for i=1:length(resi_uniq)
      org_pos_ind = find(Comp_residue == resi_uniq(i));
      X_pos_ind   = find(resi_grp == resi_uniq(i));
      
      if length(org_pos_ind) == length(X_pos_ind)
         count = count+1;
         resi_list(count) = resi_uniq(i);
         pos_index = [pos_index,atom_num(X_pos_ind)];
      end
   end
end
   