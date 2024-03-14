function [back_bone_tmp_index, back_bone_atom_map] = onlyTheBackBone(X,atom_map,Comp)
%% Intended to be called after reg module to get the indices of only the back bone.
%   Input:      X: 3 x k coordinates.
%        atom_map: global atom numbers corr to the atoms in X.
%  Output: back_bone_tmp_index: 1 x m array indicating indices in X which
%                               corr to backbone atoms.
%           back_bone_atom_map: 1 x m global atom number corr to 
%                               back_bone_tmp_index
%%
   if nargin ~=3
       error('module: onlyTheBackBone: incorrect number of arguments.');
   end
   
   [dim,n_atom] = size(X);
   
   if dim ~=3
      error('module: onlyTheBackBone: X must be 3 x k.'); 
   end
   
   if n_atom ~= length(atom_map)
      error('module: onlyTheBackBone: length of atom_map and X doesn''t match.'); 
   end   
   
%%
   back_bone_tmp_index = []; back_bone_atom_map = []; 
  
   count = 0;
   X_residue = Comp.residue(atom_map);
   X_resi    = unique(Comp.residue(atom_map));
   
   for i=1:length(X_resi)
      bck_bn_ind_resi_i   = getBackBoneResiInd(X_resi(i),Comp.atom_names,Comp.residue); 
      X_tmp_ind           = find(X_residue == X_resi(i));
      X_bck_bn_ind_resi_i = atom_map(X_tmp_ind);
      
      mem_ind = ismember(bck_bn_ind_resi_i, X_bck_bn_ind_resi_i);
      n_atoms_bk_bn = length(bck_bn_ind_resi_i);
      
      if sum(mem_ind) == n_atoms_bk_bn
          tmp_ind   = find(ismember(X_bck_bn_ind_resi_i,bck_bn_ind_resi_i));
          back_bone_tmp_index(count+1:count+n_atoms_bk_bn) = X_tmp_ind(tmp_ind);
          back_bone_atom_map(count+1:count+n_atoms_bk_bn)  = bck_bn_ind_resi_i;
          count = count + n_atoms_bk_bn;          
      end
   end
end

function bck_bn_atm_resi_i = getBackBoneResiInd(resi,Comp_atom_names,Comp_residue)
       atm_srch          = {'N','CA','C','O'};
       bck_bn_atm_resi_i = zeros(1,length(atm_srch));
       
       atms  = find(Comp_residue == resi);
       count = 1;
       
       for i=1:length(atms)
          if strcmp(Comp_atom_names(atms(i)),atm_srch(count))
             bck_bn_atm_resi_i(count) = atms(i);
             count = count+1;
          end
       end
   
end