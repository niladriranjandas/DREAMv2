function [resi_adj_mat,resi_adj_mat_eq,atom_adj_mat,atom_adj_mat_eq] = buildResiAdjMat(num_atoms,eq_cons,up_bounds,sdp_lo_bounds,Comp)
%% function to create the resi to resi constraint adjacent matrix
%  Input:      num_atoms: No. of atoms
%                eq_cons: list of the equality constraints
%              up_bounds: list of the upper bounds
%          sdp_lo_bounds: list of the sdp lo bounds (ang and noessy lo bounds
%                   Comp: the Comp structure
%
%  Output:  resi_adj_mat: #resi \times #resi matrix such that (i,j) is 1 if
%                         bounds (up or lo) exist between them
%%
       if nargin ~= 5
           error('Module:buildResiAdj: Not enough inputs.');
       end
       
%% build the atom to atom adjacency matrix first

    adj_mat = bsxfun(@or, sparse(up_bounds(:,1),up_bounds(:,2),1,num_atoms,num_atoms), ...
                       sparse(sdp_lo_bounds(:,1),sdp_lo_bounds(:,2),1,num_atoms,num_atoms));
          adj_mat  = (adj_mat + adj_mat')>0;
      %because of pseudo atoms adj_mat(i,i) may not be zero.hence make them
      %zero
          adj_mat(logical(eye(length(adj_mat)))) = zeros(1,length(adj_mat));
          atom_adj_mat = adj_mat;
                   
    adj_mat_weq = bsxfun(@or, adj_mat, ...
                            sparse(eq_cons(:,1),eq_cons(:,2),1,num_atoms,num_atoms));
          adj_mat_weq = (adj_mat_weq + adj_mat_weq')>0;
          atom_adj_mat_eq = adj_mat_weq;

%     if ~isuppertriangular(adj_mat)
%         adj_mat = triu(adj_mat,1);
%     end
%     
%     if ~isuppertriangular(adj_mat_weq)
%         adj_mat_weq = triu(adj_mat_weq,1);
%     end
    
%%
   num_res = length(Comp.residue_bias);
   
   resi_list = Comp.num_seq;
   
  resi_adj_mat    = zeros(num_res);
  resi_adj_mat_eq = zeros(num_res);
   
   for i=1:num_res
    % [st_atom_i,nd_atom_i] = getAtomFromResi(i,Comp.residue_bias,num_atoms);
     [st_atom_i,nd_atom_i] = getAtomFromResi(resi_list(i),Comp.residue);
       for j=i+1:num_res
           %[st_atom_j,nd_atom_j] = getAtomFromResi(j,Comp.residue_bias,num_atoms);
           [st_atom_j,nd_atom_j] = getAtomFromResi(resi_list(j),Comp.residue);
             if nnz(intersect(st_atom_i:nd_atom_i,st_atom_j:nd_atom_j))
                 error('something not right.');
             end
           resi_adj_mat(i,j)    = nnz(adj_mat(st_atom_i:nd_atom_i,st_atom_j:nd_atom_j));
           resi_adj_mat_eq(i,j) = nnz(adj_mat_weq(st_atom_i:nd_atom_i,st_atom_j:nd_atom_j));
       end
   end
end

function [st_atom,nd_atom] = getAtomFromResi(resi,Comp_residue)
    atm_list = find(Comp_residue==resi);
    st_atom = atm_list(1);
    nd_atom = atm_list(end);
end

% function [st_atom,nd_atom] = getAtomFromResi(resi,residue_bias,num_atoms)
% %% function to return the atom list for the atom
%   if resi == 1
%       st_atom =1;
%   else
%       st_atom = residue_bias(resi)+1;
%   end
%   
%   if resi == length(residue_bias)
%      nd_atom = num_atoms;
%   else
%      nd_atom = residue_bias(resi+1);
%   end
%      
%    
% end
   
   