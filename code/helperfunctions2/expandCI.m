function cI_expanded = expandCI(adj_mat,cI_resi,k)
%% If a residue in {all residues}\cI_resi(i) has more than k common residues
%  then it is included in cI_resi(i)
%   Input:     adj_mat: resi x resi adjacency matrix each cell containing
%                       the number of upper and lower bounds between them.
%              cI_resi: test.cI_resi
%                    k: 
%  Output: cI_expanded: expanded cI_expanded.
%%
   if nargin~=3
       error('module: expandCI: incorrect number of arguments.');
   end

   if size(adj_mat,1) ~= size(adj_mat,2)
       error('module: expandCI: adjacency matrix incorrect.');
   end
   
   if ~isuppertriangular(adj_mat)
       if issymmetric(adj_mat)
           adj_mat = (triu(adj_mat,1))>1;
       else
           error('module: expandCI: adjacency matrix neither symmetric nor upper triangular.');
       end
   end
%%
   %k = 4;
   
   grp = length(cI_resi);   
   cI_expanded = cell(1,grp);
   
   for i=1:grp
      neighs      = findNeighs(adj_mat,cI_resi{i},k);
      cI_expanded(i) = {sort(union(cI_resi{i},neighs))};
   end
   
end

function neighs = findNeighs(adj_mat,tmp_resi,k)
%% function to get the neighbours of tmp_resi having greater than k 
   if isuppertriangular(adj_mat)
       adj_mat = adj_mat+adj_mat';
   end
   
   members = setdiff(1:length(adj_mat),tmp_resi);
   
   neigh_mat = adj_mat(tmp_resi,members);
   
   tmp_index = find(sum(neigh_mat,1)>k);
   neighs = members(tmp_index);
   
end
