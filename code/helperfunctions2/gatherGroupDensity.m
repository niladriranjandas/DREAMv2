function [allgrp_density,global_density,grp_density] = gatherGroupDensity(grps,adj_mat)
%% function to calculate the groupwise graph density (defined as #(edges in each grp_i)/#(complete graph of size grp_i)
%  Input:   grps: 1 x k cells each cell 1 x k_i double containing the
%                list of vertices 
%        adj_mat: adjacency matrix corr to the vertices in the grps
% Output:   global_density: density of the full adj_mat
%              grp_density: density of the individual grps
%%
   if nargin ~=2
       error('Module gatherGroupDensity: Not enough inputs.');
   end
   
   if ~isuppertriangular(adj_mat)
       if ~issymmetric(adj_mat)
           error('Module gatherGroupDensity: adj_mat neither upper triangular nor symmetric');
       else
          adj_mat = triu(adj_mat ,1);
       end
   end
   
   adj_mat = adj_mat > 0; % change to binary adj_mat
        
%% anonymous func for calc graph density
    density = @(vertex_list,adj_mat) (nnz(adj_mat(vertex_list,vertex_list))/(length(vertex_list) * (length(vertex_list)-1)));
  
%%
   allgrp_density = density(unique([grps{1:end}]),adj_mat);
   global_density = density(1:length(adj_mat),adj_mat);
   grp_density    = zeros(1,length(grps));
   
   for i=1:length(grps)      
       grp_density(i) = density(grps{i},adj_mat);
   end

end