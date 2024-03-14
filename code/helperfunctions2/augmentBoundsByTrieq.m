function tri_bounds = augmentBoundsByTrieq(bounds,atoms)
%% function to augment the upper bounds with additional bounds using the triangle inequality
%   Input: bounds: 3*k array (i,j,d) where d is upper dist bounds between
%                  atom i and j
%            atom: total no. of atoms
%  Output: tri_bounds: upper bounds derived from the triangle inequality
%                      (i,j,d)
%%
   if nargin ~= 2
       error('Module: augmentBoundsByTrieq: incorrect no. of inputs.');
   end
   
   bounds = reshape_i_gt_j(bounds); % bounds: (i,j,d) such that i<j
%% build the adjacency matrix and square it
   adj_mat    = sparse(bounds(:,1),bounds(:,2),bounds(:,3),atoms,atoms);
   %adj_mat(logical(eye(atoms))) = 0;     % this is req coz due to pseudo atoms sometimes (i,i) 
   %                                      % may be present in bounds
   adj_mat    = adj_mat + adj_mat';  
   adj_mat_01 = adj_mat>0;
   
   adj_mat_2 = adj_mat_01^2;
   
%% loop through adj_mat_2 and get the triangles
     [vert_i,vert_j,~] = find(triu(adj_mat_2,1));
     count=0;
     
     for i=1:length(vert_i)
        if ~adj_mat(vert_i(i),vert_j(i))  %add new bounds if its not already there
          neigh_vert_i = find(adj_mat(vert_i(i),:));
          neigh_vert_j = find(adj_mat(vert_j(i),:));
          comm_vert = intersect( neigh_vert_i, neigh_vert_j);
          dist = adj_mat(vert_i(i),comm_vert(1)) + adj_mat(comm_vert(1),vert_j(i));
           for j=2:length(comm_vert)
              dist = min(dist,  adj_mat(vert_i(i),comm_vert(j)) + adj_mat(comm_vert(j),vert_j(i)));              
           end
          count = count +1;
          tri_bounds(count,:) = [vert_i(i),vert_j(i),dist];
        end
     end
end 

function bounds = reshape_i_gt_j(raw_bounds)
  count = 0;
  for i=1:size(raw_bounds,1)
     if raw_bounds(i,1) < raw_bounds(i,2)
         count  = count+1;
         bounds(count,:) = raw_bounds(i,:);
     elseif raw_bounds(i,1) > raw_bounds(i,2)
         count = count+1;
         bounds(count,:) = [raw_bounds(i,2), raw_bounds(i,1), raw_bounds(i,3)];
     end
  end
end