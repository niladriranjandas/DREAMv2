function overlap_mat = getPatchOverlap(cI)
%% function to return the overlap matrix between the cells in cI
%   Input:          cI: 1*n cells containing the atom no.(s)
%  Output: overlap_mat: n * n upper triangular matrix
%                       where (i,j) denotes the no. of common points
%                       between cell i and j.
%% check inputs
  if nargin ~=1
      error('Module: getPatchOverlap: incorrect number of inputs.');
  end
  
%%
 
  num_cells = length(cI)
  overlap_mat = zeros(num_cells);
  
   % assign the diagonal elements
   overlap_mat(1:num_cells+1:end) = cellfun(@length,cI);%cellfun(@length calculates length of each cell
     for i=1:num_cells-1
         for j=i+1:num_cells
            overlap_mat(i,j) = length(intersect(cI{i},cI{j}));
         end         
     end
end     