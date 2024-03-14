function modified_cI = changeToSeq(cI,jump)
%% function to convert cI to seq like (i,i+1,....)(j,j+1,...) for each cell in cI
%   Input:          cI: 1*n cell each cell containing atoms no.(s) in
%                       increasing order
%  Output: modified_cI: 1* n cell each cell having sequences like
%                       (i,i+1,....,i+k)(j,j+1,...). Jump between i+k and j
%                       is taken as 50
%% check input

   if nargin ~= 2
       error('Module: changeToSeq: incorrect number of inputs.');
   end
   
   if jump <=0
       error('Module: changeToSeq: jump argument must be positive integer.');
   end
   
   num_cells = length(cI);
%%  
  %jump = 25;%50;
  
    modified_cI = cell(1,num_cells);
    for i=1:num_cells
       mod_cell_row = [];
       succ_diff = cI{i}(2:end) - cI{i}(1:end-1);
       
       ind = find(succ_diff>jump);
       
       prev=1;
       
       if length(ind)==0
         mod_cell_row = cI{i}(1):cI{i}(end);
       else
         for j = 1:length(ind)
             mod_cell_row = [mod_cell_row, cI{i}(prev):cI{i}(ind(j))];
             prev = ind(j)+1;
         end       
         mod_cell_row = [ mod_cell_row, cI{i}(prev):cI{i}(end)];
       end
      modified_cI(i) = {mod_cell_row} ;
    end
end    