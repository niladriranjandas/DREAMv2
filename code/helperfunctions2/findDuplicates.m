function [flag,dup_info] = findDuplicates(patch_overlap_mat)
%% function to remove the duplicate entries from the test.cI
%  Input:  patch_overlap_mat: num_patches x num_patches matrix containing
%                             the number of overlapping points. (i,i) 
%                             contains the size of the patch-i.
% Output:               flag: 1 if no duplicates, 0 if duplicates. for 
%                             patch (i,j,k,...) debing duplicates only i is
%                             retained. (i<j<k<..)
%                   dup_info: duplicate info for each group in a cell
%%
   if ~issymmetric(patch_overlap_mat)
       if ~all(all(tril(patch_overlap_mat)==0)==0)
           error('Module: findDuplicates: incoreect patch_overlap_mat.');
       end
   else
       patch_overlap_mat = triu(patch_overlap_mat);
   end
   
   if any(diag(patch_overlap_mat)==0)
       error('Module: findDuplicates: incoreect patch_overlap_mat.');
   end
   
%%
   num_patch = length(patch_overlap_mat);
   dup_info = cell(1,num_patch);
   flag = ones(1,num_patch);
   
   for i=1:num_patch-1
      [~,col,~] = find(patch_overlap_mat(i,:)==patch_overlap_mat(i,i));
      
      tmp = [];
      if length(col)>1 % col will ofcourse contain i by default
        for j=col(2:end)
          if patch_overlap_mat(j,j) == patch_overlap_mat(i,i)
             flag(j) = 0;
             tmp = [tmp,j];
          end
        end
        dup_info(i) = {tmp};
      end
   end