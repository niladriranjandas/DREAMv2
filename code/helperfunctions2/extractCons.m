%% This module is called to extract the rows of constraints
%   equality, upper or lower bound corr to atoms in the patch.
%
%  Input- cons:        constraints matrix
%         patch_index: the index of atoms belonging to the patch
%
%  Output- extracted_cons: dim same as cons in I/P
%                          rows of cons corr to atoms in patch
%
%  Note that in the indexing of the atoms in extracted_cons is different
%  from cons. For e.g. if patch_index = {2,4,...} and cons has a 
%  row [2,4,2.67], it is re-indexed to [1,2,2.67] in extracted_cons.


function extracted_cons = extractCons(cons,patch_index)
      
      if size(patch_index,1) > size(patch_index,2)
          patch_index = patch_index';
      end
      tmp = mexExtractCons(cons(:,1:3),patch_index);
      extracted_cons = [tmp, cons(1:size(tmp,1),4:end)];
%     loc_1 = ismember(cons(:,1),patch_index);
%     loc_2 = ismember(cons(:,2),patch_index);
%     cons_row = bsxfun(@and,loc_1,loc_2);
% 
%     cons_index = find(cons_row);
%     extracted_cons = zeros(size(cons_index,1),size(cons,2));
%     for j=1:size(cons_index,1)
%        extracted_cons(j,:) = [ find(patch_index==cons(cons_index(j),1)),...
%                                find(patch_index==cons(cons_index(j),2)),...
%                                cons(cons_index(j),3:end)];
%     end

end