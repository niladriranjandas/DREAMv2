function extracted_cons = extractCons_old(cons,patch_index)
        
    loc_1 = ismember(cons(:,1),patch_index);
    loc_2 = ismember(cons(:,2),patch_index);
    cons_row = bsxfun(@and,loc_1,loc_2);

    cons_index = find(cons_row);
    extracted_cons = zeros(size(cons_index,1),size(cons,2));
    for j=1:size(cons_index,1)
       extracted_cons(j,:) = [ find(patch_index==cons(cons_index(j),1)),...
                               find(patch_index==cons(cons_index(j),2)),...
                               cons(cons_index(j),3:end)];
    end

end