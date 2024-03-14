function indices = getCompFromWhComp(Comp_atm_names, wh_Comp_atm_names, Comp_residues, wh_Comp_residues)
%% returns the position of Comp atom names (Comp.atom_names from Wh_Comp.atom_names 
%   Input: Comp.atom_names, wh_Comp.atom_names
%  Output: indices (1 x size of Comp.atom_names) if any atom name is not
%          found Nan is returned
%%
  
    if length(Comp_atm_names) >= length(wh_Comp_atm_names)
           error('Module: getCompFromWhComp: size of Comp.atom_names > wh_Comp.atom_names');
    end
    
    %check if residues are same or not
    residues = unique(Comp_residues);
    if ~ all( (residues - unique(wh_Comp_residues))==0)
        error('Module: getCompFromWhComp: resdues mismatch between Comp.residues and wh_Comp.residues')
    end
    
%%
    sum = 0;
    indices = nan(1,length(Comp_atm_names));
    for i = 1:length(residues)
        indx_Comp    = find(Comp_residues == residues(i));
        indx_wh_Comp = find(wh_Comp_residues == residues(i));
        
        indicator = ismember( wh_Comp_atm_names(indx_wh_Comp), Comp_atm_names(indx_Comp));
        pos = find(indicator);
        
        indices(sum+1:sum+length(pos)) = indx_wh_Comp(pos);
        sum = sum + length(pos);
    end
    
    if any(isnan(indices))
           error('Module: getCompFromWhComp: some indices not found');
    end
end