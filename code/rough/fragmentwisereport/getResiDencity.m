function residensity = getResiDencity(resi_list,wh_Comp,wh_up_bounds,wh_sdp_lo_bounds,wh_eq_bounds,min_res,max_res)
%%

%%
   
   resi_list_new = reindexResilist(resi_list,min_res, max_res);
   
    %[resi_adj_mat,resi_adj_mat_eq,atom_adj_mat,atom_adj_mat_eq] = buildResiAdjMat(num_atoms,eq_cons,up_bounds,sdp_lo_bounds,Comp)
    n_atoms = length(wh_Comp.residue);
    %n_atoms = length(resi_list);
   [resi_adj_mat_org,resi_adj_mat_eq,adj_mat,adj_mat_weq] = ...
      buildResiAdjMat(n_atoms,wh_eq_bounds,wh_up_bounds,wh_sdp_lo_bounds,wh_Comp);
    resi_adj_mat = resi_adj_mat_org(resi_list_new,resi_list_new);
    
    
    length_resiadj = length(resi_adj_mat);
    residensity=nnz(resi_adj_mat)/(length_resiadj * (length_resiadj-1));
end

function resilist_new = reindexResilist(resilist, minres, maxres)
   index_tmp = minres:maxres;
   resilist_new = nan(1,length(resilist));
   for i=1:length(resilist)
      resilist_new(i) = find(index_tmp == resilist(i)); 
   end
end
