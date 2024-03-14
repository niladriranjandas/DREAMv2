function [datafilesave, resi_index_i] = localizePartsWithAugBounds(resi_start, resi_end, datafile, pdbname)
%%
%
%%

if exist(datafile,'file')
	  load(datafile);     
else
	   error('No file found %s',datafile)
end	  

      atom_start_ = find(Comp.residue==resi_start);
      atom_end_   = find(Comp.residue==resi_end);
      atom_start  = atom_start_(1);
      atom_end    = atom_end_(end);

      atom_indices = atom_start:atom_end; %resi 104-112
      localize_i_eq_cons_all   = extractCons(eq_cons_all,atom_indices);
      %eq_cons_all_resi        = extractCons(eq_cons_all,test.cI_modified_org{i});
      localize_i_up_bounds     = extractCons(up_bounds,atom_indices);
      %up_bounds_resi          = extractCons(up_bounds,test.cI_modified_org{i});
      localize_i_wh_vdw_bounds = extractCons(wh_vdw_bounds,atom_indices);     
      %for H atoms needs to be mapped
      localize_i_vdw_bounds    = extractCons(vdw_bounds,atom_indices);
      %vdw_bounds_resi         = extractCons(vdw_bounds,test.cI_modified_org{i});
      localize_i_sdp_lo_bounds = extractCons(sdp_lo_bounds,atom_indices);
      %sdp_lo_bounds_resi      = extractCons(sdp_lo_bounds,test.cI_modified_org{i});

      localize_i_lo_bounds     = extractCons(lo_bounds,atom_indices);
      %lo_bounds_resi          = extractCons(lo_bounds,test.cI_modified_org{i});
      localize_i_wh_lo_bounds  = extractCons(wh_lo_bounds,atom_indices);
      %for with H atoms needs to be mapped

	  localize_i_X_noref = doSPROSatomwise_v2(atom_indices,...
      	                                  	  rand_X,Comp.Cq,...
        		                          localize_i_eq_cons_all,localize_i_up_bounds,localize_i_wh_vdw_bounds,localize_i_vdw_bounds,localize_i_sdp_lo_bounds,localize_i_lo_bounds,localize_i_wh_lo_bounds,... 
                                      		  hydrogen_omission);

	  [localize_i_X_refine,localize_i_info] = postProcessingMe(localize_i_X_noref, localize_i_eq_cons_all, localize_i_lo_bounds, full(localize_i_up_bounds), [1e1,1,2,-0.001],[10,10,10,10,10]);
	  [localize_i_include_index,localize_i_chk_coord_ref] = getRidOfPseudo(localize_i_X_refine',atom_indices,Comp.info);
      
      pdb_file = strcat(OUTPUT.local_folder_pdb, filesep, pdbname);%'/data2/nmr/our_algo/output_pdb/1n6u_resi_105_112.pdb';
      resi_index_i={writeToPDB(pdb_file,localize_i_include_index,localize_i_chk_coord_ref',Comp)}

      datafilesave = sprintf('%s%slocalize_resi_%d_%d.mat',INPUTS.protein_path,filesep,resi_start,resi_end);
%      datafilesave = strcat(INPUTS.protein_path, filesep, "localize_resi_", resi_start, "_", resi_end, ".mat");
      save(datafilesave);

end
