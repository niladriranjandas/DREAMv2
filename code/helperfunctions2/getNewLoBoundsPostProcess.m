function [ret_struct, lo_bounds] = getNewLoBoundsPostProcess(X, atom_map, eq_cons_all, up_bounds, lo_bounds, Comp, pdb_file_name)
%%
%
%%
   new_lo = getNewLoBounds(X, atom_map, Comp);
   lo_bounds_org = lo_bounds;
   lo_bounds     = [lo_bounds;new_lo];

   reg_struct.eq_cons_all   = extractCons(eq_cons_all,atom_map);
   reg_struct.up_bounds     = extractCons(up_bounds,atom_map);
   %reg_struct.reg_sdp_lo_bounds = extractCons(sdp_lo_bounds,atom_map);
   reg_struct.lo_bounds     = extractCons(lo_bounds,atom_map);

       fprintf('\n__________refine module___\n');
        [ret_struct.reg_X_refine,ret_struct.reg_info] = postProcessingMe(X,...
                                                   reg_struct.eq_cons_all, reg_struct.lo_bounds, full(reg_struct.up_bounds),...
                                                   [1e4,100,10,-0.001],[10,10,10,10,10]);
        %---------------extract the original PDB coordinates---------------------
        [ret_struct.reg_include_index,ret_struct.reg_chk_coord_ref] = getRidOfPseudo(ret_struct.reg_X_refine',atom_map,Comp.info);
                ret_struct.reg_resi_index={writeToPDB(pdb_file_name,ret_struct.reg_include_index,ret_struct.reg_chk_coord_ref',Comp)};

        %plot_filename_rama = sprintf('%s_ramachandran_run_1.fig',protein_name);
        %plot_file_rama    = strcat(local_folder_plots_postreg,filesep,plot_filename_rama);

        %[rama.out,fig_handler] = ramachandranMe(pdb_file_name,'Regions',true,'Glycine',true);

        %resi in various area in the ramachandran map
        %[rama.region,rama.part_tot_resi] = getRamachandranReigionDistri(rama.out.Angles);

 %% get percent disallowed region for ramachandran map

     %rama.disallowed = length(rama.region(end).resi) / length([rama.region.resi]);

end
