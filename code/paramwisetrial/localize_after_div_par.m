%time.div_conq_prep = toc(div_conq_prep);
   %% loop through the cells in test.modified_cI, localize each
div_conq_loc_all = tic;
[localize] = doLocalizeFragment(test.cI_modified,hydrogen_omission,rand_X,Comp.Cq, Comp.info, test.dup_flag, test.cI_modified_org, eq_cons_all, up_bounds, wh_vdw_bounds, vdw_bounds, sdp_lo_bounds, lo_bounds, wh_lo_bounds);  
time.div_conq_loc_all = toc(div_conq_loc_all);
  %% ==========write into pdb file and pymol script==========================

 %--------------pdb file---------------------------------------------------
  resi_index = cell(1,length(localize));
display('__________writing pdb files for each group _______________')
   %local_folder_pdb = '/home/niladri/Documents/Disco_etc_all_in_1/SPROS_try/localization_try/pdb_resi';   
   local_folder_pdb = OUTPUT.local_folder_pdb;
          if ~exist(local_folder_pdb,'dir')
             mkdir(local_folder_pdb);
          end
  for i=1:length(localize)
   if localize(i).method > 0
    fprintf('\n____Writing PDB for grp-%d____\n',i);
     file     = sprintf('gr_%s_%u.pdb',protein_name,i);
     pdb_file = strcat(local_folder_pdb,filesep,file);
     %resi_index(i)={writeToPDB(pdb_file,localize(i).atoms,localize(i).X_refine,Comp)};
     resi_index(i)={writeToPDB(pdb_file,localize(i).include_index,localize(i).chk_coord_ref',Comp)};
     LOG.curr_txt = sprintf('\n %s written.\n',pdb_file);
     LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
   else
      fprintf('\n____skiping PDB file for grp-%d____\n',i);
   end
  end
  
  %% report violations before and after refinement, and write logs
  
  fprintf('\n Groups created: %d / solved: %d / Duplicate: %d / Not solved: %d\n',length(localize), ...
                                                        length(find([localize.method]==1)), ...
                                                        length(find([localize.method]==-1)),...
                                                        length(find([localize.method]==0)));
  LOG.curr_txt = sprintf('\n Groups created: %d / solved: %d / Duplicate: %d / Not solved: %d\n',length(localize), ...
                                                        length(find([localize.method]==1)), ...
                                                        length(find([localize.method]==-1)),...
                                                        length(find([localize.method]==0)));
  LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
                                                      
  
  fprintf('\n=================================== Violations (Stage: Divide and conquer) =========================================');
  LOG.curr_txt = sprintf('\n=================================== Violations (Stage: Divide and conquer) =========================================');
  LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
  for i=1:length(localize)
      if localize(i).method > 0   % 1 for spros 2 for cvx; -1 for duplicate 0 for could not solve
        fprintf('\n************************************Group No. %d ***************************************************************',i);
        LOG.curr_txt = sprintf('\n************************************Group No. %d , Method: %d*************************************************',i,localize(i).method);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
          
        fprintf('\n Before refinement');
        LOG.curr_txt = sprintf('\n Before refinement');
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
        
        [violation.mean_err, violation.min_err, violation.max_err, violation.err_percent] = calcViolations(localize(i).X_noref,localize(i).eq_cons_all,3);
        fprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
        
        [violation.mean_err, ~, violation.max_err, violation.err_percent] = calcViolations(localize(i).X_noref,localize(i).up_bounds,1);
        fprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);

        [violation.mean_err, violation.min_err, ~, violation.err_percent] = calcViolations(localize(i).X_noref,localize(i).sdp_lo_bounds,2);
        fprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);

        fprintf('\n\n After refinement');
        LOG.curr_txt = sprintf('\n\n After refinement');
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
        
        [violation.mean_err, violation.min_err, violation.max_err, violation.err_percent] = calcViolations(localize(i).X_refine,localize(i).eq_cons_all,3);
        fprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
        
        [violation.mean_err, ~, violation.max_err, violation.err_percent] = calcViolations(localize(i).X_refine,localize(i).up_bounds,1);
        fprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);

        [violation.mean_err, violation.min_err, ~, violation.err_percent] = calcViolations(localize(i).X_refine,localize(i).sdp_lo_bounds,2);
        fprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
        
        fprintf('\n Time taken: %d', localize(i).time);
        LOG.curr_txt = sprintf('\n Time taken: %d', localize(i).time);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
        
      end
  end
  
  LOG.curr_txt = sprintf('Module: Divide-n-conquer: time for div-n-expand: %d \t time for conquer: %d',time.div_conq_prep,time.div_conq_loc_all);
  LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
  
  LOG.prev_txt = mexDolog(LOG.prev_txt,'\n_____________END of localization of groups __________\n',1,LOG.div_conquer);
