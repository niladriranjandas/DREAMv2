%% ---------------------- Divide and Conquer ------------------------------
%   - Model as residue graph: 
%   * break the graph: 
%   * include neighbours:
%   * expand groups over multiple iterations:
%   * include small jumps:
%   - localize each group
%   - refine (BFGS)
%%
   hydrogen_omission = PARAMS.hydrogen_omission;
   protein_name      =  INPUTS.protein_name;
   
  LOG.div_conquer = strcat(LOG.div_n_conquer,filesep,'div_n_conquer_',protein_name,'.txt');
  LOG.prev_txt = '';
%%  preprocess
div_conq_prep = tic;
  preprocess;
  
  test.n_atoms = size(rand_X,2);
  test.n_eq_cons = size(eq_cons,1);
  test.n_up_cons = size(up_bounds,1);
  test.n_lo_cons = size(sdp_lo_bounds,1);
  
   if hydrogen_omission
     wh_eq_cons_all = equality_con_former(wh_rand_X,wh_Comp,2);
     eq_cons_all    = equality_con_former(rand_X,Comp,2);  % instead of the base equality cons get all the equality cons
   else
     wh_eq_cons_all = equality_con_former(rand_X,Comp,2);
     eq_cons_all    = wh_eq_cons_all;
   end  
   
 %% augment up_bounds by triangle inequality
     take_aug_bounds =  PARAMS.aug_bounds ;%8;%7.5;
     if take_aug_bounds > 0
         up_bounds_org = up_bounds;   
         augment_bounds = augmentBoundsByTrieq(up_bounds(:,1:3),test.n_atoms);
         aug_I = find(augment_bounds(:,3)<take_aug_bounds);
         augment_bounds = [ augment_bounds(aug_I,:), zeros(length(aug_I),1)];
         up_bounds = [up_bounds; augment_bounds];
     end

     % augment up_bounds again
     up_bounds_before_multi_augbounds = up_bounds;
     if isfield(PARAMS, 'aug_bounds_again_iter')
        for i=1:PARAMS.aug_bounds_again_iter
                augment_bounds = augmentBoundsByTrieq(up_bounds(:,1:3),test.n_atoms);
                aug_I = find(augment_bounds(:,3)<PARANS.aug_bounds_again_bounds);
                augment_bounds = [ augment_bounds(aug_I,:), zeros(length(aug_I),1)];
                up_bounds = [up_bounds; augment_bounds];
        end
     end
  
 % for the adjacency matrix
 [test.resi_adj_mat,test.resi_adj_mat_eq,test.adj_mat,test.adj_mat_weq] = ...
      buildResiAdjMat(test.n_atoms,eq_cons_all,up_bounds,sdp_lo_bounds,Comp);   
  
 %% call the breaking algorithm

    test.cI_resi = preBreakByDenseGraph((test.resi_adj_mat+test.resi_adj_mat')>0,3,PARAMS.eta_lo);
 %test.cI_resi = mapToResi(test.cI_resi,Comp.num_seq);  %included after bug_report_4_03_2016
 
  %-------------grow each group to include neighbours----------------------
  k = PARAMS.include_neighbour;
  test.cI_resi_expand = expandCI(test.resi_adj_mat,test.cI_resi,k); %check the paramater k
  %-------------gather the density of the clustered groups-----------------
  [test.resi_allgrp_density,test.resi_global_density,test.resi_grp_density] = gatherGroupDensity(test.cI_resi_expand,test.resi_adj_mat);  
  test.cI_resi_expand = mapToResi(test.cI_resi_expand,Comp.num_seq);  %included after bug_report_4_03_2016
  test.cI      = changeToAtomSeq(test.cI_resi_expand,Comp.residue);
  
  %jump postponed till after expanssion
  test.cI_modified_org = test.cI;
  %jump = 25;%15;%25;%11;
  %test.cI_modified_org = changeToSeq(test.cI,jump); % see the jump param
    
  %-------------do the expansion multiple no. of times---------------------
     k                  = PARAMS.multi_expand_k2;
     size_change_cutoff = PARAMS.multi_expand_size_cutoff;
     grp_size_cutoff    = PARAMS.multi_expand_grp_min;
     atleast_incr       = PARAMS.multi_incr_min;
    [test.cI_modified1,test.multiexpand_count] = expandCI_multiple(test.adj_mat,test.cI_modified_org,k,size_change_cutoff,grp_size_cutoff,atleast_incr);
    
  %-------------do the jump step now---------------------------------------
      jump =  PARAMS.grp_expand;
      test.cI_modified = changeToSeq(test.cI_modified1,jump); % see the jump param
     
  test.patch_overlap_mat = getPatchOverlap(test.cI_modified);
  test.patch_adj = sparse(triu(test.patch_overlap_mat,1));  
  [test.dup_flag,test.dup_info] = findDuplicates(test.patch_overlap_mat);
  
  %---gather the density of the clustered groups------
  [test.atom_allgrp_density,test.atom_global_density,test.atom_grp_density] = gatherGroupDensity(test.cI_modified,test.adj_mat);
  
time.div_conq_prep = toc(div_conq_prep);
   %% loop through the cells in test.modified_cI, localize each
div_conq_loc_all = tic;
%[localize] = doLocalizeFragment(test.cI_modified,hydrogen_omission,rand_X,Comp.Cq, Comp.info, test.dup_flag, test.cI_modified_org, eq_cons_all, up_bounds, wh_vdw_bounds, vdw_bounds, sdp_lo_bounds, lo_bounds, wh_lo_bounds);

[localize] = doLocalizeFragment_v2(test.cI_modified,hydrogen_omission,rand_X,Comp.Cq, Comp.info, test.dup_flag, test.cI_modified_org, eq_cons_all, up_bounds, wh_vdw_bounds, vdw_bounds, sdp_lo_bounds, lo_bounds, wh_lo_bounds);
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
