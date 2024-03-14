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