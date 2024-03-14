%% ---------------------- Divide and Conquer ------------------------------
%   - Model as residue graph: 
%   * break the graph: 
%   * include neighbours:
%   * expand groups over multiple iterations:
%   * include small jumps:
%   - save pre localize data 
%%
   hydrogen_omission =  0; %PARAMS.hydrogen_omission;
   protein_name      =  INPUTS.protein_name;
   
  LOG.div_conquer = strcat(LOG.div_n_conquer,filesep,'div_n_conquer_',protein_name,'.txt');
  LOG.prev_txt = '';
%%  preprocess
div_conq_prep = tic;
  preprocess;

  % ---- write the raw_upls ----- %
  zz = raw_up{1};
  upl_name_raw = strcat(INPUTS.upl_file,'.raw_upl.upl');
  raw_upl_file = fopen(upl_name_raw{1},'w');
  for k = 1:length(zz)
     if k == 1
        fprintf(raw_upl_file,'%d\t%s\t%s\t%d\t%s\t%s\t%f',zz(k).sres,A(zz(k).stype).name,zz(k).satom,zz(k).tres,A(zz(k).ttype).name,zz(k).tatom,zz(k).dist);
     else
        fprintf(raw_upl_file,'\n%d\t%s\t%s\t%d\t%s\t%s\t%f',zz(k).sres,A(zz(k).stype).name,zz(k).satom,zz(k).tres,A(zz(k).ttype).name,zz(k).tatom,zz(k).dist);
     end
  end
  fclose(raw_upl_file);
  % ----------------------------- %
  
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
     %take_aug_bounds =  PARAMS.aug_bounds ;%8;%7.5;
     take_aug_bounds = 0; % testing
     if take_aug_bounds > 0
         up_bounds_org = up_bounds;   
         augment_bounds = augmentBoundsByTrieq(up_bounds(:,1:3),test.n_atoms);
         aug_I = find(augment_bounds(:,3)<take_aug_bounds);
         augment_bounds = [ augment_bounds(aug_I,:), zeros(length(aug_I),1)];
         up_bounds = [up_bounds; augment_bounds];
     end

     % augment up_bounds again
     up_bounds_before_multi_augbounds = up_bounds;
     %if isfield(PARAMS, 'aug_bounds_again_iter')
     %   for i=1:PARAMS.aug_bounds_again_iter
     %           augment_bounds = augmentBoundsByTrieq(up_bounds(:,1:3),test.n_atoms);
     %           aug_I = find(augment_bounds(:,3)<PARANS.aug_bounds_again_bounds);
     %           augment_bounds = [ augment_bounds(aug_I,:), zeros(length(aug_I),1)];
     %           up_bounds = [up_bounds; augment_bounds];
     %   end
     %end
  
 % for the adjacency matrix
 [test.resi_adj_mat,test.resi_adj_mat_eq,test.adj_mat,test.adj_mat_weq] = ...
      buildResiAdjMat(test.n_atoms,eq_cons_all,up_bounds,sdp_lo_bounds,Comp);   
  
 %% call the breaking algorithm

    [test.cI_test_resi, test.ret_densities] = preTestByDenseGraph((test.resi_adj_mat+test.resi_adj_mat')>0,3,0.1);%PARAMS.eta_lo);
