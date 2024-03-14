%% Use modeller result as start point do anchored localization
%   
%%
 clearvars fillQ Aloc fillQ_chiral
 
 anchor_loc_prep = tic;
 %fillQ.pdbStruct = pdbread('modeller/current_run_fixed_anchors/grp_5_H.pdb');  % output of modeller
 fillQ.pdbStruct = pdbread(OUTPUT.modeller_n_H);
   LOG.anchor = strcat(LOG.div_n_conquer,filesep,'anchored_loc.txt');
   LOG.prev_txt = '';
   
    local_folder_plots_postreg = OUTPUT.local_folder_plots_postreg;

 %%
fillQ.Atoms     = fillQ.pdbStruct.Model.Atom;

fillQ.X_all = [[fillQ.Atoms.X]; [fillQ.Atoms.Y];[fillQ.Atoms.Z]];
fillQ.AtomName_all =  {fillQ.Atoms.AtomName};
fillQ.residue_all  = [fillQ.Atoms.resSeq];
fillQ.residuename_all = {fillQ.Atoms.resName};
       
    if fillQ.residue_all(1) ~= org_num(1)
       error('Module doAnchoredLoc: input PDB residue-1 and seq file residue-1 mismatch'); 
    end

fillQ.indx_begin = find(fillQ.residue_all==min_res); fillQ.indx_begin = fillQ.indx_begin(1);
fillQ.indx_end   = find(fillQ.residue_all==max_res); fillQ.indx_end   = fillQ.indx_end(end);

% fillQ.X        = fillQ.X_all(:,fillQ.indx_begin:end);
% fillQ.AtomName = fillQ.AtomName_all(fillQ.indx_begin:end);
% fillQ.residue  = fillQ.residue_all(fillQ.indx_begin:end);
% fillQ.residuename = fillQ.residuename_all(fillQ.indx_begin:end);

fillQ.X        = fillQ.X_all(:,fillQ.indx_begin:fillQ.indx_end);
fillQ.AtomName = fillQ.AtomName_all(fillQ.indx_begin:fillQ.indx_end);
fillQ.residue  = fillQ.residue_all(fillQ.indx_begin:fillQ.indx_end);
fillQ.residuename = fillQ.residuename_all(fillQ.indx_begin:fillQ.indx_end);
%------- convert pdb to our algo atom compatible atom name seq ------------
fillQ.x_filled = fillInPseudo(fillQ.X, fillQ.residue, fillQ.AtomName, fillQ.residuename,wh_Comp, A);

%-------- do chirality chk, correction and refine ----------------------
 disp('------------chirality check ----------')
fillQ.x_filled_old = fillQ.x_filled;
fillQ.MAX_chiral_itter = 20;
fillQ_chiral1(1).chiral_info = nan;

  for i=2:fillQ.MAX_chiral_itter % 20 will be assigned a param say max_itter
    disp(i)
    [fillQ.x_filled, ~, fillQ_chiral1(i).chiral_info] = chiralChknCorr_v2(fillQ.x_filled, 1:length(fillQ.x_filled_old), wh_up_bounds, wh_lo_bounds, eq_cons_all , wh_Comp, 1);
    if (i>3) && ...
       (fillQ_chiral1(i).chiral_info.err_main+fillQ_chiral1(i).chiral_info.err_side == fillQ_chiral1(i-2).chiral_info.err_main + fillQ_chiral1(i-2).chiral_info.err_side  )
           break
    end
  end
  
 if (fillQ_chiral1(i).chiral_info.err_main+fillQ_chiral1(i).chiral_info.err_side) ~= 0
     [fillQ.x_filled, ~, fillQ_chiral1(i+1).chiral_info] = chiralChknCorr_v2(fillQ.x_filled, 1:length(fillQ.x_filled_old), wh_up_bounds, wh_lo_bounds, eq_cons_all , wh_Comp, 0);
 end
%[fillQ.x_filled, fillQ.chiral_chk_info] = chiralChknCorr(fillQ.x_filled_old, 1:length(fillQ.x_filled_old), wh_up_bounds, wh_lo_bounds, eq_cons_all , wh_Comp);

%% anchored localization
 fillQ.atom_map = find(~isnan(fillQ.x_filled(1,:)));
  tmp  = find(wh_Comp.info(5,:)==0);                 %pseudo atoms letf out
 fillQ.atom_map = setdiff(fillQ.atom_map, tmp);      %pseudo atoms letf out coz optimization fails (unbound) otherwise
 %fillQ.atom_map = find(wh_Comp.info(5,:)~=0);

Aloc.atom_map     = fillQ.atom_map;
Aloc.Xref_atommap = fillQ.x_filled(:,fillQ.atom_map);

 %-------------------------- anchors and gap ------------------------------
 [Aloc.back_bone_tmp_index, Aloc.back_bone_atom_map] = onlyTheBackBone(reg.chk_coord_ref',reg.include_index,Comp);
 Aloc.anchors_resi = unique(Comp.residue(Aloc.back_bone_atom_map));
        [zz2.back_bone_tmp_index, zz2.back_bone_atom_map] = onlyTheBackBone(reg2.chk_coord_ref',reg2.include_index,Comp);
        Aloc.anchors_resi = union(Aloc.anchors_resi,unique(Comp.residue(zz2.back_bone_atom_map)));
 Aloc.gaps_resi = setdiff(Comp.num_seq, Aloc.anchors_resi);
 
 %--------------------------anchors and vars index ------------------------
 Aloc.wh_Comp_residue_atommap = wh_Comp.residue(Aloc.atom_map);
  Aloc.X_indx = [];
  for i=1:length(Aloc.gaps_resi)
      Aloc.atom_index = find(Aloc.wh_Comp_residue_atommap==Aloc.gaps_resi(i));       
      Aloc.X_indx = [Aloc.X_indx, Aloc.atom_index'];
  end   
  Aloc.anchors_indx = [];
  for i=1:length(Aloc.anchors_resi)
      Aloc.atom_index = find(Aloc.wh_Comp_residue_atommap==Aloc.anchors_resi(i));
      Aloc.anchors_indx = [Aloc.anchors_indx, Aloc.atom_index'];
  end 
 
 Aloc.anchors = Aloc.Xref_atommap(:,Aloc.anchors_indx);
 Aloc.X       = Aloc.Xref_atommap(:,Aloc.X_indx);
  
 %-------------------------- extract constraints---------------------------
 % Aloc.eq_cons_atommap      = extractCons(equality_con_former(wh_rand_X,wh_Comp,2), Aloc.atom_map);
  Aloc.eq_cons_atommap      = extractCons(wh_eq_cons_all, Aloc.atom_map);
  Aloc.wh_up_bounds_atommap = extractCons(wh_up_bounds, Aloc.atom_map);
  Aloc.wh_lo_bounds_atommap = extractCons(wh_lo_bounds, Aloc.atom_map);
 
  
 %---------- get bounds for (vars, anchors) and (vars, vars) --------------
%  [Aloc.eq, Aloc.up, Aloc.lo] = getBoundsForAnchor(Aloc.X_indx, Aloc.anchors_indx, ...
%                                                  Aloc.eq_cons_atommap,...
%                                                  Aloc.wh_up_bounds_atommap,...
%                                                  Aloc.wh_lo_bounds_atommap);
 [Aloc.eq, Aloc.up, Aloc.lo] = getBoundsForAnchor_v2(Aloc.X_indx, Aloc.anchors_indx, ...
                                                     Aloc.eq_cons_atommap,...
                                                     Aloc.wh_up_bounds_atommap,...
                                                     Aloc.wh_lo_bounds_atommap);
                                             
 Aloc.up.X = Aloc.up.X(:,1:3);
 Aloc.up.A = Aloc.up.A(:,1:3);

 Aloc.lo.X = Aloc.lo.X(:,1:3);
 Aloc.lo.A = Aloc.lo.A(:,1:3);
 
 time.anchor_loc_prep = toc(anchor_loc_prep);
 %========================anchored localization============================
        disp('----------start anchored minimization----------')
        anchore_loc = tic;
        [Aloc.localize.Xsolv, Aloc.localize.info] = anchoredLocalize(Aloc.Xref_atommap, Aloc.X_indx, Aloc.anchors_indx,....
                                                         Aloc.eq, Aloc.lo, Aloc.up,...
                                                         [2,2,1,1,1,1,-1], f);
        Aloc.localize.X(:,Aloc.anchors_indx) = Aloc.anchors;
        Aloc.localize.X(:,Aloc.X_indx)       = reshape(Aloc.localize.Xsolv,3,length(Aloc.X_indx));
        time.anchore_loc = toc(anchore_loc);
        %-----------------------------logs---------------------------------
         fprintf('\n=================================== Violations (Stage: Anchored localization) =========================================');
         LOG.curr_txt = sprintf('\n=================================== Violations (Stage: Anchored localization) =========================================');
         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
         
         fprintf('\n Before refinement');
         LOG.curr_txt = sprintf('\n Before refinement');
         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
         
         [violation.mean_err, violation.min_err, violation.max_err, violation.err_percent] = calcViolations(Aloc.localize.X,Aloc.eq_cons_atommap ,3);
         fprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
         LOG.curr_txt = sprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
        
         [violation.mean_err, ~, violation.max_err, violation.err_percent] = calcViolations(Aloc.localize.X,Aloc.wh_up_bounds_atommap,1);
         fprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
         LOG.curr_txt = sprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);

         [violation.mean_err, violation.min_err, ~, violation.err_percent] = calcViolations(Aloc.localize.X,Aloc.wh_lo_bounds_atommap,2);
         fprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
         LOG.curr_txt = sprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
        %------------------------------------------------------------------
        tmp=cell2mat(Aloc.localize.info.fevalrec); figure;plot(1:length(tmp),tmp);hold on; title('anchored'); hold off;
 %=========================================================================
 %---------------------------write pdb ------------------------------------
 Aloc.filename  = sprintf('%s_%d_ancloc_nochiral.pdb',protein_name,OUTPUT.model_count);
 Aloc.file_full = strcat(OUTPUT.local_folder_pdb, filesep, Aloc.filename);
 {writeToPDB(Aloc.file_full,Aloc.atom_map,Aloc.localize.X,wh_Comp)};
 LOG.curr_txt = sprintf('\n Written file %s \n',Aloc.file_full);
 LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
 
 %---------------------------write ramachandran map -----------------------
 [Aloc.localize.rama_out,fig_handler] = ramachandranMe(Aloc.file_full,'Regions',true,'Glycine',true);
 plot_filename_rama = sprintf('%s_rama_anchor_noref_%d.fig',protein_name,OUTPUT.model_count);
 plot_file_rama    = strcat(local_folder_plots_postreg,filesep,plot_filename_rama);
 savefig(fig_handler,plot_file_rama);
 LOG.curr_txt = sprintf('\n Written file %s \n',plot_file_rama);
 LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);

%% refine

    anchor_loc_ref = tic;
    %--------------------get the constraints--------------------------------
    fillQ.wh_eq_cons_all = Aloc.eq_cons_atommap;%extractCons(wh_eq_cons_all, fillQ.atom_map);
    fillQ.wh_up_bounds   = Aloc.wh_up_bounds_atommap;%extractCons(wh_up_bounds,   fillQ.atom_map);  
    %fillQ.wh_vdw_bounds  = extractCons(wh_vdw_bounds,  fillQ.atom_map);   
    fillQ.wh_lo_bounds   = Aloc.wh_lo_bounds_atommap;%extractCons(wh_lo_bounds,   fillQ.atom_map);
        
    %--------------------call the refine module-----------------------------
    fprintf('\n__________refine module___\n');
%     [fillQ.X_refine,fillQ.info] = postProcessingMe(Aloc.localize.X,...
%                                                    fillQ.wh_eq_cons_all, fillQ.wh_lo_bounds, fillQ.wh_up_bounds,...
%                                                    [1e4,1,1,-0.001],[10,10,10,10,10]);
    [fillQ.X_refine,fillQ.info] = postProcessingMe(Aloc.localize.X,...
                                                   fillQ.wh_eq_cons_all, fillQ.wh_lo_bounds, full(fillQ.wh_up_bounds),...
                                                   [1e4,1,1,-0.001],[10,10,10,10,10]);
    time.anchor_loc_ref = toc(anchor_loc_ref);                                           
    %----------------------------logs--------------------------------------
        fprintf('\n\n After refinement');
        LOG.curr_txt = sprintf('\n\n After refinement');
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
        
        [violation.mean_err, violation.min_err, violation.max_err, violation.err_percent] = calcViolations(fillQ.X_refine,fillQ.wh_eq_cons_all,3);
        fprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
        
        [violation.mean_err, ~, violation.max_err, violation.err_percent] = calcViolations(fillQ.X_refine,fillQ.wh_up_bounds,1);
        fprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);

        [violation.mean_err, violation.min_err, ~, violation.err_percent] = calcViolations(fillQ.X_refine,fillQ.wh_lo_bounds,2);
        fprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
    %----------------------------------------------------------------------
  
   fillQ.filename  = sprintf('%s_%d_ancloc_nochiral_refine.pdb',protein_name,OUTPUT.model_count);
   fillQ.file_full = strcat(OUTPUT.local_folder_pdb, filesep, fillQ.filename);                                               
   {writeToPDB(fillQ.file_full,fillQ.atom_map,fillQ.X_refine,wh_Comp)};
   LOG.curr_txt = sprintf('\n Written file %s \n',Aloc.file_full);
   LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
   
   [fillQ.rama_out,fig_handler] = ramachandranMe(fillQ.file_full,'Regions',true,'Glycine',true);
   plot_filename_rama = sprintf('%s_rama_anchor_ref_%d.fig',protein_name,OUTPUT.model_count);
   plot_file_rama    = strcat(local_folder_plots_postreg,filesep,plot_filename_rama);
   savefig(fig_handler,plot_file_rama);
   LOG.curr_txt = sprintf('\n Written file %s \n',plot_file_rama);
   LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
   
   figure;plot(1:length(fillQ.info.obj),fillQ.info.obj);hold on; title('refine');hold off

   clearvars fig_handler
   
   LOG.curr_txt = sprintf('\n Module: Anchored Localization: prep time: %d \t anchored localization: %d \t refinement: %d \t',time.anchor_loc_prep,time.anchore_loc,time.anchor_loc_ref);
   LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
   
   LOG.prev_txt = mexDolog(LOG.prev_txt,'\n END of anchored localization of groups',1,LOG.anchor);
   
   %% ------------- chirality chk and fix ----------------------
   [fillQ.x_ref_chiral, fillQ.x_ref_chiral_info] = chiralChknCorr(fillQ.X_refine, fillQ.atom_map, fillQ.wh_up_bounds, fillQ.wh_lo_bounds, fillQ.wh_eq_cons_all , wh_Comp);
   disp('------------chirality check ----------')
   %fillQ.x_ref_chiral = fillQ.X_refine;
   fillQ.MAX_chiral_itter = 20;
   fillQ_chiral2(i).chiral_info = nan;
      
   for i=2:fillQ.MAX_chiral_itter % 20 will be assigned a param say max_itter
     disp(i)
     [fillQ.x_ref_chiral, ~, fillQ_chiral2(i).chiral_info] = chiralChknCorr_v2(fillQ.x_ref_chiral, fillQ.atom_map, fillQ.wh_up_bounds, fillQ.wh_lo_bounds, fillQ.wh_eq_cons_all , wh_Comp, 1);
     if (i>3) && ...
         (fillQ_chiral2(i).chiral_info.err_main+fillQ_chiral2(i).chiral_info.err_side == fillQ_chiral2(i-2).chiral_info.err_main + fillQ_chiral2(i-2).chiral_info.err_side  )
            break
     end
   end
  
 if (fillQ_chiral2(i).chiral_info.err_main+fillQ_chiral2(i).chiral_info.err_side) ~= 0
     [fillQ.x_ref_chiral, ~, fillQ_chiral2(i+1).chiral_info] = chiralChknCorr_v2(fillQ.x_ref_chiral, fillQ.atom_map, fillQ.wh_up_bounds, fillQ.wh_lo_bounds, fillQ.wh_eq_cons_all , wh_Comp, 0);
 end
  %------------------------------------------------------------------------
  
   fillQ.filename  = sprintf('%s_%d_ancloc_refine.pdb',protein_name,OUTPUT.model_count);
   fillQ.file_full = strcat(OUTPUT.local_folder_pdb, filesep, fillQ.filename);                                               
   {writeToPDB(fillQ.file_full,fillQ.atom_map,fillQ.x_ref_chiral,wh_Comp)};
   
   %----------------------------logs--------------------------------------
        fprintf('\n\n After chiral fix and refinement');
        LOG.curr_txt = sprintf('\n\n After chiral fix and refinement');
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
        
        [violation.mean_err, violation.min_err, violation.max_err, violation.err_percent] = calcViolations(fillQ.x_ref_chiral,fillQ.wh_eq_cons_all,3);
        fprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
        
        [violation.mean_err, ~, violation.max_err, violation.err_percent] = calcViolations(fillQ.x_ref_chiral,fillQ.wh_up_bounds,1);
        fprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);

        [violation.mean_err, violation.min_err, ~, violation.err_percent] = calcViolations(fillQ.x_ref_chiral,fillQ.wh_lo_bounds,2);
        fprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
        
LOG.prev_txt = mexDolog(LOG.prev_txt,'\n END of chiral corr post anchored localization of groups',1,LOG.anchor);        