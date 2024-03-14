%time.div_conq_prep = toc(div_conq_prep);
   %% loop through the cells in test.modified_cI, localize each
div_conq_loc_all = tic;
  for i=1:length(test.cI_modified)
   div_conq_loc_i = tic;
   display('==============================================================');
   fprintf('\n-------------------------Grp-%d-------------------------------\n',i);
      localize(i).atoms = test.cI_modified{i};
      localize(i).atoms_resi = test.cI_modified_org{i};
   if test.dup_flag(i)
      localize(i).method = 1;
      
   %=====extract the constraints===========================================
      localize(i).eq_cons_all   = extractCons(eq_cons_all,localize(i).atoms);
      eq_cons_all_resi          = extractCons(eq_cons_all,localize(i).atoms_resi);
      localize(i).up_bounds     = extractCons(up_bounds,localize(i).atoms);
      up_bounds_resi            = extractCons(up_bounds,localize(i).atoms_resi);
      localize(i).wh_vdw_bounds = extractCons(wh_vdw_bounds,localize(i).atoms); 
      %for H atoms needs to be mapped
      localize(i).vdw_bounds    = extractCons(vdw_bounds,localize(i).atoms);
      vdw_bounds_resi           = extractCons(vdw_bounds,localize(i).atoms_resi);
      localize(i).sdp_lo_bounds = extractCons(sdp_lo_bounds,localize(i).atoms);
      sdp_lo_bounds_resi        = extractCons(sdp_lo_bounds,localize(i).atoms_resi);
    
      localize(i).lo_bounds     = extractCons(lo_bounds,localize(i).atoms);
      lo_bounds_resi            = extractCons(lo_bounds,localize(i).atoms_resi);
      localize(i).wh_lo_bounds  = extractCons(wh_lo_bounds,localize(i).atoms); 
      %for with H atoms needs to be mapped
   %=======================================================================
      
   %=======================do SPROS========================================
   fprintf('\n_______Localization of grp-%d_____\n',i);
        try
                localize(i).X_noref = doSPROSatomwise_v2(localize(i).atoms,...
                                      rand_X,Comp.Cq,...
                                      localize(i).eq_cons_all,localize(i).up_bounds,localize(i).wh_vdw_bounds,localize(i).vdw_bounds,localize(i).sdp_lo_bounds,localize(i).lo_bounds,localize(i).wh_lo_bounds,... 
                                      hydrogen_omission);     
        catch err       
                disp(err.message);       
                warning('Grp-%d could NOT be localized by SPROS',i);
       
                disp('---try to solve without SPROS-----')
          % localize a particular group
               display('-------begin cvx-------')               
             try
                 localize(i).method = 2;
                 [X_tmp,eigens,slacks] = solve_by_cvx(length(localize(i).atoms),localize(i).eq_cons_all,localize(i).up_bounds,localize(i).sdp_lo_bounds,[],[10,10,10,10,10]);
                
                 localize(i).X_noref = X_tmp(1:3,:);
              catch err2
                 localize(i).method = 0;
                 warning('Grp-%d CVX error. Could NOT be localized.',i)
              continue
             end       
        end   
      
        %=======================refinement======================================
        %------violation before refinement----
        disp('------violations before refinement--');     
        if localize(i).method > 0                                           % ---- changed on 24-June-2021 solv_by_cvx gives empty --- %
            tmp_atm_index  = find(ismember(localize(i).atoms, localize(i).atoms_resi));
            X_noref_resi = localize(i).X_noref(:,tmp_atm_index);
            %[localize(i).CompSeq,localize(i).p_report] = ...
            %       doChecknCorrect(localize(i).atoms_resi,Comp,X_noref_resi,lo_bounds_resi,up_bounds_resi,eq_cons_all_resi,1);
            %[localize(i).CompSeq,localize(i).p_report] = ...
            %     doChecknCorrect(localize(i).atoms_resi,Comp,localize(i).X_noref,localize(i).lo_bounds,localize(i).up_bounds,localize(i).eq_cons_all,1);
     
            %-------------------------------------
            fprintf('\n________Refinement of grp-%d________\n',i);
      
            %[localize(i).X_refine,localize(i).info] = postProcessingMe(localize(i).X_noref, localize(i).eq_cons_all, localize(i).lo_bounds, localize(i).up_bounds, [1e1,1,2,-0.001],[10,10,10,10,10]);
            [localize(i).X_refine,localize(i).info] = postProcessingMe(localize(i).X_noref, localize(i).eq_cons_all, localize(i).lo_bounds, full(localize(i).up_bounds), [1e1,1,2,-0.001],[10,10,10,10,10]);
       
            disp('------violations after refinement-----');
            X_ref_resi = localize(i).X_refine(:,tmp_atm_index);
         
            %==============extract the original PDB coordinates=====================
            [localize(i).include_index,localize(i).chk_coord_ref] = getRidOfPseudo(localize(i).X_refine',localize(i).atoms,Comp.info);
            [include_index_noref,chk_coord_noref] = getRidOfPseudo(localize(i).X_noref',localize(i).atoms,Comp.info);      
           %localize(i).include_index = include_index;
        else                                                               % ---- changed on 24-June-2021 solv_by_cvx gives empty --- %     
            localize(i).X_refine      = [];
            localize(i).info          = [];
            localize(i).include_index = [];
            localize(i).chk_coord_ref = [];           
        end
   %=======================================================================
    localize(i).time = toc(div_conq_loc_i);
   else
      localize(i).method = -1; 
   end
  end
  
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
