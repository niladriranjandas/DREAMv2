save_path="localize_tmp";

%%
num_groups       = length(grps);
 localize_X_noref = cell(1,num_groups);
 localize_time    = cell(1,num_groups);

 for i=1:length(grps)
      div_conq_loc_i = tic;
      disp('==============================================================');
      fprintf('\n-------------------------Grp-%d-------------------------------\n',grps(i));
      if test_dup_flag(grps(i))
        %=======================do SPROS========================================
        fprintf('\n_______Localization of grp-%d_____\n',grps(i));
        try
                %localize_X_noref{i} = doSPROSatomwise_v3(localize_atoms{grps(i)},...
                %                                         rand_X_grp{grps(i)},Comp_Cq_grp{grps(i)},...
                %                                         localize_eq_cons_all{grps(i)},localize_up_bounds{grps(i)},localize_wh_vdw_bounds{grps(i)},localize_vdw_bounds{grps(i)},localize_sdp_lo_bounds{grps(i)},localize_lo_bounds{grps(i)},localize_wh_lo_bounds{grps(i)},... 
                %                      hydrogen_omission); 
                localize_X_noref{i} = doSPROSatomwise_v4(localize_atoms{grps(i)},...
                                                         rand_X_grp{grps(i)},Comp_Cq_grp{grps(i)},...
                                                         localize_eq_cons_all{grps(i)},localize_up_bounds{grps(i)},localize_wh_vdw_bounds{grps(i)},localize_vdw_bounds{grps(i)},localize_sdp_lo_bounds{grps(i)},localize_lo_bounds{grps(i)},localize_wh_lo_bounds{grps(i)},... 
                                                         hydrogen_omission);                
        catch err       
                disp(err.message);       
                warning('Grp-%d could NOT be localized by SPROS',i);
       
                disp('---try to solve without SPROS-----')
          % localize a particular group
               display('-------begin cvx-------')               
             try
                 localize_method(i) = 2;
                 [X_tmp,eigens,slacks] = solve_by_cvx(length(localize_atoms{grps(i)}),localize_eq_cons_all{grps(i)},localize_up_bounds{grps(i)},localize_sdp_lo_bounds{grps(i)},[],[10,10,10,10,10]);
                
                 localize_X_noref{i} = X_tmp(1:3,:);
              catch err2
                 localize_method{i} = 0;
                 warning('Grp-%d CVX error. Could NOT be localized.',i)
              continue
             end       
        end   
        localize_time{i} = toc(div_conq_loc_i);
      else
        localize_method{i} = -1;         
        localize_X_noref{i} = [];
        localize_time{i} = -1;
      end
 end
    
 disp("===== now refine ====")
 localize_X_refine = cell(1,num_groups);
 localize_info     = cell(1,num_groups);
 for i=1:length(grps)
      %if test_dup_flag(grps(i))
      if test_dup_flag(grps(i)) && localize_method{i} > 0 %%              % ---- changed on 24-June-2021 solv_by_cvx gives empty --- %
        %-------------------------------------
        fprintf('\n________Refinement of grp-%d________\n',grps(i));
        [localize_X_refine{i},localize_info{i}] = postProcessingMe(localize_X_noref{i}, localize_eq_cons_all{grps(i)}, localize_lo_bounds{grps(i)}, full(localize_up_bounds{grps(i)}), [1e1,1,2,-0.001],[10,10,10,10,10]);
        disp('------violations after refinement-----');
      else
        continue;
        %localize_method{i} = -1; 
      end
 end         
  
 
 localize_include_index = cell(1,num_groups);
 localize_chk_coord_ref = cell(1,num_groups);
 for i=1:length(grps)
      %if test_dup_flag(grps(i))       
      if test_dup_flag(grps(i)) && localize_method{i} > 0 %%              % ---- changed on 24-June-2021 solv_by_cvx gives empty --- %
        %==============extract the original PDB coordinates=====================
        [localize_include_index{i},localize_chk_coord_ref{i}] = getRidOfPseudo(localize_X_refine{i}',localize_atoms{grps(i)},Comp_info);
        [include_index_noref,chk_coord_noref] = getRidOfPseudo(localize_X_noref{i}',localize_atoms{grps(i)},Comp_info);      
        %=======================================================================
        %localize_time{i} = toc(div_conq_loc_i);
      else
        continue;
        %localize_method{i} = -1; 
      end
 end


  %%
   for i=1:length(grps)
      localize(i).atoms         =  localize_atoms{grps(i)};
      localize(i).atoms_resi    = localize_atoms_resi{grps(i)};
      
      localize(i).method        = localize_method{i};
      localize(i).eq_cons_all   =  localize_eq_cons_all{grps(i)};
      localize(i).up_bounds     = localize_up_bounds{grps(i)};
      localize(i).wh_vdw_bounds = localize_wh_vdw_bounds{grps(i)};
      localize(i).vdw_bounds    = localize_vdw_bounds{grps(i)};
      localize(i).sdp_lo_bounds = localize_sdp_lo_bounds{grps(i)};
      localize(i).lo_bounds     = localize_lo_bounds{grps(i)};
      localize(i).wh_lo_bounds  =  localize_wh_lo_bounds{grps(i)};
            
      localize(i).X_noref       = localize_X_noref{i};
      localize(i).time          = localize_time{i};
      
      localize(i).X_refine      = localize_X_refine{i};
      localize(i).info          = localize_info{i};
      
      localize(i).include_index = localize_include_index{i};
      localize(i).chk_coord_ref = localize_chk_coord_ref{i};
      
   end

%%

