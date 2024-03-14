function [localize] = doLocalizeFragment(test_cI_modified,hydrogen_omission,rand_X,Comp_Cq, Comp_info, test_dup_flag, test_cI_modified_org, eq_cons_all, up_bounds, wh_vdw_bounds, vdw_bounds, sdp_lo_bounds, lo_bounds, wh_lo_bounds)
%%
%
%%
 for i=1:length(test_cI_modified)
   localize_atoms{i} = test_cI_modified{i};
   localize_atoms_resi{i} = test_cI_modified_org{i}; 
   
      if test_dup_flag(i)
          localize_method{i} = 1;
      
    %=====extract the constraints===========================================
          localize_eq_cons_all{i}   = extractCons(eq_cons_all,localize_atoms{i});
          eq_cons_all_resi          = extractCons(eq_cons_all,localize_atoms_resi{i});
          localize_up_bounds{i}     = extractCons(up_bounds,localize_atoms{i});
          up_bounds_resi            = extractCons(up_bounds,localize_atoms_resi{i});
          localize_wh_vdw_bounds{i} = extractCons(wh_vdw_bounds,localize_atoms{i}); 
          %for H atoms needs to be mapped
          localize_vdw_bounds{i}    = extractCons(vdw_bounds,localize_atoms{i});
          vdw_bounds_resi           = extractCons(vdw_bounds,localize_atoms_resi{i});
          localize_sdp_lo_bounds{i} = extractCons(sdp_lo_bounds,localize_atoms{i});
          sdp_lo_bounds_resi        = extractCons(sdp_lo_bounds,localize_atoms_resi{i});
    
          localize_lo_bounds{i}     = extractCons(lo_bounds,localize_atoms{i});
          lo_bounds_resi            = extractCons(lo_bounds,localize_atoms_resi{i});
          localize_wh_lo_bounds{i}  = extractCons(wh_lo_bounds,localize_atoms{i}); 

          
      else
          localize_method{i} = -1;
          
          localize_eq_cons_all{i}   = [];
          eq_cons_all_resi          = [];
          localize_up_bounds{i}     = [];
          up_bounds_resi            = [];
          localize_wh_vdw_bounds{i} = []; 
          %for H atoms needs to be mapped
          localize_vdw_bounds{i}    = [];
          vdw_bounds_resi           = [];
          localize_sdp_lo_bounds{i} = [];
          sdp_lo_bounds_resi        = [];
    
          localize_lo_bounds{i}     = [];
          lo_bounds_resi            = [];
          localize_wh_lo_bounds{i}  = [];          
      end
   
 end
 disp("-- extraction done--")

 %%
 num_groups       = length(localize_method);
 localize_X_noref = cell(1,num_groups);
 localize_time    = cell(1,num_groups);
 
 parfor i=1:length(test_cI_modified)
        div_conq_loc_i = tic;
        disp('==============================================================');
        fprintf('\n-------------------------Grp-%d-------------------------------\n',i);
      if test_dup_flag(i)
        %=======================do SPROS========================================
        fprintf('\n_______Localization of grp-%d_____\n',i);
        try
                localize_X_noref{i} = doSPROSatomwise_v2(localize_atoms{i},...
                                      rand_X,Comp_Cq,...
                                      localize_eq_cons_all{i},localize_up_bounds{i},localize_wh_vdw_bounds{i},localize_vdw_bounds{i},localize_sdp_lo_bounds{i},localize_lo_bounds{i},localize_wh_lo_bounds{i},... 
                                      hydrogen_omission);     
        catch err       
                disp(err.message);       
                warning('Grp-%d could NOT be localized by SPROS',i);
       
                disp('---try to solve without SPROS-----')
          % localize a particular group
               display('-------begin cvx-------')               
             try
                 localize_method(i) = 2;
                 [X_tmp,eigens,slacks] = solve_by_cvx(length(localize_atoms{i}),localize_eq_cons_all{i},localize_up_bounds{i},localize_sdp_lo_bounds{i},[],[10,10,10,10,10]);
                
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
 parfor i=1:length(test_cI_modified)
      if test_dup_flag(i)
        %=======================refinement======================================
        %------violation before refinement----
        %disp('------violations before refinement--');     
        %tmp_atm_index  = find(ismember(localize_atoms(i), localize_atoms_resi(i)));
        %localize_X_noref_i = localize_X_noref(i)
        %X_noref_resi = localize_X_noref_i(:,tmp_atm_index);
        %[localize(i).CompSeq,localize(i).p_report] = ...
        %       doChecknCorrect(localize_atoms_resi(i),Comp,X_noref_resi,lo_bounds_resi,up_bounds_resi,eq_cons_all_resi,1);
        %[localize(i).CompSeq,localize(i).p_report] = ...
        %     doChecknCorrect(localize_atoms_resi(i),Comp,localize_X_noref(i),localize_lo_bounds(i),localize_up_bounds(i),localize_eq_cons_all(i),1);
     
        %-------------------------------------
        fprintf('\n________Refinement of grp-%d________\n',i);
      
        %[localize_X_refine(i),localize_info(i)] = postProcessingMe(localize_X_noref(i), localize_eq_cons_all(i), localize_lo_bounds(i), localize_up_bounds(i), [1e1,1,2,-0.001],[10,10,10,10,10]);
        [localize_X_refine{i},localize_info{i}] = postProcessingMe(localize_X_noref{i}, localize_eq_cons_all{i}, localize_lo_bounds{i}, full(localize_up_bounds{i}), [1e1,1,2,-0.001],[10,10,10,10,10]);
       
        disp('------violations after refinement-----');
        %localize_X_refine_i = localize_X_refine{i};
        %X_ref_resi = localize_X_refine_i(:,tmp_atm_index);
      else
        continue;
        %localize_method{i} = -1; 
      end
 end         
  
 localize_include_index = cell(1,num_groups);
 localize_chk_coord_ref = cell(1,num_groups);
 for i=1:length(test_cI_modified)
      if test_dup_flag(i)
        
        %==============extract the original PDB coordinates=====================
        [localize_include_index{i},localize_chk_coord_ref{i}] = getRidOfPseudo(localize_X_refine{i}',localize_atoms{i},Comp_info);
        [include_index_noref,chk_coord_noref] = getRidOfPseudo(localize_X_noref{i}',localize_atoms{i},Comp_info);      
        %localize_include_index(i) = include_index;
        %=======================================================================
        %localize_time{i} = toc(div_conq_loc_i);
      else
        continue;
        %localize_method{i} = -1; 
      end
 end


  %%
   for i=1:length(test_cI_modified)
      localize(i).atoms      =  localize_atoms{i};
      localize(i).atoms_resi = localize_atoms_resi{i};
      
      localize(i).method        = localize_method{i};
      localize(i).eq_cons_all   =  localize_eq_cons_all{i};
      localize(i).up_bounds     = localize_up_bounds{i};
      localize(i).wh_vdw_bounds = localize_wh_vdw_bounds{i};
      localize(i).vdw_bounds    = localize_vdw_bounds{i};
      localize(i).sdp_lo_bounds = localize_sdp_lo_bounds{i};
      localize(i).lo_bounds     = localize_lo_bounds{i};
      localize(i).wh_lo_bounds  =  localize_wh_lo_bounds{i};
            
      localize(i).X_noref       = localize_X_noref{i};
      localize(i).time          = localize_time{i};
      
      localize(i).X_refine      = localize_X_refine{i};
      localize(i).info          = localize_info{i};
      
      localize(i).include_index = localize_include_index{i};
      localize(i).chk_coord_ref = localize_chk_coord_ref{i};
      
   end

end
