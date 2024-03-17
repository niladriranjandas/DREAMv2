%%
%Copyright (c) 2024 niladriranjandas
%make the localize data structure after substructures has been modeled from fragments separately
%%

function localize = assimilate_after_seplocalize(oppfile, test_dup_flag)

%%
    for j=1:length(test_dup_flag)
      % if ~(test_dup_flag)
            localize(j).atoms         = [];
            localize(j).atoms_resi    = [];
            localize(j).method        = -1;
            localize(j).eq_cons_all   = [];
            localize(j).up_bounds     = [];
            localize(j).wh_vdw_bounds = [];
            localize(j).vdw_bounds    = [];
            localize(j).sdp_lo_bounds = [];
            localize(j).lo_bounds     = [];
            localize(j).wh_lo_bounds  = [];
    
            localize(j).X_noref       = [];
            localize(j).time          = [];
            localize(j).X_refine      = [];
            localize(j).info          = [];
            localize(j).include_index = [];
            localize(j).chk_coord_ref = []; 
      % end        
    end
 
%%
	fid = fopen(oppfile,'rt');
	while true
	  thisline = fgetl(fid);
	  if ~ischar(thisline)
          break; 
      end  %end of file
	  localize_grp_i = load(thisline,'localize','grps');
	
	  for j = 1:length(localize_grp_i.localize)
			localize(localize_grp_i.grps(j)).atoms         = localize_grp_i.localize(j).atoms ;       
			localize(localize_grp_i.grps(j)).atoms_resi    = localize_grp_i.localize(j).atoms_resi;   
			localize(localize_grp_i.grps(j)).method        = localize_grp_i.localize(j).method;       
			localize(localize_grp_i.grps(j)).eq_cons_all   = localize_grp_i.localize(j).eq_cons_all ; 
			localize(localize_grp_i.grps(j)).up_bounds     = localize_grp_i.localize(j).up_bounds ;   
			localize(localize_grp_i.grps(j)).wh_vdw_bounds = localize_grp_i.localize(j).wh_vdw_bounds;
			localize(localize_grp_i.grps(j)).vdw_bounds    = localize_grp_i.localize(j).vdw_bounds ;  
			localize(localize_grp_i.grps(j)).sdp_lo_bounds = localize_grp_i.localize(j).sdp_lo_bounds;
			localize(localize_grp_i.grps(j)).lo_bounds     = localize_grp_i.localize(j).lo_bounds  ;  
			localize(localize_grp_i.grps(j)).wh_lo_bounds  = localize_grp_i.localize(j).wh_lo_bounds ;
	
			localize(localize_grp_i.grps(j)).X_noref       = localize_grp_i.localize(j).X_noref  ;    
			localize(localize_grp_i.grps(j)).time          = localize_grp_i.localize(j).time  ;       
			localize(localize_grp_i.grps(j)).X_refine      = localize_grp_i.localize(j).X_refine ;    
			localize(localize_grp_i.grps(j)).info          = localize_grp_i.localize(j).info ;        
			localize(localize_grp_i.grps(j)).include_index = localize_grp_i.localize(j).include_index;
			localize(localize_grp_i.grps(j)).chk_coord_ref = localize_grp_i.localize(j).chk_coord_ref;	
	  end
  	end

  	fclose(fid);
    
    %% check
    for j=1:length(test_dup_flag)
       if test_dup_flag
          if isempty(localize(j).X_refine)
             error("ERROR: module: assimilate_after_seplocalize: Grp-%d is empty.",j);
          end
       end
    end
    
    %% write the pdb files 
end
