function [filestorun, oppfile, filename_grps] = doLocalizeFragment_sepprocess_pre(protname,test_cI_modified,hydrogen_omission,rand_X,Comp_Cq, Comp_info, test_dup_flag, test_cI_modified_org, eq_cons_all, up_bounds, wh_vdw_bounds, vdw_bounds, sdp_lo_bounds, lo_bounds, wh_lo_bounds)
%%
%
%%

 NUM_PROCESSES = 3;
 save_path = sprintf("localize_tmp/%s",protname);

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

	      rand_X_grp{i}                 = rand_X(:,localize_atoms{i});       
	  %% for Comp_Cq{i}
          [I,J,~] = find(Comp_Cq(localize_atoms{i},:)); % [I,J] are sorted in ascending order J wise
          cq_grp = unique(J);
          J_tmp = zeros(length(J),1);
          for ii=1:length(J)
                J_tmp(ii) = find(cq_grp==J(ii));
          end
          Comp_Cq_grp{i} = sparse(I,J_tmp,1,length(localize_atoms{i}),length(cq_grp)); 
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
		
	  rand_X_grp{i}             = [];
	  Comp_Cq_grp{i}            = [];     
      end
   
 end
 disp("-- extraction done--")

%%
 filename = sprintf("%s/%s_prelocalize.mat",save_path,protname);
 oppfile = sprintf("%s/%s_after_localize.txt",save_path,protname);
 filestorun = sprintf("%s/%s_run_to_localize.txt",save_path,protname);
 save(filename);
 
%% create separate scripts %%
 localize_grp_sizes = nan(1,length(localize_atoms));
 for kk = 1:length(localize_atoms)
	if test_dup_flag(kk)
		localize_grp_sizes(kk) = length(localize_atoms{kk});
	end
 end
 
 grps = divideIntoGrps(localize_grp_sizes, localize_atoms, test_dup_flag, NUM_PROCESSES);

 opp_fid = fopen(oppfile,'w');
 run_fid = fopen(filestorun,'w');
 if opp_fid > 0 && run_fid > 0
    for i=1:length(grps)
        filename_grps = sprintf("%s/run_%s_part_%d.m",save_path,protname,i);
        savefile_grps = sprintf("%s/postloc_%s_grp_%d.mat",save_path,protname,i);
	
        fprintf(opp_fid,"%s\n",savefile_grps);
        fprintf(run_fid,"%s\n",filename_grps);
    
    	fid = fopen(filename_grps,'w');
        if (fid > 0)
            fprintf(fid,"\ngrps = [");
            write_grps = grps{i};
            for j = 1:length(grps{i})			
    			if j == length(grps{i})
        			fprintf(fid,"%d];\n",write_grps(j));
                else
                	fprintf(fid,"%d,",write_grps(j));
                end
            end
            fprintf(fid,"\ncd ..");
            fprintf(fid,"\naddpath(genpath(pwd));");
            fprintf(fid,"\ncd code;");
            fprintf(fid,"\nload %s",filename);
            fprintf(fid,"\ndoLocalizeFragment_sepprocess");
            fprintf(fid,"\nsave('%s','localize','grps','-v7.3')",savefile_grps);
            fprintf(fid,"\nexit");
            fclose(fid);
        else
            disp('ERROR');
        end
    end
    fclose(opp_fid);
    fclose(run_fid);
 end
 
end

function tmp_cell = divideIntoGrps(localize_grp_sizes, localize_atoms, test_dup_flag, num_process)
	
%%
	consider_grps = [];
	count = 0;
	for kk = 1:length(localize_atoms)
		if test_dup_flag(kk)
			count = count +1;
			consider_grps(count) = kk;
		end
    end
	
    localize_grp_sizes_consider = localize_grp_sizes(consider_grps);
    
    [sorted_grpsize, sorted_grp_indx ] = sort(localize_grp_sizes_consider);
    consider_grps_sorted = consider_grps(sorted_grp_indx);
%%
    begin_p = 1;
    end_p   = length(consider_grps_sorted);

    tmp_cell = cell(1,num_process);
    toggle = 1;
    while begin_p <= end_p
        if toggle
            for i=1:num_process
                tmp_cell{i} = [tmp_cell{i}, consider_grps_sorted(begin_p)];
                begin_p = begin_p + 1;
                if begin_p > end_p
                    break; 
                end
            end
            toggle = 0;
        else
            for i=1:num_process
                tmp_cell{i} = [tmp_cell{i}, consider_grps_sorted(end_p)];
                end_p = end_p -1;
                if end_p < begin_p
                    break; 
                end
            end
            toggle = 1;
        end
    end
end


% function grps = divideIntoGrps(localize_grp_sizes, localize_atoms, test_dup_flag, num_processes)
% 	
% %%
% 	consider_grps = [];
% 	count = 0;
% 	for kk = 1:length(localize_atoms)
% 		if test_dup_flag(kk)
% 			count = count +1;
% 			consider_grps(count) = kk;
% 		end
% 	end
% 	
% %%
% 	
% 	size_per_grp = ceil(length(consider_grps)/num_processes) ;    
% 	grps = cell(1,num_processes);
%     count = 0;
% 	for i = 1:size_per_grp:length(consider_grps)
%         count = count +1;
%         if i+size_per_grp <= length(consider_grps)
%             grps{count} = consider_grps(i:i+size_per_grp-1);	
%         else
%             grps{count} = consider_grps(i:length(consider_grps));
%         end
% 	end
% end

