function gaprankmatrix = getGapArray(pdbfile, gaps)
%%
%  Input:  pdbfile: pdbfile
%	      gaps: gaps[1].resi = array of gap-1
%                   gaps[2].resi = array of gap-2 etc
% Output:  gaprankmatrix: gaprankmatrix[1].flag = 1 for allowed 2 for partial or disallowed
%                          gaprankmatrix[1].dist = distance from furthest allowed region
%%

	%[rama.out,fig_handler] = ramachandranMe(pdbfile,'Regions',true,'Glycine',true);
	[rama.out,fig_handler] = ramachandranMe(pdbfile,'Regions',true,'Glycine',true, 'Plot','None'); 
	[rama.region,rama.part_tot_resi] = getRamachandranReigionDistri(rama.out.Angles);

	
	n_gaps = length(gaps);
	region_info = feval('ramachandranRegions');
	
	for i=1:n_gaps
	    curr_gap_flag = 1;
	    curr_gap_dist = -1;
		resi = gaps(i).resi;
        
        curr_dist_partial = [];
        curr_dist_disallow = [];
        
		for j = 1:length(resi)
	              curr_resi = resi(j);
	              curr_torsion = rama.out.Angles(curr_resi,:);
	              curr_region = getRamachandranReigion(curr_torsion);
		          phi = curr_torsion(1);
		          psi = curr_torsion(2);
                  
                  if (isnan(phi) || isnan(psi))
                     continue
                  end
		          
                  if (curr_region>3 && curr_region<10 && curr_gap_flag ~= -1)
                    curr_gap_flag = 0;
			        tmp = region_info(curr_region).Patch;
			        xv = tmp(1,:); yv = tmp(2,:);
                    [curr_d_min_tmp, x_d_min, y_d_min, is_vertex, idx_c, xc, yc, is_in_seg, Cer, Ppr] = p_poly_dist(phi, psi, xv, yv);
                    curr_dist_partial = [ curr_dist_partial, curr_d_min_tmp];
		          end
	              if curr_region == 10
	              	    curr_gap_flag = -1;
	              	
	              	    if(phi<0 & psi <0)
	              			     tmp = region_info(2).Patch;
	              	             xv = tmp(1,:); yv = tmp(2,:);
                  	             [curr_d_min, x_d_min, y_d_min, is_vertex, idx_c, xc, yc, is_in_seg, Cer, Ppr] = p_poly_dist(phi, psi, xv, yv);                  	              	
                        elseif (phi<0 & psi >0)
                  	    		 tmp = region_info(1).Patch;
	              	             xv = tmp(1,:); yv = tmp(2,:);
                  	             [curr_d_min1, x_d_min, y_d_min, is_vertex, idx_c, xc, yc, is_in_seg, Cer, Ppr] = p_poly_dist(phi, psi, xv, yv);
                  	    		 tmp = region_info(2).Patch;
	              	             xv = tmp(1,:); yv = tmp(2,:);
                  	             [curr_d_min2, x_d_min, y_d_min, is_vertex, idx_c, xc, yc, is_in_seg, Cer, Ppr] = p_poly_dist(phi, psi, xv, yv);
                  	             curr_d_min = min([curr_d_min1, curr_d_min2]);
                        elseif (phi >0 & psi < 0)
                  	    		tmp = region_info(2).Patch;
	              	             xv = tmp(1,:); yv = tmp(2,:);
                  	             [curr_d_min1, x_d_min, y_d_min, is_vertex, idx_c, xc, yc, is_in_seg, Cer, Ppr] = p_poly_dist(phi, psi, xv, yv);
                  	    		 tmp = region_info(3).Patch;
	              	             xv = tmp(1,:); yv = tmp(2,:);
                  	             [curr_d_min2, x_d_min, y_d_min, is_vertex, idx_c, xc, yc, is_in_seg, Cer, Ppr] = p_poly_dist(phi, psi, xv, yv);
                  	            curr_d_min = min([curr_d_min1, curr_d_min2]);
                        elseif (phi > 0 & psi > 0)
                  	    		tmp = region_info(1).Patch;
	              	            xv = tmp(1,:); yv = tmp(2,:);
                  	           	[curr_d_min1, x_d_min, y_d_min, is_vertex, idx_c, xc, yc, is_in_seg, Cer, Ppr] = p_poly_dist(phi, psi, xv, yv);
                  	    		tmp = region_info(2).Patch;
	              	            xv = tmp(1,:); yv = tmp(2,:);
                  	           	[curr_d_min2, x_d_min, y_d_min, is_vertex, idx_c, xc, yc, is_in_seg, Cer, Ppr] = p_poly_dist(phi, psi, xv, yv);
                  	    		tmp = region_info(3).Patch;
	              	            xv = tmp(1,:); yv = tmp(2,:);
                  	           	[curr_d_min3, x_d_min, y_d_min, is_vertex, idx_c, xc, yc, is_in_seg, Cer, Ppr] = p_poly_dist(phi, psi, xv, yv);
			                    curr_d_min = min([curr_d_min1, curr_d_min2, curr_d_min3]);
                        end 
                        curr_dist_disallow = [ curr_dist_disallow, curr_d_min];
                  end
                  
        end	
        gaprankmatrix(i).flag = curr_gap_flag;
        if curr_gap_flag == 0
            gaprankmatrix(i).dist = min(curr_dist_partial);  
        elseif curr_gap_flag == -1
            gaprankmatrix(i).dist = max(curr_dist_disallow);
        elseif curr_gap_flag == 1
            gaprankmatrix(i).dist = -1;  
        end
    end

end
