function [density_i, density_j, density_ij, density_ij_diff] = getSubgraphDensity_v2(upl_i, upl_j, densities, resi_densities)
%%

%%
	len_resi_densities = length(resi_densities);
	len_densities = length(densities);

	if len_resi_densities ~= len_densities
		error('Module: getSubgraphDensity: lengths of resi_densities and densities missmatch');
	end

%%

	density_i       = 0;
	density_j       = 0;
	density_ij      = 0;
	density_ij_diff = 0;  

       flag_i       = 0;
       flag_j       = 0;
       flag_ij      = 0;
       flag_ij_diff = 0;
	for i=len_resi_densities:-1:1
		curr_resi_cells = resi_densities{i};
		for j = 1:length(curr_resi_cells)
			if flag_i == 0
				if ismember(upl_i, curr_resi_cells{j})
					density_i = densities(i);
					flag_i = 1;
				end
			end
			if flag_j == 0
				if ismember(upl_j, curr_resi_cells{j})
					density_j = densities(i);
					flag_j = 1;
				end
			end
			if flag_ij == 0
				if ismember(upl_i, curr_resi_cells{j}) && ismember(upl_j, curr_resi_cells{j})
					density_ij = densities(i);
					flag_ij = 1;
				end
			end
			if flag_i==1 && flag_j==1 && flag_ij==1
				break;
			end
		end
	end

	for i=len_resi_densities-1:-1:1
	    curr_resi_cells = resi_densities{i};	    
	    all_resis = [curr_resi_cells{:}];

	    if flag_ij_diff == 0
               if ismember(upl_i, all_resis) && ismember(upl_j, all_resis)
                      density_ij_diff = densities(i);
                  flag_ij_diff = 1;
               end
	    end
	    if flag_ij_diff == 1
	    	break;
	    end
	end
