function new_lo = getNewLoBounds(X, atm_map, Comp)
%%
%  go residue by residue and find new lo bounds
%  Input:
%           X : 3*n 
%     atm_map : 1*n atom index
%  Output:
%       new_lo: (i,j, dist) i and j are mapped as per atm_map
%%

    dist_cutoff = 2;
    penetration = 0.9; % as per vdw_maker.m
    resis = Comp.residue(atm_map);

    resi_i_old = resis(1);
    start_indx = 1;

    count=1;
	for i=2:length(resis)
		resi_i_new = resis(i);

		if resi_i_new ~= resi_i_old || i==length(resis)
          	     end_indx = i-1;

	             source.X   = X(:,start_indx:end_indx);
          	     source.indx= atm_map(start_indx:end_indx);
           	     target.X   = X(:,end_indx+1:end);
	             target.indx= atm_map(end_indx+1:end);

          	     [Idx_mapped, atomIndices,D] = findAtomsNearGrp(source, target, dist_cutoff);
            
	             for j = 1:length(source.indx)
          	 	   curr_j = Idx_mapped{j};
           		   for k = 1:length(curr_j)                            
                               if Comp.info(5,source.indx(j)) ~= 0 && Comp.info(5,curr_j(k)) ~=0  % excluding pseudo atoms as per vdw_maker.m
           	   		   vdw_i_j = Comp.info(2,source.indx(j)) + Comp.info(2,curr_j(k));
		                   new_lo(count,:) = [source.indx(j), curr_j(k), vdw_i_j*penetration, 1];
                		   count = count+1;
                               end
               		   end
           	     end
                     start_indx = i;
		end
		resi_i_old = resi_i_new;
	end

end
