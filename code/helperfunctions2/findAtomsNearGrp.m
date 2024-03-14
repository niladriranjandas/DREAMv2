function [Idx_mapped, atomIndices,D] = findAtomsNearGrp(source_pts, target_pts, dist_cutoff)
%%
% Given an array of atoms find atoms within a distance cutoff. If needed will modify with atom radii
%  Input:	source_pts.indx	:	atom indices
%			   source_pts.X :	3*|pts| 3D coordinates
%           target_pts.indx :   atom_indices
%              target_pts.X :   3*|pts| 3D coordinates
% Find points in target_pts.X such that they lie at a dist < dist_cutoff from any atom in source_pts.X
%%

	[r,c] = size(source_pts.X);
	if r ~= 3
		error('Module: findAtomsNearGrp: size error: source_pts.indx should be 3xn');
	elif c ~= length(source_pts.indx)
		error('Module: findAtomsNearGrp: size mismatch error: no. of indices in source_pts.pts and source_pts.indx should match');
	end

	[r1,c1] = size(target_pts.X);
	if r ~= 3
		error('Module: findAtomsNearGrp: size error: source_pts.indx should be 3xn');
	elif c ~= length(target_pts.indx)
		error('Module: findAtomsNearGrp: size mismatch error: no. of indices in source_pts.pts and source_pts.indx should match');
	end

	%%

	[Idx,D] = rangesearch(target_pts.X', source_pts.X', dist_cutoff);
	Idx_mapped = cell(1,length(Idx));

	for i=1:length(Idx)
		curr_Idx = Idx{i};
		mapped_idx = zeros(1,length(curr_Idx));
		for j=1:length(curr_Idx)
		    mapped_idx(j) = target_pts.indx(curr_Idx(j));
		end
		Idx_mapped{i} = mapped_idx;
	end

    
	atomIndices = unique([Idx_mapped{:}]);

end

function [Idx_mod, D_mod] = findDistWithVanDerWaalsRadii(Idx,D,atom_names, A_info)
 	
 	for i=1:length(Idx)
 		curr_Idx = Idx{i};
 		curr_D   = D{i};

 		
 	end


end
