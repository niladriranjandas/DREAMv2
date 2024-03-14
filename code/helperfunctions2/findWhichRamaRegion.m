function [indxarr, region_name] = findWhichRamaRegion(resilist, rama_region)
%%
%  rama_region is a  10*1 structure
%   rama_region(1).Name = 'Core Beta'  rama_region(1).resi = resi_list 1*n double
%   etc as per ramachandranregion in bioinfo toolbox or 
%%

	indxarr = zeros(1,length(resilist));
	region_name = strings(1,length(resilist));

		for j=1:length(rama_region)
		 	pos = ismember(resilist, rama_region(j).resi);
			indx = find(pos);

			indxarr(indx) = j;
			region_name(indx) = rama_region(j).Name;
		end

end
