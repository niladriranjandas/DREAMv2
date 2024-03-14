function resis_outliers = getOutliers(pdbfile)
%% just get the outliers from matlab's bioinfo package 
	resis_outliers = findOutliers(pdbfile)
end

function outliers = findOutliers(pdbfile)
%%
%
%%

	if ~isfile(pdbfile)
		error('%s File do not exist',pdbfile);
	end
	
	 %[rama.out,fig_handler] = ramachandranMe(pdbfile,'Regions',true,'Glycine',true);
	 [rama.out,fig_handler] = ramachandranMe(pdbfile,'Regions',true,'Glycine',true,'Plot','None');
         %resi in various area in the ramachandran map
         [rama.region,rama.part_tot_resi] = getRamachandranReigionDistri(rama.out.Angles);

         %outliers = rama.region(end).resi;
         outliers_ = rama.region(end).resi;
         outlierstmp = [outliers_(1)];
         for i=2:length(outliers_)-1
               outlierstmp = [outlierstmp, outliers_(i)-1, outliers_(i), outliers_(i+1)];
         end
         outlierstmp = [outlierstmp, outliers_(end)];
         outliers = unique(outlierstmp);
end

