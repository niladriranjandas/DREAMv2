function gaps = giveGaps(pdbfile, stage1mat, gapopfile)
	
      %% ---------------------------------------------------
	DIFF_TOL = 2;
	load(stage1mat, 'test')
	resis_ = test.cI_resi_expand;
	resis  = unique(cat(2,resis_{:}));

   %% ------------------------------------------------
    outliers = findOutliers(pdbfile);
    % include before and after
    
    
	resis_secondary = findSecondary(pdbfile);
	
	org_gaps = setdiff(unique([resis,outliers]), resis_secondary);
    
    %%
    load(stage1mat, 'num');
    load(stage1mat,'org_num');
    count = 0;
    %% this is a quick fix %%
    if org_num(1) <= 0
	org_num_tmp = org_num(org_num>0);
        org_num = org_num_tmp;
    end
    %%
    if org_num(1) < num(1)
       for i=org_num(1):num(1)-1
            count = count+1;
            gaps(count) = i;
       end
    end
    %%
    count = count+1;
    gaps(count) = org_gaps(1);
	
    for i=2:length(org_gaps)
        diff = org_gaps(i) - org_gaps(i-1);
   	if (diff>1)
	     if (diff <= DIFF_TOL)
		fillup = org_gaps(i-1)+1:org_gaps(i)-1;
		gaps = [gaps,fillup];
             end
        end	
        gaps = [gaps,org_gaps(i)];
    end
    % -- end few resis -- %
    if org_num(end) > num(end)
       for i=num(end)+1:org_num(end)
          gaps = [gaps,i]; 
       end
    end
    

    fgap = fopen(gapopfile, "w");
    if fgap > 0
	    for i=1:length(gaps)
		if i ~= 1
			fprintf(fgap, ",%d",gaps(i));
		else
			fprintf(fgap, "%d",gaps(i));
		end	
	    end
    end
    fclose(fgap);
end


function outliers = findOutliers(pdbfile)
%%
%
%%

	if ~isfile(pdbfile)
		error('%s File do not exist',pdbfile);
	end
	
	 %[rama.out,fig_handler] = ramachandranMe(pdbfile,'Regions',true,'Glycine',true);
	 [rama.out,fig_handler] = ramachandranMe(pdbfile,'Regions',true,'Glycine',true,  'Plot','None');
         %resi in various area in the ramachandran map
         [rama.region,rama.part_tot_resi] = getRamachandranReigionDistri(rama.out.Angles);

         outliers = rama.region(end).resi;
         %outliers_ = rama.region(end).resi;
	 %outlierstmp = [outliers_(1)];
         %for i=2:length(outliers_)-1
	 %	outlierstmp = [outlierstmp, outliers_(i)-1, outliers_(i), outliers_(i+1)];
	 %end
	 %outlierstmp = [outlierstmp, outliers_(end)];
	 %outliers = unique(outlierstmp);
end

function resis = findSecondary(pdbfile)
	resis = [];
	%cmd = strcat('./resisSecondary.sh',{' '},pdbfile);
	cmd = strcat('gap_correct/resisSecondary.sh',{' '},pdbfile);
	[shellstatus, shellout] = system(cmd{1});

	if shellstatus == 0
		c=0;
		resi_arrs = strsplit(shellout,',');
		for i=1:length(resi_arrs)
			if ~isempty(resi_arrs{i})
				c = c+1;
				resis(c) = str2num(resi_arrs{i});
			end
		end
	end	
end
