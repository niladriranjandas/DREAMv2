function pdbnames = findWhichToChoose(pdbfiles,opfile)
%%
%   Off the files from multialign find which ones to choose
%   first sort by increasing order of outliers then break any
%   ties by seeing one having greater number of seconday structure
%
%   Output	pdbnames: {'name1','name2',...}
%    Input	pdbfiles: {'name1','name2',...}
%%
	numfiles = length(pdbfiles);
	pdbnames = cell(1,numfiles);
	alloutliers = zeros(1,numfiles);
	resissecond = zeros(1,numfiles);
	
	for i=1:numfiles
		alloutliers(i) = length(findOutliers(pdbfiles{i}));
		resissecond(i) = length(findSecondary(pdbfiles{i}));
    end
	
    new_mat = [alloutliers',resissecond'];
    [sorted_new_mat, index_new_mat] = sortrows(new_mat,[1,-2]);

	fid = fopen(opfile,'w');
	if fid>0
		for i=1:numfiles
			fprintf(fid,'%s\n',pdbfiles{i});
            pdbnames{i} = pdbfiles{i};
		end
	else
		error('ERROR: Module findWhichToChoose:could not open file for op');
	end

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

         outliers = rama.region(end).resi;
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
