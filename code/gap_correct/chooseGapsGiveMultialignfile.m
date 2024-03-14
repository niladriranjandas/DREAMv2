function chooseGapsGiveMultialignfile(pdbfile, seqfile, gaps, op_multialign, op_gaps)
%%
%  Input: pdbfile={'pdb1.pdb','pdb2.pdb',...,'pdb20.pdb'}
%             gaps=[i,j...]
% 
%%

        %% --- choose gap
	for i=1:length(pdbfile)
		gaprankmatrix = getGapArray(pdbfile{i},gaps);
		gaparr(i,:) = [gaprankmatrix.flag];
		distarr(i,:) = [gaprankmatrix.dist];
		close all;
	end
	
	[modelgaps, gap_fill_struct] = findWhoseGap(gaparr, distarr, gaps);
	
	%% --- write file for multi-align
	gaps = [];
	fgaps = fopen(op_gaps,'w');
	fseq = fopen(seqfile);
	lines = textscan(fseq,'%s %s');
	fclose(fseq);
	col2 = str2double(lines{2});
	fid = fopen(op_multialign,'w');
	for i=1:length(gap_fill_struct)		
		fprintf(fid,"\n%s",pdbfile{gap_fill_struct(i).models})
		resis=gap_fill_struct(i).resi;
		first_resi = resis;
		%first_resi = setdiff(col2,resis);
		for j=1:length(resis)
			if j==1
				fprintf(fid," %d",resis(j))			        	
			else
		        	fprintf(fid,",%d",resis(j))
			end	
			gaps = [gaps,resis(j)];		
		end
		fprintf(fid," %d",first_resi(1))
	end
	fclose(fid);

	% -- write the gaps file -- %
	gaps_sort = sort(gaps);
	for i=1:length(gaps_sort)
		if i==1
			fprintf(fgaps,"%d",gaps_sort(i));
		else
			fprintf(fgaps,",%d",gaps_sort(i));
		end
	end
	fclose(fgaps);
end
