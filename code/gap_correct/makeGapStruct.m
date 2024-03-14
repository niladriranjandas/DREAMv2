function gapstruct = makeGapStruct(gapfile)
%%

%%
	fid = fopen(gapfile);	
	if fid>0
		c=0
		tline = fgetl(fid);
		while ischar(tline)
		    c=c+1;
		    resis = split(tline, ",");
		    gapstruct(c).resi = str2num(resis{1}):str2num(resis{2});
		    tline = fgetl(fid);
		end
	end
	fclose(fid);
end
