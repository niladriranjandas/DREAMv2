function giveFirstResi(gapsfile, seqfile, firstresifile)
%%

%%
        fseq = fopen(seqfile);
        lines = textscan(fseq,'%s %s');
        fclose(fseq);
        col2 = str2double(lines{2});

        fgaps = fopen(gapsfile);
        strgaps = textscan(fgaps,'%s');
        fclose(fgaps);
        tmp1 = strgaps{1};
        tmp2 = split(tmp1,',');
        gaps = str2double(tmp2);

	fid = fopen(firstresifile, "w");
	if fid>0
		residiff = setdiff(col2,gaps);
		fprintf(fid,"%d",residiff(1));
	else
		error("MODULE: giveFirstResi: O/P file cannot be opened")
	fi
	fclose(fid);
end
