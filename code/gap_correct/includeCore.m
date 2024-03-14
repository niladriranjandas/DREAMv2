function includeCore(multifile, gapsfile, seqfile, corepdb, newmultifile)
%%

%%
        fseq = fopen(seqfile);
        lines = textscan(fseq,'%s %s');
        fclose(fseq);
        col2 = str2double(lines{2});
        
        fgaps = fopen(gapsfile);
        strgaps = textscan(fseq,'%s');
        fclose(fseq);
        tmp1 = strgaps{1}
        tmp2 = split(tmp1,',');
        gaps = str2double(tmp2);
        
        cores = setdiff(col2, gaps);
        fnewmulti = fopen(newmultifile,'w')
        fprintf(fnewmulti,'%s ',corepdb);
        for i=1:length(cores)
        	if i==1
        		fprintf(fnewmulti,"%d",cores(i));
        	else
        		fprintf(fnewmulti,",%d",cores(i));
        	end
        end
	fprintf(fnewmulti," %d",cores(1));
        
        fid = fopen(multifile);
        if fid>0
                c=0
                tline = fgetl(fid);
                while ischar(tline)
		    if ~isempty(tline)
			    fprintf(fnewmulti,"\n%s",tline);
		    end
                    tline = fgetl(fid);
                end
        end
  fclose(fid);
  fclose(fnewmulti);
end
