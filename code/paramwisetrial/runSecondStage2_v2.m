function runSecondStage2_v2(matfiles, opfile)
%%
%
%%
    
    xdata = [];
    count1 = 1;
    fid = fopen(opfile,'w');
    if (fid > 0)
        count = 0;
        files = split(matfiles, ',');    
        for j = 1:length(files)
            load(strcat('../protein/',files{j}));%, structure_done, protein_name);
            if exist('structure_done','var') == 1
               if structure_done == 1
                    paramsprot.name(count1,:) = protein_name;
                    paramsprot.cover(count1) =  double(length(reg.resi_index{1}))/double(max_res - min_res +1);                   
                    xdata = [xdata, paramsprot.cover(count1)];
                     count1 = count1+1;
                    %fprintf(fid,'%s\n',protein_name);
                    clearvars -except fid count files j paramsprot xdata count1
               end
            end
        end
        [m,s] = normfit(xdata);
        for k = 1:count1-1
        	if paramsprot.cover(k) >= m
        		fprintf(fid,'%s\n',paramsprot.name(k,:));
        	end
        end
        fclose(fid);
    else
        fprintf('\n%s:Error opening file\n',opfile)
    end
end
