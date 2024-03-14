function runSecondStage2(matfiles, opfile)
%%
%
%%
    fid = fopen(opfile,'w');
    if (fid > 0)
        count = 0;
        files = split(matfiles, ',');    
        for j = 1:length(files)
            load(strcat('../protein/',files{j}));%, structure_done, protein_name);
            if exist('structure_done','var') == 1
               if structure_done == 1
                    fprintf(fid,'%s\n',protein_name);
                    clearvars -except fid count files j
               end
            end
        end
        fclose(fid);
    else
        fprintf('\n%s:Error opening file\n',opfile)
    end 