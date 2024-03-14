function writeValidModelsOnly(paramname, matfile,filename)
%%

%%
    fid = fopen(filename ,'a+');
    if fid>0
       load(matfile, is_valid_structure);
       if is_valid_structure
           fprintf(fid, '%s\n', paramname);
       end
       fclose(fid);
    end    
end