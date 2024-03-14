%% check for valid structure based on ramachandran and set flag is_valid_structure

 %[rama.out,fig_handler] = ramachandranMe(pdb_file,'Regions',true,'Glycine',true);		 %19-09-20 no fig run
  [rama.out,~] = ramachandranMe(pdb_file,'Regions',true,'Glycine',true,'Plot','None'); %19-09-20 no fig run
  %resi in various area in the ramachandran map
 [rama.region,rama.part_tot_resi] = getRamachandranReigionDistri(rama.out.Angles);
        
 %% get percent disallowed region for ramachandran map
 
 rama.disallowed = length(rama.region(end).resi) / length([rama.region.resi]);
 
  if rama.disallowed < 0.4
      is_valid_structure = 1;
  else    
      is_valid_structure = 0;
  end
