function reg_struct = reflectBylocalreflect(resi_start, resi_end, X,atom_map, Comp, pdbnameprefix)
%% call lcocalreflect func with x, y, z. choose 1 depending on ramachandran disallowed
%
%%
  tmp_pdb_dir = '../output_pdb';

  tmp_pdb_I11 = strcat(pdbnameprefix,'_I11');
  tmp_pdb_I22 = strcat(pdbnameprefix,'_I22');
  tmp_pdb_I33 = strcat(pdbnameprefix,'_I33');

  	reflect_try_I11 = localreflect(X, resi_start, resi_end, Comp, atom_map, tmp_pdb_I11, tmp_pdb_dir,1);
	pdb_file_name = strcat(tmp_pdb_dir, filesep, tmp_pdb_I11, '_reflect.pdb')
 	[rama11.out,fig_handler] = ramachandranMe(pdb_file_name,'Regions',true,'Glycine',true);
 	%resi in various area in the ramachandran map
 	[rama11.region,rama11.part_tot_resi] = getRamachandranReigionDistri(rama11.out.Angles);
 	%% get percent disallowed region for ramachandran map
 	rama11.disallowed = length(rama11.region(end).resi) / length([rama11.region.resi]);
	[indx_I11, region_name_I11] = findWhichRamaRegion([resi_start:resi_end],rama11.region)
    nt_allow_I11 = sum(indx_I11 >= 7); % partially allowed in 1st and 4th quadrant (allowed 4,5 , 6) and disallowed
        close all;

	reflect_try_I22 = localreflect(X, resi_start, resi_end, Comp, atom_map, tmp_pdb_I22, tmp_pdb_dir,2);
	pdb_file_name = strcat(tmp_pdb_dir, filesep, tmp_pdb_I22, '_reflect.pdb')
 	[rama22.out,fig_handler] = ramachandranMe(pdb_file_name,'Regions',true,'Glycine',true);
 	%resi in various area in the ramachandran map
 	[rama22.region,rama22.part_tot_resi] = getRamachandranReigionDistri(rama22.out.Angles);
 	%% get percent disallowed region for ramachandran map
 	rama22.disallowed = length(rama22.region(end).resi) / length([rama22.region.resi]);
	[indx_I22, region_name_I22] = findWhichRamaRegion([resi_start:resi_end],rama22.region)
    nt_allow_I22 = sum(indx_I22 >= 7); % partially allowed in 1st and 4th quadrant (allowed 4,5 , 6) and disallowed
        close all;

  	reflect_try_I33 = localreflect(X, resi_start, resi_end, Comp, atom_map, tmp_pdb_I33, tmp_pdb_dir,3);
	pdb_file_name = strcat(tmp_pdb_dir, filesep, tmp_pdb_I33, '_reflect.pdb')
 	[rama33.out,fig_handler] = ramachandranMe(pdb_file_name,'Regions',true,'Glycine',true);
 	%resi in various area in the ramachandran map
 	[rama33.region,rama33.part_tot_resi] = getRamachandranReigionDistri(rama33.out.Angles);
 	%% get percent disallowed region for ramachandran map
 	rama33.disallowed = length(rama33.region(end).resi) / length([rama33.region.resi]);
	[indx_I33, region_name_I33] = findWhichRamaRegion([resi_start:resi_end],rama33.region)
    nt_allow_I33 = sum(indx_I33 >= 7); % partially allowed in 1st and 4th quadrant (allowed 4,5 , 6) and disallowed
        close all;

    [val_nt_allow, pos_nt_allow] = min([nt_allow_I11, nt_allow_I22, nt_allow_I33]);
    [val_disallow, pos_disallow] = min([rama11.disallowed, rama22.disallowed, rama33.disallowed]);

    switch pos_nt_allow
           case 1
                if nt_allow_I11 == nt_allow_I33 
                       if rama11.disallowed < rama33.disallowed
                            reg_struct = reflect_try_I11;
                            disp('I11')
                       	else
                            reg_struct = reflect_try_I33;
                            disp('I33')
                       	end
                elseif nt_allow_I11 == nt_allow_I22
                       if rama11.disallowed < rama22.disallowed
                            reg_struct = reflect_try_I11;
                            disp('I11')
                       	else
                            reg_struct = reflect_try_I22;
                            disp('I22')
                       	end
                else
                        reg_struct = reflect_try_I11;
                        disp('I11')    
                end
           case 2
                if nt_allow_I22 == nt_allow_I33 
                       if rama22.disallowed < rama33.disallowed
                            reg_struct = reflect_try_I22;
                            disp('I22');
                       	else
                            reg_struct = reflect_try_I33;
                            disp('I33')
                       	end
                elseif nt_allow_I22 == nt_allow_I11
                       if rama22.disallowed < rama11.disallowed
                            reg_struct = reflect_try_I22;
                            disp('I22')
                       	else
                            reg_struct = reflect_try_I11;
                            disp('I11')
                       	end
                else
                        reg_struct = reflect_try_I22;
                        disp('I22')    
                end
           case 3
                if nt_allow_I33 == nt_allow_I11 
                       if rama33.disallowed < rama11.disallowed
                            reg_struct = reflect_try_I33
                            disp('I33')
                       	else
                            reg_struct = reflect_try_I11
                            disp('I11')
                       	end
                elseif nt_allow_I33 == nt_allow_I22
                       if rama33.disallowed < rama22.disallowed
                            reg_struct = reflect_try_I33;
                            disp('I33')
                       	else
                            reg_struct = reflect_try_I22;
                            disp('I22')
                       	end
                else
                        reg_struct = reflect_try_I33;
                        disp('I33')    
                end
    end

end
