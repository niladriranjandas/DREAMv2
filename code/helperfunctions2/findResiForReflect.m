function collect_resis_sort = findResiForReflect(pdb_file)
%%
% regions Allowed 4,5 and 6 as per ramachandranregions.m in bioinfo toolbox
%%
	[rama.out,fig_handler] = ramachandranMe(pdb_file,'Regions',true,'Glycine',true);
	[rama.region,rama.part_tot_resi] = getRamachandranReigionDistri(rama.out.Angles);

        collect_resis=[];
	for i=1:length(rama.region)
		if strcmpi(rama.region(i).Name , 'allowed 4') ...
		   || strcmpi(rama.region(i).Name , 'allowed 5') ...
                   || strcmpi(rama.region(i).Name , 'allowed 6') ...
		   || strcmpi(rama.region(i).Name , 'core L-Alpha')
			collect_resis = [collect_resis, rama.region(i).resi];
		end
		
		if strcmpi(rama.region(i).Name , 'disallowed')
			for j=1:length(rama.region(i).resi)
				resi_loc = find(rama.out.ResidueNum == rama.region(i).resi(j));
				torsion_angle = rama.out.Angles(resi_loc,:);
				phi_resi = torsion_angle(1);
				psi_resi = torsion_angle(2);
				
				if phi_resi >=0 && phi_resi <= 180
					collect_resis = [collect_resis, rama.region(i).resi(j)];
				end
			end
		end
	end

	collect_resis_sort = sort(collect_resis);

	zflag=1;
	for i=1:length(collect_resis_sort)-1
		if(collect_resis_sort(i) +1 == collect_resis_sort(i+1))
			if zflag==1
				fprintf('\n%d',collect_resis_sort(i));
				zflag=0;
			end
			fprintf(' %d',collect_resis_sort(i+1));
		else
			zflag = 1;
		end
	end

end
