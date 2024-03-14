function flag = checkGrpReflectedAfterReg(x_n_index, rama_angles, Comp)
%%
%  check if reg(i) were is disallowed region or not after cosolidate. if so do reflect
%      
%%
	tol_frac = 0.35;

	[rama.region,~] = getRamachandranReigionDistri(rama_angles);
	flag = zeros(1,length(x_n_index));
	for i =1:length(x_n_index)
		resis =  unique(Comp.residue(x_n_index(i).ind));   % check if its Comp.residue or something else
		 [indx, region_name] = findWhichRamaRegion(resis,rama.region);

		 disallowed = sum(indx==10);
		 coreLalpha = sum(indx==3);
		 allowed4   = sum(indx==7);
		 allowed5   = sum(indx==8);
		 allowed6   = sum(indx==9);

		 frac_      = sum(disallowed + coreLalpha + allowed4 + allowed5 + allowed6)/length(resis);

		 if frac_ >= tol_frac
		 	flag(i) = -1;
		 else
		 	flag(i) = 1;
		 end
	end

end
