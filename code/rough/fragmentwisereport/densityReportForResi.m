% mkdssp -i 2M4K_0001.pdb | awk '/#  RESIDUE/,/^$/' | awk '{print $1 " " substr($0,17,1)}' | awk '{gsub( "  "," -" ); printf "\n%d '%s'", $1 , $2}' | egrep 'E|B|S|T|H|I|G' | awk '{print $1}'| tr '\n' ','
%resi_list=[3,4,5,7,8,9,10,11,12,13,15,16,17,18,19,20,21,22,23,24,25,26,28,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,47,48,49,50,51,52,53,54,55,56,57,58,59,60,61,62,63,65,66,67,68,69,71,72,73,74,75,76,77,78,79,80]


function density_resi = densityReportForResi(wh_Comp,wh_up_bounds, resi_list)

    %natoms = length(wh_Comp.atom_names);
    natoms = 0;
    for i = 1:length(wh_Comp.residue)
       if( ismember(wh_Comp.residue(i), resi_list))
           natoms = natoms + 1;
       end
    end
    num_bounds = length(wh_up_bounds);
    count=0;
    for i = 1:num_bounds
        atm_i = wh_up_bounds(i,1);
        atm_j = wh_up_bounds(i,2);
        
        resi_i = wh_Comp.residue(atm_i);     
        resi_j = wh_Comp.residue(atm_j);
        
        if( ismember(resi_i, resi_list) && ismember(resi_j, resi_list))
           count = count+1; 
        end
    end
    
    density_resi = 2*count/(natoms * (natoms-1));
end

