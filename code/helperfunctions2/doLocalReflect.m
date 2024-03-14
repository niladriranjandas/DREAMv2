function reg_struct = doLocalReflect(pdbname, X,atm_map, Comp)
%% first call findResiForReflect then call reflectBylocalreflect
%
%

       %% determine seq of resis to reflect
	resilist = findResiForReflect(pdbname);
  close all;     
        begin_resi = []; end_resi = [];
        diff_resi = resilist(2:end) - resilist(1:end-1);
        if resilist(1)+1 == resilist(2)
           diff_resi_tmp = [1, diff_resi];
        else
           diff_resi_tmp = [0, diff_resi];
        end
        diff_resi_ = (diff_resi_tmp==1);	
		count1 = 1; count2 =1;
		for i=1:length(diff_resi_)-1
		    if i==1
		        if diff_resi_(i)==1 && diff_resi_(i+1) == 1
		           begin_resi(count1) = resilist(i);
		           count1 = count1+1;
		        end
		    end
		    if diff_resi_(i) == 0 && diff_resi_(i+1) == 1
		        begin_resi(count1) = resilist(i);
		        count1 = count1+1;
		    end
		    if diff_resi_(i) == 1 && diff_resi_(i+1) == 0
		        end_resi(count2) = resilist(i);
		        count2 = count2+1;
		    end   
		    if i == length(diff_resi_)-1
		         if diff_resi_(i+1) == 1
		            end_resi(count2) = resilist(i+1);
		         end
		     end                       
		end        
    
     %% do the reflection    
     if length(begin_resi) ~= length(end_resi)
        error('Module: doLocalReflect: begin_resi end_resi length mismatch');
     end

    reg_struct.reg_X_noref = X;
    reg_struct.reg_atom_map = atm_map; 
    for i =1:length(begin_resi)
         if end_resi(i) - begin_resi(i) > 1
             tmp1 = split(pdbname,filesep);
             tmp = split(tmp1{end},'.pdb');
             pdbprefix = strcat(tmp{1},'_',string(begin_resi(i)),'_',string(end_resi(i)));
             reg_struct = reflectBylocalreflect(begin_resi(i), end_resi(i), reg_struct.reg_X_noref, reg_struct.reg_atom_map, Comp, pdbprefix);
         end
     end
end
