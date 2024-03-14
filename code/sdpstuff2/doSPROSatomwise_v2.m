function X_noref = doSPROSatomwise_v2(atom_grp,...
                                      rand_X,Comp_Cq,...
                                      eq_cons_all,up_bounds,wh_vdw_bounds,vdw_bounds,sdp_lo_bounds,lo_bounds,wh_lo_bounds,... 
                                      hydrogen_omission)

rand_X_grp = rand_X(:,atom_grp);

[I,J,~] = find(Comp_Cq(atom_grp,:)); % [I,J] are sorted in ascending order J wise

cq_grp = unique(J);

J_tmp = zeros(length(J),1);

 for i=1:length(J)
    J_tmp(i) = find(cq_grp==J(i)); 
 end
 
 Comp_Cq_grp = sparse(I,J_tmp,1,length(atom_grp),length(cq_grp));
 
 %[U, cliq_dims] = reducer(X, Comp, type)
 
    switch hydrogen_omission
      case 0
        [U, cliq_dims_grp] = reducer_Me(rand_X_grp,Comp_Cq_grp);        
      case 1
        dont_compute_U = 1;
        [U, cliq_dims_grp] = reducer_Me(rand_X_grp,Comp_Cq_grp);
     %   [wh_U, wh_Comp.cliq_dims] =
     %   reducer(wh_rand_X,wh_Comp,dont_compute_U);
    end    

%% Solve the SDP
%==========================================================================
%==========================================================================
display('----solve SDP- SPROS---')
f = [10,10,10,10,10];
[rawX, eV, eS, Vec_cons_1st_eq] = solve_sdpt3(U,eq_cons_all,up_bounds,sdp_lo_bounds,[],f);

orig_rawX = rawX;
X_noref = rawX(1:3,:);

end