function X_noref = doSPROSatomwise_v4(atom_grp,...
                                      rand_X_grp,Comp_Cq_grp,...
                                      eq_cons_all,up_bounds,wh_vdw_bounds,vdw_bounds,sdp_lo_bounds,lo_bounds,wh_lo_bounds,... 
                                      hydrogen_omission)
 
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
[~,size_atoms] = size(rand_X_grp);
if size_atoms > 200 %150
	[rawX, eV, eS, Vec_cons_1st_eq] = solve_sdpt3(U,eq_cons_all,up_bounds,sdp_lo_bounds,[],f);
else
	[rawX, eigens_cvx_spros, slacks_cvx_spros] = solve_by_cvx_SPROS(U, eq_cons_all, up_bounds, sdp_lo_bounds, [], f);
end

orig_rawX = rawX;
X_noref = rawX(1:3,:);

end
