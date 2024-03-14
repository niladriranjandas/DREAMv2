% pool = gcp;
% 
% for i=1:2  %length(test_cI_modified)
%     
%    f(i) = parfeval(pool, @doSPROSatomwise_v3, 1, localize_atoms{i},...
%                              rand_X_grp{i},Comp_Cq_grp{i},...
%                              localize_eq_cons_all{i},localize_up_bounds{i},localize_wh_vdw_bounds{i},localize_vdw_bounds{i},localize_sdp_lo_bounds{i},localize_lo_bounds{i},localize_wh_lo_bounds{i},... 
%                              hydrogen_omission);
% end
% 
% disp('here')

%%

pool = gcp;

[row, col] = find(test_dup_flag);
for i=1:length(col)
    
   f(i) = parfeval(pool, @doSPROSatomwise_v3, 1, localize_atoms{col(i)},...
                             rand_X_grp{col(i)},Comp_Cq_grp{col(i)},...
                             localize_eq_cons_all{col(i)},localize_up_bounds{col(i)},localize_wh_vdw_bounds{col(i)},localize_vdw_bounds{col(i)},localize_sdp_lo_bounds{col(i)},localize_lo_bounds{col(i)},localize_wh_lo_bounds{col(i)},... 
                             hydrogen_omission);
end

disp('here')

