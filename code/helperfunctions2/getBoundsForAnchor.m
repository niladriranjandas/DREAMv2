function [Aeq, Aup, Alo] = getBoundsForAnchor(vars, anchors, eq, up, lo)
%% Function to prepare arrays (eq/upper/lower bounds) for anchor localization 
%    Input:   vars: (1 x no.) atom numbers specifying the variables
%          anchors: (1 x no.) atoms numbers specifying the anchors
%               eq: (i,j,eij) equality bounds (H-included)
%               up: (i,j,eij, type) upper bounds (H-included)
%               lo: (i,j,eij, type) lower bounds (H-included)
%   Output:    Aeq: Aeq.ex  (x_i, x_j, e_ij)
%                   Aeq.ea  (x_i, a_j, e_ij)
%              Aup: Aup.ex  (x_i, x_j, up_ij, type)
%                   Aup.ea  (x_i, a_j, up_ij, type)  
%              Alo: Alo.ex  (x_i, x_j, lo_ij, type)
%                   Alo.ea  (x_i, a_j, lo_ij, type)
%%
     if nargin ~= 5
        error('Module: getBoundsForAnchor: incorrect no. of inputs');
     end

     [~, col_eq] = size(eq);
     [~, col_up] = size(up);
     [~, col_lo] = size(lo);
     
     if col_eq ~= 3
         error('Module: getBoundsForAnchor: incorrect eq bound matrix dim.');
     end
     if col_up ~= 4
         error('Module: getBoundsForAnchor: incorrect upper bound matrix dim.');
     end
     if col_lo ~= 4
         error('Module: getBoundsForAnchor: incorrect eq bound matrix dim.');
     end
%%
   %---------- equality bound ------------------------     
   [ind_eq_i, ind_eq_j] = anchorInOneNotOther(eq(:,1:2),anchors);
   ind_eq_xij = varBothIniandj(eq,vars);
   
   Aeq.X = eq(ind_eq_xij,:);
   Aeq.A = eq(ind_eq_i,:);
   tmp      = [eq(ind_eq_j,2), eq(ind_eq_j,1), eq(ind_eq_j,3)];
   Aeq.A = [Aeq.A; tmp];
   
   %---------- upper bound ---------------------------
   [ind_up_i, ind_up_j] = anchorInOneNotOther(up(:,1:2),anchors);
   ind_up_xij = varBothIniandj(up,vars);
   
   Aup.X = up(ind_up_xij,:);
   Aup.A = up(ind_up_i,:);
   tmp      = [up(ind_up_j,2), up(ind_up_j,1), up(ind_up_j,3), up(ind_up_j,4)];
   Aup.A = [Aup.A; tmp];   
   
   %---------- lower bound ----------------------------
   [ind_lo_i, ind_lo_j] = anchorInOneNotOther(lo(:,1:2),anchors);
   ind_lo_xij = varBothIniandj(lo,vars);   
   
   Alo.X = lo(ind_lo_xij,:);
   Alo.A = lo(ind_lo_i,:);
   tmp      = [lo(ind_lo_j,2), lo(ind_lo_j,1), lo(ind_lo_j,3), lo(ind_lo_j,4)];
   Alo.A = [Alo.A; tmp];   
   
end

function [indx_xi_aj, indx_ai_xj] = anchorInOneNotOther(bound_indx_mat, anchors)
%% function to return pair (x_i, a_j) or (a_j, x_i) from bounds matrix
%
%%       
   r_i=[];  r_j=[];
   for i=1:length(anchors)
      [r_i_,~,~] = find(bound_indx_mat(:,1)==anchors(i));
      [r_j_,~,~] = find(bound_indx_mat(:,2)==anchors(i));
      r_i = union(r_i,r_i_);  r_j = union(r_j,r_j_);
   end
  
   indx_ai_xj = setdiff(r_i,r_j);
   indx_xi_aj = setdiff(r_j,r_i);

end

function indx = varBothIniandj(bound_indx_mat, vars)
   r_i=[];  r_j=[];
   for i=1:length(vars)
      [r_i_,~,~] = find(bound_indx_mat(:,1)==vars(i));
      [r_j_,~,~] = find(bound_indx_mat(:,2)==vars(i));
      r_i = union(r_i,r_i_);  r_j = union(r_j,r_j_);
   end
  
   indx = intersect(r_i,r_j);
end