function [Aeq, Aup, Alo] = getBoundsForAnchor_v2(vars, anchors, eq, up, lo)
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
   if size(anchors,1) > size(anchors,2)
        anchors = anchors';
   end
   if size(vars,1) > size(vars,2)
        vars = vars';
   end

   %---------- equality bound ------------------------    
   [Aeq.A, Aeq.X] = mexExtractAnchorsVars(eq(:,1:3),anchors,vars);
   
   %---------- upper bound ------------------------    
   [Aup.A, Aup.X] = mexExtractAnchorsVars(up(:,1:3),anchors,vars);

   %---------- lower bound ------------------------    
   [Alo.A, Alo.X] = mexExtractAnchorsVars(lo(:,1:3),anchors,vars);   
end