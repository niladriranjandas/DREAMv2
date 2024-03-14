function [X, info] = anchoredLocalize(X, X_indx,A_indx, eq_bounds, lo_bounds, up_bounds, w, f)
%% function which calls the bfgs module to go the anchored localition
%   Input:      X: 3 x no. of atoms ( X0 and anchors together (X0 choosen by modeller) )
%                A: 3 x no. of anchors
%      eq_bounds.X: (i,j,eq_ij)
%      eq_bounds.A: (i,j,eq_ij)
%
%      lo_bounds.X: (i,j,lo_ij)
%      lo_bounds.A: (i,j,lo_ij)
%
%      up_bounds.X: (i,j,up_ij)
%      up_bounds.A: (i,j,up_ij)
%
%                w: 1 x 4 (weight for eq, lo and up bounds)
%                f: 1 x   (weight for diff kinds of bounds e.g. H, disulphide bond etc.)
%
%  Output:       X: 3 x no. of atoms (returned by bfgs)
%            info: struct (optimization params - returned by bfgs)
%%
   if nargin ~= 8
       error('Module:anchoredLocalize: incorrect no. of inputs.');
   end
   
   [dim, n_atoms] = size(X);
   if dim ~= 3
       error('Module:anchoredLocalize: incorrect dim for X (3 x no. of atoms expected).');
   end
   
   if length(X_indx) + length(A_indx) ~= n_atoms
       error('Module:anchoredLocalize: (no. anchors + no. dim for vars) ~= X.');
   end
   
   if size(eq_bounds.X,2) ~= 3
       error('Module:anchoredLocalize: incorrect size for eq_bounds (i,j,eq_ij) expected.');
   end
   if size(eq_bounds.A,2) ~= 3
       error('Module:anchoredLocalize: incorrect size for eq_bounds (i,j,eq_ij) expected.');
   end
   
   if size(lo_bounds.X,2) ~= 3
       error('Module:anchoredLocalize: incorrect size for lo_bounds (i,j,lo_ij) expected.');
   end
   if size(lo_bounds.A,2) ~= 3
       error('Module:anchoredLocalize: incorrect size for lo_bounds (i,j,lo_ij) expected.');
   end
   
   if size(up_bounds.X,2) ~= 3
       error('Module:anchoredLocalize: incorrect size for up_bounds (i,j,up_ij) expected.');
   end
   if size(up_bounds.A,2) ~= 3
       error('Module:anchoredLocalize: incorrect size for up_bounds (i,j,up_ij) expected.');
   end
   
   if length(w) ~= 7
       error('Module:anchoredLocalize: w must be 1 x 8.');
   end
%% set params for bfgs call
%  obj  = mexFindObjAnchored(X, X_indx, eq_bounds.X, eq_bounds.A, ...
%                                   up_bounds.X, up_bounds.A, ...
%                                   lo_bounds.X, lo_bounds.A,w);
%   obj  = mexFindObjAnchored_omp(X, X_indx, eq_bounds.X, eq_bounds.A, ...
%                                 up_bounds.X, up_bounds.A, ...
%                                 lo_bounds.X, lo_bounds.A,w);     
  obj  = mexFindObjAnchored_omp_v2(X, X_indx, eq_bounds.X, eq_bounds.A, ...
                                   up_bounds.X, up_bounds.A, ...
                                   lo_bounds.X, lo_bounds.A,w);
                              
  w(7) = w(7)*obj/(25*trace(X(:,X_indx)'*X(:,X_indx)));   % see this further.
  
  pars.nvar          = numel(X(:,X_indx));%length(X_indx);%numel(X0);
  pars.fgname        = 'objngradatx';
  pars.X             = X;
  pars.X_indx        = X_indx;
  pars.A_indx        = A_indx;
  pars.equality_cons = eq_bounds;
  pars.up_bounds     = up_bounds;
  pars.lo_bounds     = lo_bounds;
  pars.dim           = dim;
  pars.w             = w;
  %pars.f            = f;
  
  tmp              = X(:,X_indx);
  options.x0       = tmp(:);
  options.prtlevel = 0;
  options.maxit    = 50;%250;
  options.normtol  = 1e-15;
    
 %[x,f,d,H,     iter,info,     X,G,w,fevalrec,     xrec,     Hrec]
  [X,~,~,~,info.iter,info.stop,~,~,~,info.fevalrec,info.xrec,info.Hrec] = bfgs(pars,options);
  
end
  
  