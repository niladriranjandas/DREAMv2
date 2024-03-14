function [X, info] = postProcessingMe(X0,eq_bounds,lo_bounds,up_bounds,w,f)
%% This function performs post processing using the result from SDP localization
%  Same as hansp_post_processing in the SPROS paper's implementation except
%  equality bounds are not calculated again.
%%


%%

equality_cons = eq_bounds;%equality_con_former(ref_X,Comp,2);

pars.nvar = numel(X0);
pars.fgname = 'fgcalc';
% 
% f_hb  = f[0];
% f_ta  = f[1];
% f_vdw = f[2];
%f = [10 10 1];

%%
      lo_bounds = [lo_bounds, ones(size(lo_bounds,1),1)];
      up_bounds = [up_bounds, zeros(size(up_bounds,1),1)];

%%
ow = w;
w(4) = 0.0;

pars.equality_cons = equality_cons;
pars.up_bounds = up_bounds;
pars.lo_bounds = lo_bounds;
pars.dim = size(X0,1);

obj  = objfunmex(X0,equality_cons,lo_bounds,up_bounds,w,f);

w(4) = ow(4)*obj/(25*trace(X0'*X0));
pars.w = w;
pars.f = f;

options.x0 = X0(:);
options.prtlevel = 0;
options.maxit = 250;
options.normtol = 1e-15;%10^-9;
%options.quitLSfail = 0;
%[x f] = bfgs(pars,options);
[x, ~, ~, ~, ~, ~, ~, ~, ~, info.obj] = bfgs(pars,options);

max_i = numel(info.obj);
temp = nan(max_i,1);
for i = 1:max_i
   temp(i) = mean(info.obj{i}); 
end
info.obj = temp;

     
X = reshape(x,size(X0,1),size(X0,2));

end