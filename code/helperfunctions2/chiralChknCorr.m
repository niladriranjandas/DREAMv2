function [x, x_chiral, ref_info] = chiralChknCorr(x_chk, atm_map, up, lo, eq, wh_Comp)
%% module for chirality check and correction
%  Input:  x_chk: 3 * n
%        atm_map: 1*n the atoms index corr to x_chk
%             up:upper bound    lo:lower bound     eq:equality bound
%        wh_Comp:
% Output:  x : 3 * n  after chiral correction followed by refinement
%    x_chiral: 3 * n  after chiral correction
%    ref_info: info struct returned by post
%%
   if nargin ~= 6
       error('Module: chiralChknCorr: incorrect no. of inputs.');
   end
   
   [dim, n_pts] = size(x_chk);   
   if dim ~= 3 
       error('Module: chiralChknCorr: incorrect no. of inputs.');
   end
   
   if length(atm_map) ~= n_pts
       error('Module: chiralChknCorr: atom map size doesn''t matches size of x_chk.');
   end
   
   if max(up(:,1)) > n_pts || max(lo(:,1)) > n_pts || max(eq(:,1)) > n_pts
        error('Module: chiralChknCorr: supplied bounds matrix should be reindexed corr to atm_map: (module: extractCons).');
   end
%%
   %create the wh_Comp for chirality_check
   tmp_wh_Comp.seq        = wh_Comp.seq;
   tmp_wh_Comp.num_seq    = wh_Comp.num_seq;
   tmp_wh_Comp.atom_names = wh_Comp.atom_names(atm_map);
   tmp_wh_Comp.residue    = wh_Comp.residue(atm_map);
   tmp_wh_Comp.atom_types = wh_Comp.atom_types(atm_map);
   tmp_wh_Comp.residue_bias = [0, find(tmp_wh_Comp.residue(2:end) - tmp_wh_Comp.residue(1:end-1))'];
   
   %chiral check
   chiral_rep = chirality_check(x_chk, tmp_wh_Comp, 1);%0);
   
   if any(~chiral_rep)
       x_chiral = chirality_correction(x_chk,tmp_wh_Comp,chiral_rep);              
       % call the refine module
%          [x,ref_info] = postProcessingMe(x_chiral,...
%                                          eq, lo, up,...
%                                          [1e4,1,1,-0.001],[10,10,10,10,10]);  
         [x,ref_info] = postProcessingMe(x_chiral,...
                                         eq, lo, full(up),...
                                         [1e4,1,1,-0.001],[10,10,10,10,10]);    
   else
      x = x_chk;
      x_chiral = x_chk;
      ref_info = nan;
   end
end