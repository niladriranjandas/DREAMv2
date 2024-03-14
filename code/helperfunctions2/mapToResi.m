function cI_resi_mapped = mapToResi(cI_resi,Comp_num_seq)
%% The residue sequence in the protein seq file (.seq ) may not begin from 1.
%  Thus the need for rectification. (bug_report_4_03_2016)
%   Input:        cI_resi: cells each block contains residue numbers
%            Comp_num_seq: Comp.num_seq
%  Output: cI_resi_mapped: mapped cI_resi
%%
   if nargin ~=2
       error('Module: mapToResi: not enough inputs.');
   end
   
%%
   grps           = length(cI_resi);
   cI_resi_mapped = cell(1,grps);
   for i=1:grps             
      grp_i = cI_resi{i};
      tmp   = zeros(size(grp_i));
       for j=1:length(grp_i)
          tmp(j) = Comp_num_seq(grp_i(j));
       end
     cI_resi_mapped(i) = {tmp};
   end 
           