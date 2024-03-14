function [cI_expanded,count] = expandCI_multiple(adj_mat,cI_resi,k,percent_cutoff,grp_size_cutoff,atleast_incr)
%% Call module expandCI multiple times until number of groups in which size change occurs
%  is less than "percent_cutoff" or the size of a group exceeds "grp_size_cutoff"
%  
%   Input:     adj_mat: resi x resi adjacency matrix each cell containing
%                       the number of upper and lower bounds between them.
%              cI_resi: test.cI_resi
%                    k: parameter to be supplied to expandCI module
%       percent_cutoff: lies between (0,1) - as explained in module defn
%      grp_size_cutoff: lies between (0,1) - as explained in module defn
%         atleast_incr: in each itter count only groups which increased
%                       size by atleast_incr
%
%  Output: cI_expanded: cI expanded.
%%
   if nargin~=6
       error('module: expandCI_multiple: incorrect number of arguments.');
   end

   if size(adj_mat,1) ~= size(adj_mat,2)
       error('module: expandCI_multiple: adjacency matrix incorrect.');
   end
   
   if ~isuppertriangular(adj_mat)
       if issymmetric(adj_mat)
           adj_mat = (triu(adj_mat,1))>0;
       else
           error('module: expandCI_multiple: adjacency matrix neither symmetric nor upper triangular.');
       end
   end
   
   if (percent_cutoff<=0) || (percent_cutoff>1)
       error('module: expandCI_multiple: param 0< percent_cutoff <1.');
   end
   
   if(grp_size_cutoff<=0)
       error('module: expandCI_multiple: grp_size_cutoff >0.');
   end
   
%%
   grp_incr = 1;   grp_size_min = 0; count=0;
   cI_expanded = cI_resi;
   num_cI      = length(cI_resi);
   
   while (grp_incr >= percent_cutoff && grp_size_min <= grp_size_cutoff)
      old_length  = cellfun(@length, cI_expanded);
      cI_expanded = expandCI(adj_mat,cI_expanded,k);
      
      curr_length  = cellfun(@length,cI_expanded);
      length_diff  = curr_length - old_length;
      grp_incr     = sum(length_diff > atleast_incr) / num_cI;
      grp_size_min = min(curr_length);
      count = count+1;      
   end
end   
    
   