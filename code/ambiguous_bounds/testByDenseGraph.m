function [cI_test, ret_densities] = testByDenseGraph(adj_mat,eta_lo)

%% use densest subgraph to break the graph
%  Paper: Dense Subgraph Extraction with Application to community detection
%  Link:  http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=5677532
%
%  This code iterates through param values lying within a range and gathers
%  the clustered groups accordingly.
%  The code requires the file param_dense_subgraph.mat file having group
%  size and param value pair, that has been determined apriori.
%%***********************************************************************

%%
   %%adjacency matrix has to be upper triangular
   
   if ~isuppertriangular(adj_mat)
       adj_mat = triu(adj_mat,1);
   end
   
   %if the graph is disconnected break algo doesn't work
   %(requires matgraph package
   g=graph; sparse(g);
   set_matrix(g,(adj_mat+adj_mat')>0);
    if ~isconnected(g)
       error('Module: breakByDenseGraph:Graph corr to adjacency matrix is not connected.'); 
    end
    
%% iterate through the param values and take the groups accordingly

param_start = eta_lo;
param_end   = 1;

param_req_size = 3; %4 for localization 3 for density testing for ambiguous bounds
count_test = 0;    
   
    %iniialize cI (otherwise err occurs in case no dense group is found
    cI_test = cell(1,size(param_start:0.01:param_end,2));

%% iterate through the param_start to param_end and gather the groups   
   [row,col] = find(adj_mat);
   for i=param_start:0.01:param_end
      tmp_cI = mexdse(row,col,length(adj_mat),i); 
      tmp_cI_test= [];
      count = 0;                 
      for j=1:length(tmp_cI)
           if length(tmp_cI{j}) >= param_req_size
               fprintf('\n --------(param_val: %f, grp-size: %d',i,length(tmp_cI{j}));
               count = count+1;
               tmp_cI_test{count} = sort(tmp_cI{j});
           end
      end
      count_test = count_test+1;
      cI_test{count_test} = tmp_cI_test;
   end
   
   ret_densities = param_start:0.01:param_end;
end
           