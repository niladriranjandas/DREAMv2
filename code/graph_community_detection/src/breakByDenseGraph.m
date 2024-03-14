function cI = breakByDenseGraph(adj_mat,eta_lo)

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
   
    %iniialize cI (otherwise err occurs in case no dense group is found
    cI = [];
%% iterate through the param values and take the groups accordingly

param_start = eta_lo;%0.8;%0.6;%0.65;%0.45;%0.55;%0.2;%0.5;
param_end   = 1;%0.9;%0.8;

param_req_size = 4;%5;

%% load the group-size and param-value pair from param_test_global table
%    if ~exist('param_dense_subgraph.mat','file')
%        error('Module breakByDenseGraph: param file not found');
%    end
%    
%    tmp_load  = load('param_dense_subgraph.mat','param_test_global');
%    %tmp_load  = load('param_dense_subgraph.mat','param_test_local');
%    
%    if isempty(fieldnames(tmp_load))
%        error('Module breakByDenseGraph: param value not found');
%    end
%    
%    param_test_global = tmp_load.param_test_global;
%    %param_test_global  = tmp_load.param_test_local;
%    
%    %%find the group-size param-value pair lying within the range and round
%    %%them to the nearest place after decimal
%    [row,~,~] = find(param_test_global(:,2)>=param_start & ...
%                        param_test_global(:,2)<param_end+0.1);
%    param_val = [param_test_global(row,1),floor(param_test_global(row,2)*10)/10];
   
%% iterate through the param_start to param_end and gather the groups   
   count = 0;
   [row,col] = find(adj_mat);
%   for i=param_start*10:param_end*10       % param_start:0.1:param_end gives error in find
%       i = i/10;
    for i=param_start:0.01:param_end
       tmp_cI = mexdse(row,col,length(adj_mat),i);
       
    %%get the requried group size for param value 'i'
    %[tmp_row,~,~] = find(param_val(:,2)==i);
    
    %param_req_size = 5;%param_val(tmp_row(1),1); %------------------chenged----
       for j=1:length(tmp_cI)
           if length(tmp_cI{j}) >= param_req_size
               fprintf('\n --------(param_val: %f, grp-size: %d',i,length(tmp_cI{j}));
               count = count+1;
               cI{count} = sort(tmp_cI{j});
           end
       end
   end
   
end
           