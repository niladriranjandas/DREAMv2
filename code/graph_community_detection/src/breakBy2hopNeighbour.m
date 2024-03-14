function cI = breakBy2hopNeighbour(adj_mat)
%% function to search the members in the 2 hop neighbour of each node for dense graph
%  Input  adj_mat: adjacency matrix
%  output      cI: cells of dense groups
%%
   if (size(adj_mat,1) ~= size(adj_mat,2))
       error('Module breakBy2hopNeighbour: adjacency matrix incorrect.');
   elseif ~issymmetric(adj_mat)
     if ~isuppertriangular(adj_mat)
       error('Module breakBy2hopNeighbour: adjacency matrix not symmetric.');
     end
   end
   
   %check if adj_mat 0 binary
   validateattributes(adj_mat,{string(class(adj_mat))},{'binary'},'breakBy2hopNeighbour');

%%
param_start = 0.7;%0.2;
param_end   = 1;%0.9;%0.8;
%%
  if ~isuppertriangular(adj_mat)
      adj_mat = adj_mat+adj_mat';
  end
  
  adj_mat_hop2 = adj_mat^2;  %gotto implemenet fast method for this since adj_mat is sparse
  
  adj_mat_2_n_1 = bsxfun(@or,adj_mat,adj_mat_hop2>0);
  adj_mat_2_n_1(logical(eye(length(adj_mat_2_n_1)))) = 0; %making the leading diag zero;
  
 %% iterate through each node & find dense graph in its 2hop neighbourhood
    count=0;
    for i= 1:length(adj_mat_2_n_1)
       [~,neigh_nodes,~] = find(adj_mat_2_n_1(i,:));
       neigh_nodes = sort(neigh_nodes,i);
       adj_mat_tmp = adj_mat(neigh_nodes,neigh_nodes);
       
       [row,col] = find(adj_mat_tmp);
       if (length(row)>2) && (length(col)>2)
          for j=param_start*10:param_end*10       % param_start:0.1:param_end gives error in find            
              tmp_cI = mexdse(row,col,length(adj_mat),j/10)
       
         %%get the required group size for param value 'i'
           %[tmp_row,~,~] = find(param_val(:,2)==i);
          req_size = 5;%param_val(tmp_row(1),1); %------------------chenged----         
           for k=1:length(tmp_cI)
              if length(tmp_cI{k}) >= req_size
                  fprintf('\n --------(param_val: %f, grp-size: %d',i,length(tmp_cI{k}));
                  count = count+1;
                  cI{count} = sort(neigh_nodes(tmp_cI{k}));
              end
           end          
          end   
       end
    end
    
end    