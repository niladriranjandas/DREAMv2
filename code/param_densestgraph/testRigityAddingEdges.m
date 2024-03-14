function eta = testGlobalRigityAddingEdges(n)
%%  find the number of edges in a graph required for density param to be eta
%   Input:  n  : number of vertices in the graph.
%   Output  eta: value of parameter for dense subgraph algorithm.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 

%% generate the rgg
 r = 0.8;
 dims = 3;
 runs = 3; %for gotler rigidity
 [ adj_mat,~] = genrateRGG(n,r);
 
 %% get the spanning tree
 if ~isuppertriangular(adj_mat)
     adj_mat = triu(adj_mat,1);
 end
 
 [T,~] = computespanningtree(adj_mat);

 %% add random edges to the spanning tree and test for rigidity
 is_rigid = 0;
 
 %% find the edges absent in tree T
 tmp = triu(ones(n),1);
 edge_absent = bsxfun(@and,tmp,~T);
 [row_,col_,~] = find(edge_absent);
 no_edge    = length(row_);
 
 %% add edge randomly and see it its globally_rigid
     while ~is_rigid && no_edge > 0
         edge_to_add = randi([1 no_edge],1);
         T(row_(edge_to_add),col_(edge_to_add)) = 1;
       %%removing the edge_to_add from list of indices having no edges
           row_ = row_([1:edge_to_add-1,edge_to_add+1:end]);
           col_ = col_([1:edge_to_add-1,edge_to_add+1:end]);
         no_edge = no_edge - 1;
         
         if ~isuppertriangular(T)
             T = triu(T,1);
         end
         %T = T + T';
         %%call the gotler rigidity module
         [ is_rigid, is_locally_rigid, GL_dof, WW ] = ...
                                   SVD_glob_rigidity_many_times_G(T+T',dims,runs);
                    
     end

     if is_rigid 
       eta = nnz(T) * 2 / (n * (n-1));
     else
       eta = 0;
     end
end