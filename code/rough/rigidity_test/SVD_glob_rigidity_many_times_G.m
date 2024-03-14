function [ is_globally_rigid, is_locally_rigid, GL_dof, WW ] = ... 
                           SVD_glob_rigidity_many_times_G(G,d,runs)

% Given a graph G, and dimension d, check if G is 
%             Globally Rigid in R^d

%   load crazy_one_graph_svdGR
%   load good_one_graph_svdGR
sum(G);
n = length(G);
E = find_ALL_K2(G);

for i=1:runs
    
    
[ is_globally_rigid, is_locally_rigid, GL_dof ] = ...
                matlab_SVD_rank_global_rigidity_dofs( n,d,E);
                % matlab_rank_global_rigidity_dofs( n,d,E);
    W(i,:) = [ is_globally_rigid, is_locally_rigid, GL_dof ];
end
WW=W;
W = mode(W,1);

is_globally_rigid = W(1);
is_locally_rigid = W(2);
GL_dof = W(3);


end
