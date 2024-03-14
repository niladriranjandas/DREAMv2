
function [E, nr_edges, E_global ]= find_ALL_K2(Block, labels)
    [col_list , row_list] = find(Block); %%  switch the order to generate triangles in lexicographic order
    % take only entries above the diagonal to avoind double-counting
    ind_list = find(row_list <= col_list);
    row_list = row_list(ind_list);
    col_list = col_list(ind_list);
    % [row_list  col_list];
    
    nr_edges = length(row_list);
    E = [ row_list  col_list ];
    
    E_global =[];
    
    if nargin==2
        labels = reshape(labels, length(labels), 1);
        E_global = [labels(row_list)  labels(col_list)];
    end
    
    %
    % G_this_GRC = sprandsym(6, 0.4)
    % G_this_GRC = full(G_this_GRC ~=0);
    % for i=1:6; G_this_GRC(i,i)=0; end;
    % G_this_GRC
    % this_GRC_verts = [1 3 5 8 9 11]
    %
    %  [E_block_local, nr_edges, E_block_global]=find_ALL_K2(G_this_GRC, this_GRC_verts)
    %
 
end
