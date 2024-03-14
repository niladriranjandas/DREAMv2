function bgr_val = getRGBValsByDensity(adj_mat, div_cells)

if ~isuppertriangular(adj_mat)
    adj_mat = triu(adj_mat);
end
%% get density of each divided fragment
n_vert = length(adj_mat);
density_full = nnz(adj_mat)/(n_vert * (n_vert-1));
density = density_full * ones(1,n_vert);
bgr_val = density_full * ones(n_vert,3);

for i=1:length(div_cells)
   vertices  = div_cells{i};
   frag_adj  = adj_mat(vertices, vertices);
   frag_size = length(vertices);
   % frag_adj should be upper triangular
   tmp  = 2 * nnz(frag_adj) / (frag_size * (frag_size-1));
    for j=1:length(vertices)
      density(vertices(j)) = max(density(vertices(j)),tmp);
    end
end

f_r = 255/1; f_b = 255/min(density); f_g = 255/min(density);
for i=1:n_vert
    bgr_val(i,:) = [ 1, 1, round(density(i)*f_r)];
end

end
