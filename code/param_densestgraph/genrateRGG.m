function [ adj_mat,time_mod ] = genrateRGG(n,r)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%GENRATERGG : Generate random geometric graph
%   
%     INPUT : no. of points, sensing radius
%    OUTPUT : adjacency matrix
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
tic;
if nargin ~= 2
    error('Module: genrateRGG: 2 arguments needed');
end

%% generate n random points inside a unit square
      data = GenData3(n);
      
      %%using pdist
     dist = squareform(pdist(data));
      
      %%using compute distance module instead of pdist.
%       [d_row,d_col,dist_d] = computeDistance(data');
%       dist = zeros(n);
%       dist(sub2ind(size(dist),d_row,d_col)) = dist_d;
      
   adj_mat = dist < r;
   time_mod=toc;
end

function data = GenData3(no_of_points)
%generate random points inisde a unit square


if nargin <= 0 
    error('incorrect no. of inputs');
end

x_strt = randn;
y_strt = randn;
x = x_strt + 1*rand(no_of_points,1);
y = y_strt + 1*rand(no_of_points,1);

data = horzcat(x,y);

end