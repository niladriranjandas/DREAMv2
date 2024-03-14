function [x_cord,indices,rot,L_inv,B] = doGretSDP(cell_x,index_x,embed_dim,atm_map)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% do Global registration using GRET_SDP:
%  link: http://arxiv.org/pdf/1310.8135v3.pdf
%  Input:  cell_x : (1 x #grps) cells each cell having x_coordinates (3 * no. of points)
%          index_x: (1 x #grps) cells each cell having atom indices
%  Output: x_cord : global co-ordinates
%          indces : atom index of the globaly registered atoms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% check the inputs
   if nargin ~=4
       error('module: doGretSDP: Not enough inputs ');
   end
      
   num_patch = length(cell_x);
   if num_patch ~= length(index_x)
       error('module: doGretSDP: mismatch: no. of groups to be registed != no. of index cell supplied');
   end
   
   %no. of points in x_cord{i} should be equal to length(index_x{i})
   uniq_index = [];
   for i=1:num_patch
       [dim,num_points] = size(cell_x{i});
       if dim ~=embed_dim
           error('module: doGretSDP: cell_x{i} should be 3 * no. of points');       
       elseif num_points ~= length(index_x{i})
           fprintf('\n.... Group Index: %d \n',i);
           error('module: doGretSDP: no. of points of cell_x{i} and index_x{i} length mismatch');
       else
           uniq_index = union(uniq_index, index_x{i});
       end
   end
   num_pts = length(uniq_index);   
   
 %% form the BL matrix
   [B,L,Adj] = formBL(cell_x,index_x,num_pts,dim);
   
 %% check if Adj:adjacency matrix of membership graph is connected (req matgraph libr)
   g = graph;
   sparse(g);
   set_matrix(g,Adj);
   if ~isconnected(g)
      error('module: doGretSDP: membership graph is not connected.');
   end
 %% 
   B_Ldag_Bt = B * inverse(L) * B';
   
   if ~issymmetric(B_Ldag_Bt)
      B_Ldag_Bt = 0.5 .* (B_Ldag_Bt + B_Ldag_Bt'); 
   end
   
   Md = num_patch * dim;
   I  = eye(dim);
   G  = eye(Md);  % initialize with I_d
   
 %% call cvx module
    cvx_precision best;
    cvx_begin
     variable G(Md,Md) semidefinite     
     maximize (trace(B_Ldag_Bt*G))               %objective function
     subject to
              for i=1:num_patch
                  G(dim*(i-1)+1:dim*i,dim*(i-1)+1:dim*i) == I;
              end
    cvx_end
    
 %% check if G is symmetric   
  if ~issymmetric(G)
      fprintf('Module Gret SDP: G is not symmtric.');
      G = 1/2 .* (G+G');
  end
 
 %% deterministic rounding
  fprintf('\n---- Module: doGretSDP: rank(G): %d -----\n', rank(G,0.0001));
  [V,S] = eigs(G,dim,'LA');
  S=real(S); V=real(V);
 
  W = sqrt(S) * V';
  
 %% extract the rotation 
   ortho_rot = zeros(dim,Md);
   rot = cell(1,num_patch);
   for i=1:num_patch
      [U,~,V] = svd( W(:,dim*(i-1)+1:i*dim));
      ortho_rot(:,dim*(i-1)+1:i*dim) = U * V';
      rot{i} = ortho_rot(:,dim*(i-1)+1:i*dim);
   end
   
 %% extract the co-ordinates
 L_inv = pinv(L);
  x_cord = ortho_rot * B * L_inv;
  x_cord = x_cord(:, 1:num_pts);
  
  indices = uniq_index;
   
       
end