function [B,L,Adj] = formBL(cell_x,index_x,n,dim)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function to form the B and L matrix for the GRET SDP
%  Input:  cell_x : cell of x_coordinates (3 * no. of points)
%          index_x: atom indices of the groups
%          n      : total no. of atoms
%          dim    : embedding dimenssion
%  Output: x_cord : global co-ordinates
%          indces : atom index of the globaly registered atoms
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

 %No input check is requried coz it's already in module doGretSDP
 
 %
   m  = length(cell_x);
   Md = m*dim;
   
 %% initialize the matrices
   %B   = zeros(Md,n+m);
   B = [];
   L   = zeros(n+m);
   Adj = zeros(n+m);
 %% create B and L
   count=0;
   for i=1:m       
       Btmp = zeros(dim,n+m);
       Btmp(:,index_x{i}) = cell_x{i};
       sumtmp = sum(cell_x{i},2);
       Btmp(:,n+i) = -1 .* sumtmp;
       B=[B;Btmp];
      % B(dim*count+1:dim*count+dim,index_x{i}) = cell_x{i};
       
       for j=1:length(index_x{i})
          k = index_x{i}(j);
          eij = zeros(n+m);
          eij(k,k)     = 1;
          eij(n+i,n+i) = 1;
          eij(k,n+i)   = -1;
          eij(n+i,k)   = -1;
          L = L + eij;          
       end
Adj(index_x{i},n+i) = 1;
Adj(n+i,index_x{i}) = 1;
   end
   
 %%degree matrix
 deg = diag(sum(Adj));
 L_adj = Adj - deg;
 
 norm(L - L_adj,'fro')
   
end
   
 
   