
function [ is_globally_rigid, is_locally_rigid, GL_dof ] = ...
                        matlab_SVD_rank_global_rigidity_dofs(n,d,E)

                    
%%{
% disp(['n= ' int2str(n)]);
Z = rand(n,d);
m = length(E);
GL_dof = 999;

if m <= 2; is_globally_rigid=0; is_locally_rigid=0; GL_dof=999; return; end

if m==3; 
    verts = unique( [ E(:,1)'  E(:,2)' ]);
    if length(verts) ~= 3
        is_globally_rigid=0; is_locally_rigid=0; GL_dof=999; 
    else  % it's a triangle
         is_globally_rigid=1; is_locally_rigid=1; GL_dof=0; 
    end
return;
end

if nchoosek(n,2) ==m
   % disp('Complete graph..it is globally rigid'); 
   is_globally_rigid=1; is_locally_rigid=1; GL_dof=0; 
   return;
end


%% Build the rigidity matrix
R = zeros(m , n*d);
for j=1:m
    u = E(j,1);
    v = E(j,2);
    for t=1:d
        R(j, (u-1)*d + t ) = Z(u,t) - Z(v,t);
        R(j, (v-1)*d + t ) = Z(v,t) - Z(u,t);
    end
end

% is_red_rig = is_redundantly_rigid (n,R)

%% Fix bug! make sure the number of edges is > n*d - d(d+1)/2
% if m < n*d - d*(d+1)/2
%      is_globally_rigid=0; is_locally_rigid=0; GL_dof=999; 
%      return;
% end


%% Check for local rigidity
rk = rank(R, 10^(-5));
S  = svd(R);
% disp('===================');
%% NOV 3: 
rigid_rank = n*d - d*(d+1)/2; %% if the graph were rigid, this would be the rank
%rigid_rank = min(n*d,m) - d*(d+1)/2; %% if the graph were rigid, this would be the rank
dof_rank = rigid_rank - rk;
if dof_rank == 0
   % disp('Locally rigid');
   is_locally_rigid=1;
else
   % disp('NOT Locally rigid');
   is_locally_rigid=0;
   is_globally_rigid=0;
   return;
end



%%
%{
% velocities are the trivial d*(d-1)/2 vectors in the null space
% corresponding to orthogonal transformations
%         +  the d translations
velocities = zeros(d*n, d*(d-1)/2);
index = 0;
for coord1=1:d;
    for coord2=(coord1+1):d;
        index = index+1;
        velocities([coord1:d:d*N], index) = R_points(coord2,:)'; 
        velocities([coord2:d:d*N], index) = -R_points(coord1,:)';
    end;
end;
%}

%% Check for global rigidity
% disp('Checking for global rigidity...');

% first check if the number of edges is less than t
t = n*d - nchoosek(d+1,2);
if m < t
    disp('The graph has too few edges to be globally rigid...');
    is_globally_rigid = 0;     %    is_locally_rigid=0;
    GL_dof = t-m;
    return;
end



%% Step 1: find a random stress vector w_ij in null(R^T)

if rk ==m
    % disp... Left null space had dimension 0
    is_globally_rigid=0;
    GL_dof=999;
   return;
end

% Use lsqr to solve R'*x=0 with a random initialization x_0
 maxit = 500000;
% disp('Step 1: finding a random stress vector w_ij');
%}

% save freeezee

% load freeezee

% %  [w,flag,relres,iter] = lsqr(R', zeros(d*n, 1), 1e-15,maxit , [], [], rand(m, 1));
% %  w = w / norm(w)
% %  lsqr_norm = norm(R'*w)
% %  flag
% %  norm(w)
% %  size(w)
% % R

%%{
[U,S,V] = svd(R');
size(R);
S;
sing = diag(S);
poszeros = find(abs(sing) < 0.000000001);
nr_poszeros = length(poszeros);
%if nr_poszeros ~= d * (d+1) /2
    % save  graph_svd_prob n d E
    % input('Something is wrong...the number of nnzero singular values is not d(d+1)/2 !');
%end
randcoef = rand(1,length(poszeros));
randcoef_mtx = repmat(randcoef,length(V),1);
vect = V(:,poszeros);
size_V = size(V);
lincomb = vect .* randcoef_mtx;
lincomb  = sum( lincomb , 2);
lincomb = lincomb / norm(lincomb);
w = lincomb;
norm(w);
svd_w_norm = norm(R'*w);


% if (lsqr_norm < 1e-15)
 if (svd_w_norm < 1e-15)
       % disp(['lsqr converged after ' num2str(iter) ' iterations, residual ' num2str(lsqr_norm)]);
else
    if (flag == 1)
       % disp(['lsqr did not converge after ' num2str(maxit) ' iterations, residual ' num2str(lsqr_norm)]);
    elseif (flag == 0)
       % disp(['lsqr converged after ' num2str(iter) ' iterations, BIG RESIDUAL ' num2str(lsqr_norm)]);
    else
       % disp(['lsqr STAGNATED after ' num2str(iter) ' iterations, residual ' num2str(lsqr_norm)]);
    end;
end;
%}
% disp(['Found a random vector w_ij with norm(R^T*w) = ' num2str(norm(R'*w))]);

% if lsqr_norm > 10^(-6)
if svd_w_norm > 10^(-6)
    disp('Large svd_w_norm ...NOT globally rigid');
    is_globally_rigid = 0;
    return;
end

%% Step 2: Form the stress matrix W and calculate its rank
% disp(['Step 2: Form the stress matrix W and calculate its rank']);
W = sparse( E(:,1), E(:,2), w, n, n);
%full(W);
W = W + W'; % symmetrize W by adding its transpose
%full(W);

% make the matrix row stochastic
row_sums = sum(W);
diago = (1:n);
lin_ind_diags = sub2ind(size(W),diago, diago);
W(lin_ind_diags) = - row_sums;
%full(W);

% full(W);
s = n-d-1;
W_rank = rank(full(W), 10^(-10));


if W_rank==s
    is_globally_rigid =1;
    GL_dof=0;
    % disp('Globally rigid');
else
    is_globally_rigid =0;
    GL_dof = s - W_rank;
    %disp(['MATLAB Rank of W: ' num2str(W_rank)]);
    %sing_values = svd(full(W))'
    
    % disp('NOT Globally rigid');
    % input(' press key...');
end

%{
%% Check nullspace of stress matrix using lsqr
% Clearly, W*Z' = 0
W*[Z ones(n,1)];
% Determine if the null space of W is larger by the usual randomized scheme
% Pick a random vector b
%disp(' ................... ');
tic;
b = randn(n,1);

% Project b on the orthogonal subspace to R_points\
[Q_points,R] = qr([Z ones(n,1)],0);
tmp = Q_points'*b;
tmp = Q_points *tmp;
b = b - tmp;
b = b / norm(b);

% Try to solve W*x=b
tol = 1e-5/sqrt(n)
% maxit = 100000;
maxit = 10000000;
[x,flag,relres,iter,resvec] = lsqr(W,b,tol,maxit);
t = toc;
%%disp(['  ' num2str(t) ' - elapsed time']);
if (relres < tol)
    disp(['lsqr converged after ' num2str(iter) ' iterations, residual ' num2str(relres)]);
    is_globally_rigid = 1;
    disp(['Globally rigid']);
else
    if (flag == 1)
        disp(['lsqr did not converge after ' num2str(maxit) ' iterations, residual ' num2str(relres)]);
    elseif (flag == 0)
        disp(['lsqr converged after ' num2str(iter) ' iterations, BIG RESIDUAL ' num2str(relres)]);
    else
        disp(['lsqr STAGNATED after ' num2str(iter) ' iterations, residual ' num2str(relres)]);
    end;
    is_globally_rigid = 0;
    disp(['NOT Globally rigid']);
end;

if is_globally_rigid_rank == is_globally_rigid
    %%disp('The two methods agree');
else
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('The two methods for global rigidity DO NOT agree');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp('!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!');
    disp(['MATLAB Rank of W: ' num2str(W_rank)]);
    format long
    %sing_values
    format short
    input('different...'); 
end
%}

end
