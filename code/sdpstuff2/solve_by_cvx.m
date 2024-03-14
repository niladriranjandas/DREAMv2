function [X, eigens, slacks] = solve_by_cvx(num_atoms,equality_cons,upper_bounds,lower_bounds,vdw_bounds,f)

hydrogen_factor = f(1);
torsionu_factor = f(2);
torsionl_factor = f(5);
vdw_factor      = f(4);

num_upper_bounds    = size(upper_bounds,1)
num_equality_cons   = size(equality_cons,1)
num_lower_bounds    = size(lower_bounds,1)
num_vdw_bounds      = size(vdw_bounds,1)

% num_up_slacks  = 2*num_upper_bounds;
% num_lo_slacks  = 2*num_lower_bounds;
% num_vdw_slacks = 2*num_vdw_bounds;

gamma  = 50*(num_atoms)/num_upper_bounds
lambda = 0*(num_atoms)/num_upper_bounds;
%num_slacks = num_up_slacks + num_lo_slacks + num_vdw_slacks;

%num_cons = num_equality_cons + num_upper_bounds + num_lower_bounds + num_vdw_bounds + 1;
num_cons = num_equality_cons + num_upper_bounds + num_lower_bounds + num_vdw_bounds;

Q = eye(num_atoms);
I = eye(num_atoms);
e = ones(num_atoms,1);

display('-----CVX called--------');
  cvx_precision best;
  cvx_begin sdp
     variable Q(num_atoms,num_atoms) semidefinite            %define variables
     variable epsilon(num_upper_bounds) nonnegative
     variable zeta_l(num_lower_bounds)  nonnegative
     variable zeta_vdw(num_vdw_bounds)  nonnegative
     
     minimize -1*gamma* trace(I*Q) + 10 * sum(epsilon(1:num_upper_bounds)) + 10 * sum(zeta_l(1:num_lower_bounds)) + 10 * sum(zeta_vdw(1:num_vdw_bounds))
     subject to
     %equality constraint
       for k=1:num_equality_cons
           ti = equality_cons(k,1);
           tj = equality_cons(k,2);
           Mij = zeros(num_atoms);
           Mij(ti,tj) = -1;
           Mij(tj,ti) = -1;
           Mij(ti,ti) = 1;
           Mij(tj,tj) = 1;           
           trace(Mij*Q) == equality_cons(k,3)^2;
       end
       
      %upper bound constraint
       for k=1:num_upper_bounds
           ti = upper_bounds(k,1);
           tj = upper_bounds(k,2);
           Mij = zeros(num_atoms);
           Mij(ti,tj) = -1;
           Mij(tj,ti) = -1;
           Mij(ti,ti) = 1;
           Mij(tj,tj) = 1; 
           %trace(Mij*Q) <= upper_bounds(k,3)^2;
           trace(Mij*Q) - epsilon(k) <= upper_bounds(k,3)^2;
       end
       
      %lower bound constraint
       for k=1:num_lower_bounds
           ti = lower_bounds(k,1);
           tj = lower_bounds(k,2);
           Mij = zeros(num_atoms);
           Mij(ti,tj) = -1;
           Mij(tj,ti) = -1;
           Mij(ti,ti) = 1;
           Mij(tj,tj) = 1;
           %trace(Mij*Q) >= lower_bounds(k,3)^2;
           trace(Mij*Q) + zeta_l(k) >= lower_bounds(k,3)^2;
       end
       
      %van Der Waals lower constraint
       for k=1:num_vdw_bounds
           ti = vdw_bounds(k,1);
           tj = vdw_bounds(k,2);
           Mij = zeros(num_atoms);
           Mij(ti,tj) = -1;
           Mij(tj,ti) = -1;
           Mij(ti,ti) = 1;
           Mij(tj,tj) = 1;
           %trace(Mij*Q) >= vdw_bounds(k,3)^2;
           trace(Mij*Q) + zeta_vdw(k) >= vdw_bounds(k,3)^2;
       end
       
       %constraint Z1 = 0
         eeT = ones(num_atoms);
         trace(Q * eeT) == 0;
  cvx_end

  
  if ~issymmetric(Q)
      fprintf('The gram matrix found in NOT symmetric. Making it symmetric');
      Q = 1/2 .* (Q+Q');
  end

%   nDimensions = 3;
%          fprintf('\n Rank of the gram matrix',
%   [X,eigens]    = factorizematrix(Q,nDimensions);
  
Z = Q;

[VZ, LZ] = eigb(Z);
LZ = real(LZ);
LZ(LZ < 0) = 0;
X = LZ^0.5*VZ';
eigens = diag(LZ);

dLZ = diag(LZ);
ind = dLZ/max(dLZ); 
fprintf('\tRank: %d\n',sum(ind > 1e-6));

slacks=[];
% slacks(1) = epsilon;
% slacks(2) = zeta_l;
% slacks(3) = zeta_vdw;

end
  