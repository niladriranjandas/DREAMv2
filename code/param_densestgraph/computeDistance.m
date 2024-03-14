%%*****************************************************************
%% this function has been adopted from DISCO code ../DISCO/setup/computedistances.m
%% find all the pairs of atoms with distances < radius. 
%%
%%*****************************************************************

  function [Di,Dj,Dv] = computeDistances(A)

  [nDimensions,nAtoms] = size(A);

  Di = [];
  Dj = [];
  Dv = [];
  for i = 1:nAtoms
     m = nAtoms - i;
     AiAj = A * ...
           (sparse(i,1:m,1,nAtoms,m) - sparse((i+1):nAtoms,1:m,1,nAtoms,m));
     d = sqrt(sum(AiAj.^2))';
     idx = find(d);
     %idx = find(d < radius); 
     Di = [Di; i * ones(length(idx),1)];
     Dj = [Dj; i + idx];
     Dv = [Dv; d(idx)];
  end
%%*****************************************************************

