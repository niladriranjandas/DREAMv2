function [X,atom_map,rot,L_inv,B] = doRegForGroup(grp_nums, x_n_index)
%% function to do GRET_SDP for groups no.(s) given by grp_nums
%   Input:  grp_nums: groups no.(s) corr to the "localize" structure to be
%                     registered
%            x_n_index: structure array of size #grp_nums  containing members
%                  x_n_index(i).x: 3xk coordinates for grp-i
%                x_n_index(i).ind: 1xk containing the index of the atoms in
%                                  grp-i
%  Output:         X: Registered coordinates 3 x no. of atoms(n)
%           atom_map: 1 x n array containing the atom index (original)
%                     indexed by temporary index. (e.g. [18 19 20 70 71..]
%                     temp index of atom-18 is 1 as so on.
%%
   if nargin ~=2 
       error('Module: doRegForGroup: not enough inputs.');
   end
   
   if length(grp_nums) ~= length(x_n_index)
       error('Module: doRegForGroup: grp_nums and structure x_n_index size incompatible.');
   end
   
   dim =3;
   n_grp = length(grp_nums);
%% assign tmp index to the atoms
   atom_map=[];
    for i=1:n_grp
       atom_map = union(atom_map,x_n_index(i).ind); 
    end   
  
%% build the structure (cell_x,index_x) as required by doGretSDP.m
   cell_x = cell(1,n_grp);   index_x = cell(1,n_grp);
   
   for i=1:n_grp
       cell_x(i)  = {x_n_index(i).x};
       index_x(i) = {mapToTmpIndex(x_n_index(i).ind,atom_map)};
   end
   
%% call global registration module
  [X,indices,rot,L_inv,B] = doGretSDP(cell_x,index_x,dim,[]);
  
end

function tmp_index = mapToTmpIndex(global_index,atom_map)
 
  n_atom = length(global_index);
  tmp_index = zeros(1,n_atom);
  for i=1:n_atom
     tmp_index(i) = find(atom_map==global_index(i)); 
  end

end

   
   