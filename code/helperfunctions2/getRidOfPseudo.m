function [map,chk_coord] = getRidOfPseudo(coord,atom_seq,Compinfo)
%% filter out the pesudo atoms and return the rest along with the mapped index

 if nargin ~= 3
     error('Module: getRidOfPseudo: incorrect no. of inputs.')
 end
 
 if size(coord,1) ~= 3
   if size(coord,2) ~=3
       error('Module :getRidOfPseudo: incorrect coordinates inputs.');
   else
      num_atoms = size(coord,1);
   end
 else
    coord = coord'; 
    num_atoms = size(coord,1);
 end
 
 if num_atoms ~= length(atom_seq)
     error('Module :getRidOfPseudo: Incompatible coord and atom_seq.');
 end
     
%%
   skip_index_    = find(Compinfo(5,atom_seq)==0);
   skip_index     = atom_seq(skip_index_);
   map            = setdiff(atom_seq,skip_index);
   include_index  = setdiff(1:length(atom_seq),skip_index_);
   
   chk_coord = coord(include_index,:);   
   
end