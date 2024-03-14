function pdb_coord = getCorrPDBcoord(PDBstruct_atom, chk_resi, chk_atoms)
%% function which returns the pdb coordinates given the
%  Input:  PDBStruct_atom: PDBStruct.model.atom   PDBStruct=pdbread('..')
%                 chk_res: list of the residue no. of the atoms to be
%                          checked
%                chk_atom: Seq of atoms identifiers such as N, CA, C,O etc
%                          (will no include pesudo atoms)
% Output:       pdb_coord: no. of atoms * 3 matrix having (x,y,z) coords
%%

  if nargin ~= 3
       error('Module: getCorrPDBcoord: incorrect number of inputs.');
  end
  
  if length(chk_resi) ~= length(chk_atoms)
      error('Module: getCorrPDBcoord: chk_resi and chk_atoms not compatible.');
  end
  
%%
  num_atoms = length(chk_resi);
  pdb_coord = zeros(num_atoms,3);
  
  for i=1:num_atoms
      pdb_atoms = find([PDBstruct_atom.resSeq]== chk_resi(i));
      
      for j=pdb_atoms
         if strcmp(PDBstruct_atom(j).AtomName,chk_atoms(i))
            pdb_coord(i,:) =  [PDBstruct_atom(j).X, PDBstruct_atom(j).Y, PDBstruct_atom(j).Z];
            break;
         end
      end
  end

end