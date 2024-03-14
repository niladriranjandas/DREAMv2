function [resi_present, resi_abscent] = findMissingResi(reg_atom_index,Comp_residue,...
                                     Comp_atom_names)
%% Function to get the gaps in the reconstructed structure after GRET
%   Input: reg_atom_index: include index for the atom numbers in the 
%                           reconstructed structure.
%            Comp_residue: 
%          Comp_atom_names:
%  Output: resi_present: residues in the reconstructed structure
%          resi_abscent: missing esidues in the reconstructed structure.
%%
  if length(Comp_residue) ~= length(Comp_atom_names)
      error('Module:findMissingResi: size of Comp.residue and Comp.atom_type does not match');
  end
  
%%
   all_resi_list = unique(Comp_residue);
   
   
   %% return resi having full backbone atoms in the reconstructed structure
   count = 0;
   backbone_atom = {'N','CA','C','O'};     % can be changed to include 'CB'
     recon_comp_residue = Comp_residue(reg_atom_index);
     recon_resi         = unique(recon_comp_residue);
      for i=1:length(recon_resi)
          ind = find(recon_comp_residue == recon_resi(i));
          recon_atom_name = Comp_atom_names(reg_atom_index(ind));
          flag=1;
          for j=1:length(backbone_atom)
              tmp = strcmp(backbone_atom(j),recon_atom_name);
              if ~any(tmp)
                 flag=0;
                  break
              end
          end
          if flag
              count = count+1;
              resi_present(count) = recon_resi(i);              
          end
      end
   
   %%
   
   resi_abscent = setdiff(all_resi_list, resi_present);
   resi_abscent = resi_abscent';
   
end   
   
   