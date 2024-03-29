% This function checks correct chirality
%
% ILE: res: 10 atoms = {CB, CA, CG1, CG2}
% THR: res: 17 atoms = {CB, CA, OG1, CG2}

function out = chirality_check_sidechain(X,Comp,print_flag)

if ~exist('print_flag','var')
    print_flag = 0;
end

any_error = 0;

seq = Comp.seq;
num_seq = Comp.num_seq;
atom_names = Comp.atom_names;

num_res = length(seq) - 1;
out = true(num_res,1);

for i = 1:num_res
   res = seq(i);
   if res == 10 || res == 17
      if res == 10
          chiral_atoms = {'CB','CA','CG1','CG2'};
      else
          chiral_atoms = {'CB','CA','OG1','CG2'};
      end
      chiral_atoms_index = nan(1,4);
      
      temp_atoms = Comp.residue == num_seq(i);
      temp_X = X(:,temp_atoms);
      temp_atom_names = atom_names(temp_atoms);
      for j = 1:4
          for k = 1:length(temp_atom_names)
             if strcmp(chiral_atoms{j},temp_atom_names{k})
                 chiral_atoms_index(j) = k;
                 break;
             end
          end          
      end     
      chir_X = temp_X(:,chiral_atoms_index);
      chir_angle = dicalc(chir_X);
      if chir_angle < 0
          out(i) = false;
          any_error = 1;
          if print_flag
             fprintf('- residue %3d has incorrect chirality\n',num_seq(i)); 
          end
      end
   end   
end

% if ~any_error
%     fprintf('\tNone\n');
% end