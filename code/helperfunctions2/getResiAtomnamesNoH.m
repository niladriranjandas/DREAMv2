function resi_noH_details = getResiAtomnamesNoH(A)
%% One time run code to get the atom names without H
%   Output:  resi_noH_details: cell type array 1 x 21
%            e.g.   resi_noH_details(1) = {'ALA' , {'N','CA','OB','CB','C','O'}}
%           (ordering of amino acid residues are same as 'A')
%    Input:  A : rotamer libr (same as CYANA 2.1  $cyanalib
%%   

   if nargin ~= 1
       error('Module: getResiAtomnamesNoH: incorrect no. of inputs.');
   end
%%
   num_resi = length(A);   
   resi_noH_details = cell(num_resi,2);
   reg_exp = '^[H|Q|h|q]';  %no H or pseudoatom
   
   
    tmp = cell.empty;
    for i=1:num_resi
        atom_names_resi_i = {A(i).atom.name};
        resi_noH_details(i,1) = {A(i).name};
        tmp = getnameNoH(atom_names_resi_i, reg_exp);
        resi_noH_details{i,2} = tmp;
    end
end

function atom_names = getnameNoH(cell_atom_names, reg_exp)
   len = length(cell_atom_names);
   atom_names = cell.empty;
     for i=1:len
        if isempty(regexp(char(cell_atom_names(i)),reg_exp))
            atom_names(i) = cell_atom_names(i);
        end
     end
    atom_names = atom_names(~cellfun('isempty',atom_names));
end