%function densityreport = densityReport(file_path)
function densityreport = densityReport(wh_Comp,Comp,localize,protein_name,wh_up_bounds,hydrogen_omission)
%load(file_path);
%load('stage2_1.mat');

densityreport.protein_name = protein_name;

natoms = length(wh_Comp.atom_names);
n_wh_up_bounds = length(wh_up_bounds);

densityreport.fulldensity = 2*n_wh_up_bounds/(natoms * (natoms-1));

index_frgs = find([localize.method]==1);
densityreport.density_frags = nan(3,length(index_frgs));

    for i = 1:length(index_frgs)
        if hydrogen_omission
           resi_frags = findResidues(localize(index_frgs(i)).atoms, Comp.residue);
           resi_frags_noexp_nojump = findResidues(localize(index_frgs(i)).atoms_resi, Comp.residue);
        else
            resi_frags = findResidues(localize(index_frgs(i)).atoms, wh_Comp.residue);
            resi_frags_noexp_nojump = findResidues(localize(index_frgs(i)).atoms_resi, wh_Comp.residue);
        end
    
        atoms_frags = findAtomsFromResi(resi_frags, wh_Comp.residue);
        atoms_frags_noexp_nojump = findAtomsFromResi(resi_frags_noexp_nojump, wh_Comp.residue);
    
        frags_wh_up_bounds = extractCons(wh_up_bounds, atoms_frags);
        frags_wh_up_bounds_noexp_nojump = extractCons(wh_up_bounds, atoms_frags_noexp_nojump);
    
        n_atoms_frags = length(atoms_frags);
        n_wh_up_bounds_frags = length(frags_wh_up_bounds);
        n_atoms_frags_noexp_nojump = length(atoms_frags_noexp_nojump);
        n_wh_up_bounds_frags_noexp_nojump = length(frags_wh_up_bounds_noexp_nojump);
    
        densityreport.density_frags(1,i) = index_frgs(i);
        densityreport.density_frags(2,i) = n_wh_up_bounds_frags*2/(n_atoms_frags*(n_atoms_frags-1));
        densityreport.density_frags(3,i) = n_wh_up_bounds_frags_noexp_nojump*2/(n_atoms_frags_noexp_nojump*(n_atoms_frags_noexp_nojump-1));
    end

end

function resi = findResidues(atoms, Comp_residue)
     resi = unique(Comp_residue(atoms));
end


function atoms = findAtomsFromResi(residues, wh_Comp_residue)
    atoms = [];
    for i = 1:length(residues)
       indx = find(wh_Comp_residue == residues(i));
       atoms = [atoms;indx];
    end
end
