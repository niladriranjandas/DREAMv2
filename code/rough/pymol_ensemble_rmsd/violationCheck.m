%%
clearvars checkPDB

%%
cd ..
addpath(genpath(pwd));
cd code;
%------INPUTS and PARAMETERS-------

 INPUTS.protein_name = '2mc2';
 INPUTS.seq_file = '2mc2.seq';
 INPUTS.upl_file = {'2mc2.dist.1.upl'};
 INPUTS.hbond_file = '2mc2.hBond.1.upl';
 INPUTS.ang_file = '2mc2.dihed.1.aco';
 INPUTS.protein_path = '../protein/2mc2';

checkPDB.filename = '2m4k_1_ancloc_refine.pdb'
%checkPDB.filename = '2m4k_sa.pdb';

%%
preprocess

checkPDB.pdbStruct   = pdbread(checkPDB.filename);
checkPDB.Atoms       = checkPDB.pdbStruct.Model.Atom;
checkPDB.X           = [[checkPDB.Atoms.X]; [checkPDB.Atoms.Y]; [checkPDB.Atoms.Z]];
checkPDB.AtomsName   = {checkPDB.Atoms.AtomName};
checkPDB.residue     = [checkPDB.Atoms.resSeq];
checkPDB.residuename = {checkPDB.Atoms.resName};

%% check of 1st residue number is that of the min_res and max resi is max_res

%%

checkPDB.x_filled = fillInPseudo(checkPDB.X, checkPDB.residue, checkPDB.AtomsName, checkPDB.residuename, wh_Comp, A);
checkPDB.atom_map = find(~isnan(checkPDB.x_filled(1,:)));   

checkPDB.x_chk = checkPDB.x_filled(:,checkPDB.atom_map);
%% extract constraints
 checkPDB.eq_cons_atommap      = extractCons(wh_eq_cons_all, checkPDB.atom_map);
 checkPDB.wh_up_bounds_atommap = extractCons(wh_up_bounds, checkPDB.atom_map);
 checkPDB.wh_lo_bounds_atommap = extractCons(wh_lo_bounds, checkPDB.atom_map);
 
%%
 [violation.mean_err, violation.min_err, violation.max_err, violation.err_percent] = calcViolations(checkPDB.x_chk, checkPDB.eq_cons_atommap ,3);
 fprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        
 [violation.mean_err, ~, violation.max_err, violation.err_percent] = calcViolations(checkPDB.x_chk, checkPDB.wh_up_bounds_atommap,1);
 fprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
 
 [violation.mean_err, violation.min_err, ~, violation.err_percent] = calcViolations(checkPDB.x_chk, checkPDB.wh_lo_bounds_atommap,2);
 fprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
 