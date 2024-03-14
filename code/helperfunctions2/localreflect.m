function reflect2 = localreflect(X_refine, resi_start, resi_end, Comp, atom_map, pdbnameprefix,pdb_op_dir,xyz_flag)
%%
%
%%
        if nargin==6
		xyz_flag=1;
	end
      
	beginpos_ = find(Comp.residue == resi_start);
	endpos_   = find(Comp.residue == resi_end);

	beginpos  = beginpos_(1);
	endpos    = endpos_(end);
	
	begin_pos_in_X  = find(atom_map == beginpos);
	end_pos_in_X    = find(atom_map == endpos);
%% --------------------------------------------------------------
	reflect2.x_n_index(1).x   = X_refine(:,[1:begin_pos_in_X+3,end_pos_in_X-3:end]);
	reflect2.x_n_index(1).ind = atom_map([1:begin_pos_in_X+3,end_pos_in_X-3:end]);

	reflect2.x_n_index(2).x   = X_refine(:,begin_pos_in_X:end_pos_in_X);
	reflect2.x_n_index(2).ind = atom_map(begin_pos_in_X:end_pos_in_X);

	[reflect2.reg_X_noref, reflect2.reg_atom_map, reflect2.reg_rot, reflect2.reg_L_inv, reflect2.reg_B] = doRegForGroup([1,2],reflect2.x_n_index);
	reflect2.reg_rot{1} = eye(3); 
    reflect2.reg_rot{2} = eye(3);
    reflect2.reg_X_noref = [reflect2.reg_rot{1},reflect2.reg_rot{2}]  * reflect2.reg_B * reflect2.reg_L_inv;
    reflect2.reg_X_noref = reflect2.reg_X_noref(:,1:length(reflect2.reg_atom_map));
    [reflect2.reg_include_index,reflect2.reg_chk_coord_ref] = getRidOfPseudo(reflect2.reg_X_noref',reflect2.reg_atom_map,Comp.info);                              
	pdb_file_name = strcat(pdb_op_dir,filesep,pdbnameprefix,'_tmp.pdb');%'../output_pdb/grp_registration_1owa_localreg2.pdb';
        reflect2.reg_resi_index={writeToPDB(pdb_file_name,reflect2.reg_include_index,reflect2.reg_chk_coord_ref',Comp)};

	[rama.out,fig_handler] = ramachandranMe(pdb_file_name,'Regions',true,'Glycine',true);

	%save ../protein/1owa/1owa_local_reflect2_afterconsolidate.mat

	reflect2.I = eye(3); reflect2.I(xyz_flag, xyz_flag) = -1; %reflect2.I(3,3) = -1; %reflect2.I(2,2) = -1;   %%% VVIP relect.I(1,1) = -1 give wrong 
	reflect2.reg_rot_new = [ reflect2.reg_rot{1}, reflect2.I * reflect2.reg_rot{2} ];

	reflect2.reg_tmp = size(reflect2.reg_X_noref,2);                                                   
	reflect2.reg_X_noref = reflect2.reg_rot_new * reflect2.reg_B * reflect2.reg_L_inv;
	reflect2.reg_X_noref = reflect2.reg_X_noref(:,1:reflect2.reg_tmp);

        [reflect2.reg_include_index,reflect2.reg_chk_coord_ref] = getRidOfPseudo(reflect2.reg_X_noref',reflect2.reg_atom_map,Comp.info);                              
	%pdb_file_name = '../output_pdb/grp_registration_1owa_localreg_reflect2.pdb';
        pdb_file_name = strcat(pdb_op_dir,filesep,pdbnameprefix,'_reflect.pdb');
        reflect2.reg_resi_index={writeToPDB(pdb_file_name,reflect2.reg_include_index,reflect2.reg_chk_coord_ref',Comp)};

	[rama.out,fig_handler] = ramachandranMe(pdb_file_name,'Regions',true,'Glycine',true);

	%save ../protein/1owa/1owa_local_reflect2_noref.mat

% ---------------------------------------------------------------	
end
