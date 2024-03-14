z_seq_st = 78; z_seq_nd = 90;

z_Comp = Comp;
[z_Comp,z_rand_X] = extractCompnRandX([z_seq_st,z_seq_nd],Comp,rand_X);

z_seq = seq(z_seq_st:z_seq_nd);
z_num = num(z_seq_st:z_seq_nd);

z_phi = phi(z_seq_st:z_seq_nd);
z_psi = psi(z_seq_st:z_seq_nd);

z_phi_cons = phi_cons(z_seq_st:z_seq_nd);
z_psi_cons = psi_cons(z_seq_st:z_seq_nd);

%%%%%%%%%%%%%%%%%%%%%%%
[z.rand_X, z.Comp] = ibuildprot(z_seq,z_num,z_phi,z_psi,A);
%%%%%%%%%%%%%%%%%%%%%%%%%

[z_U, z.Comp.cliq_dims] = reducer(z_rand_X,z.Comp);

        z_wh_Comp = z_Comp;
        patch_wh_eq_cons = equality_con_former(z_rand_X,z.Comp);
        patch_eq_cons = patch_wh_eq_cons;
        
        
atom_start = Comp.residue_bias(z_seq_st)+1;
atom_end   = Comp.residue_bias(z_seq_nd+1);   
 
patch_up_bounds = extractCons(up_bounds,atom_start:atom_end);     
patch_sdp_lo_bounds = extractCons(sdp_lo_bounds,atom_start:atom_end);
patch_lo_bounds = extractCons(lo_bounds,atom_start:atom_end);

[z_rawX, z_patch_eV, z_patch_eS] = solve_sdpt3(z_U,patch_eq_cons,...
                                                      patch_up_bounds,patch_sdp_lo_bounds,[],f);
        

                                                  

