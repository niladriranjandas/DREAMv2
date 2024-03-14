Aloc.tmp_wh_Comp.seq        = wh_Comp.seq;
Aloc.tmp_wh_Comp.num_seq    = wh_Comp.num_seq;
Aloc.tmp_wh_Comp.atom_names = wh_Comp.atom_names(Aloc.atom_map);
Aloc.tmp_wh_Comp.residue    = wh_Comp.residue(Aloc.atom_map);
Aloc.tmp_wh_Comp.atom_types = wh_Comp.atom_types(Aloc.atom_map);
Aloc.tmp_wh_Comp.residue_bias = [0, find(Aloc.tmp_wh_Comp.residue(2:end) - Aloc.tmp_wh_Comp.residue(1:end-1))'];

Aloc.chiral_rep = chirality_check(fillQ.X_refine, Aloc.tmp_wh_Comp, 1);%0);
Aloc.x_chiral   = chirality_correction(fillQ.X_refine,Aloc.tmp_wh_Comp,Aloc.chiral_rep);

%% calc violations
fprintf('\n=================================== Violations (Stage: Post anchored localization): chiral corr =========================================');
%         LOG.curr_txt = sprintf('\n=================================== Violations (Stage: Anchored localization) =========================================');
%         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
         
%         fprintf('\n Before refinement');
%         LOG.curr_txt = sprintf('\n Before refinement');
%         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
         
         [violation.mean_err, violation.min_err, violation.max_err, violation.err_percent] = calcViolations(Aloc.x_chiral,Aloc.eq_cons_atommap ,3);
         fprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
%         LOG.curr_txt = sprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
%         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
        
         [violation.mean_err, ~, violation.max_err, violation.err_percent] = calcViolations(Aloc.x_chiral,Aloc.wh_up_bounds_atommap,1);
         fprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
%         LOG.curr_txt = sprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
%         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);

         [violation.mean_err, violation.min_err, ~, violation.err_percent] = calcViolations(Aloc.x_chiral,Aloc.wh_lo_bounds_atommap,2);
         fprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
%         LOG.curr_txt = sprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
%         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
disp('')

%% refine

    [Aloc.x_chiral_ref,Aloc.info] = postProcessingMe(Aloc.x_chiral,...
                                                   fillQ.wh_eq_cons_all, fillQ.wh_lo_bounds, fillQ.wh_up_bounds,...
                                                   [1e4,1,1,-0.001],[10,10,10,10,10]);
                                               
figure;plot(1:length(Aloc.info.obj),fillQ.info.obj);hold on; title('refine');hold off                                               
                                               
%% calc violations
fprintf('\n=================================== Violations (Stage: Post anchored localization): chiral corr =========================================');
%         LOG.curr_txt = sprintf('\n=================================== Violations (Stage: Anchored localization) =========================================');
%         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
         
%         fprintf('\n Before refinement');
%         LOG.curr_txt = sprintf('\n Before refinement');
%         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
         
         [violation.mean_err, violation.min_err, violation.max_err, violation.err_percent] = calcViolations(Aloc.x_chiral_ref,Aloc.eq_cons_atommap ,3);
         fprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
%         LOG.curr_txt = sprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
%         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
        
         [violation.mean_err, ~, violation.max_err, violation.err_percent] = calcViolations(Aloc.x_chiral_ref,Aloc.wh_up_bounds_atommap,1);
         fprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
%         LOG.curr_txt = sprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
%         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);

         [violation.mean_err, violation.min_err, ~, violation.err_percent] = calcViolations(Aloc.x_chiral_ref,Aloc.wh_lo_bounds_atommap,2);
         fprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
%         LOG.curr_txt = sprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
%         LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.anchor);
disp('')

%% write pdb
{writeToPDB(sprintf('%s_anchor_%d.pdb',protein_name,OUTPUT.model_count),fillQ.atom_map,Aloc.x_chiral_ref,wh_Comp)};

%%

z(1).x_ref_chiral = fillQ.x_ref_chiral;
z(1).x_ref_chiral_info = nan;

for i=2:20 % 20 will be assigned a param say max_itter
    disp(i)
  [z(i).x_ref_chiral, z(i).x_chiral, z(i).chiral_info] = chiralChknCorr_rough(z(i-1).x_ref_chiral, fillQ.atom_map, fillQ.wh_up_bounds, fillQ.wh_lo_bounds, fillQ.wh_eq_cons_all , wh_Comp);
  if (i>3) && (z(i).chiral_info.err_main+z(i).chiral_info.err_side == z(i-2).chiral_info.err_main + z(i-2).chiral_info.err_side  )
     break
  end
end

for i=2:length(z)
       z_main(i-1) = z(i).chiral_info.err_main; z_side(i-1) = z(i).chiral_info.err_side;
end

plot(1:10,z_main,'k');hold on;plot(1:10,z_side,'y')
plot(1:10,z_main,'k');hold on;plot(1:10,z_side,'b')
hold on;
plot(1:10,z_main+z_side,'r');hold off
