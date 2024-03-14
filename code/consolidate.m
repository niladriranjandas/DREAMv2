% Consolidate step following divide-and-conquer
%   Globally register the fragments
%   requires package matgraph, bioinformatics toolbox in matlab
%%  prep the data to call doRegForGroup
  
clearvars reg x_n_index rama %index 

     local_folder_plots_postreg = OUTPUT.local_folder_plots_postreg;
          if ~exist(local_folder_plots_postreg,'dir')
             mkdir(local_folder_plots_postreg);
          end
     reg.frag_graph_name = strcat(protein_name,'_fragment_graph');
     
  LOG.div_conquer = strcat(LOG.div_n_conquer,filesep,'div_n_conquer_',protein_name,'.txt');
  LOG.prev_txt = '';

  % --- changed on 24-Jun-2021 -- t avoid error due to fragments not localized --- %
  index_chk = find([localize.method]==1);

  for chk_i = 1:length(index_chk)
     [chk_r,~] = size(localize(index_chk(chk_i)).info.obj );
     if chk_r == 1
          localize(index_chk(chk_i)).method = -3;
     end
  end  
  % -------------------------------------------------------------------------------%
  
  
index = find([localize.method]==1);
   %reg.frag_biograph = biograph((test.patch_adj(index,index)>3).*(test.patch_adj(index,index)),cellstr(num2str(index(:))),'ShowWeights','on');   % 19-19-20 : no fig run
   %frag_tmp_ = biograph.bggui(reg.frag_biograph); % 19-19-20 : no fig run
   %reg.frag_fig = figure; % 19-19-20 : no fig run
   %copyobj(frag_tmp_.biograph.hgAxes,reg.frag_fig); % 19-19-20 : no fig run
   %print(reg.frag_fig,'-depsc','-r300',strcat(local_folder_plots_postreg,filesep,reg.frag_graph_name)); % 19-19-20 : no fig run
   
%%
   %--------- which components to register in the patch graph-----------------
reg.g=graph; sparse(reg.g);
reg.patch_graph = (test.patch_adj(index,index)>3).*(test.patch_adj(index,index));  % consider only patches/fragments with >3 common points 
set_matrix(reg.g,(reg.patch_graph+reg.patch_graph')>0);
label(reg.g,cellstr(num2str(index(:))));
reg.conc_compo  = components(reg.g);
reg.parts_compo_tmp = parts(reg.conc_compo);                                       % connected components in the patch/fragment graph
reg.parts_compo     = cell(1,length(reg.parts_compo_tmp));


for i=1:length(reg.parts_compo_tmp)
    reg.parts_compo(i) = {index(reg.parts_compo_tmp{i})};
end

[~,reg.I] = max(cellfun(@length, reg.parts_compo));
reg.grp_nums = reg.parts_compo{reg.I};    % MUST INCLUDE LOGIC TO CHOOSE THIS (perhaps for all the connected components)
%reg.grp_nums = [1,3,4,6];

% [reg.grp_nums',[localize(reg.grp_nums).rmsd_ref]',[localize(reg.grp_nums).rmsd_noref]']
disp('---------------');
    
  for i=1:length(reg.grp_nums)
      x_n_index(i).x   = localize(reg.grp_nums(i)).X_refine;
      x_n_index(i).ind = localize(reg.grp_nums(i)).atoms;
  end

%% logs 
    fprintf('\n Connected components %d\n',length(reg.parts_compo));
    for i=1:length(reg.parts_compo)
       fprintf('\n \t size of component %d : %d',i,length(reg.parts_compo{i}));
    end
    fprintf('\n');
%% global register
time_consolidate = tic;

fprintf('\n__________global registration module___\n');  
LOG.curr_txt = sprintf('\n__________global registration module___\n');  
LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
  [reg.X_noref, reg.atom_map, reg.rot, reg.L_inv, reg.B] = doRegForGroup(reg.grp_nums,x_n_index);
           for i=1:length(reg.rot)     %CHECK sometimes gotto multiply with -1 the det of R's
              det(reg.rot{i})          %CHECK sometimes gotto multiply with -1 the det of R's
           end                         %CHECK sometimes gotto multiply with -1 the det of R's
  
%% refine
   %--------------------get the constraints--------------------------------
   reg.eq_cons_all   = extractCons(eq_cons_all,reg.atom_map);
   reg.up_bounds     = extractCons(up_bounds,reg.atom_map);
   reg.wh_vdw_bounds = extractCons(wh_vdw_bounds,reg.atom_map);
   reg.vdw_bounds    = extractCons(vdw_bounds,reg.atom_map);
   reg.sdp_lo_bounds = extractCons(sdp_lo_bounds,reg.atom_map);
    
   reg.lo_bounds     = extractCons(lo_bounds,reg.atom_map);
   reg.wh_lo_bounds  = extractCons(wh_lo_bounds,reg.atom_map);
   
   reg.count = 1;
   while 1
      if reg.count > 2
          break
      end
        %--------------------call the refine module-----------------------------
        fprintf('\n__________refine module___\n');
%         [reg.X_refine,reg.info] = postProcessingMe(reg.X_noref,...
%                                                    reg.eq_cons_all, reg.lo_bounds, reg.up_bounds,...
%                                                    [1e4,1,1,-0.001],[10,10,10,10,10]);
        [reg.X_refine,reg.info] = postProcessingMe(reg.X_noref,...
                                                   reg.eq_cons_all, reg.lo_bounds, full(reg.up_bounds),...
                                                   [1e4,1,1,-0.001],[10,10,10,10,10]);
        %---------------extract the original PDB coordinates---------------------
        [reg.include_index,reg.chk_coord_ref] = getRidOfPseudo(reg.X_refine',reg.atom_map,Comp.info);
        [reg.include_index_noref,reg.chk_coord_noref] = getRidOfPseudo(reg.X_noref',reg.atom_map,Comp.info);                                          

%% write logs and violations

        fprintf('\n=================================== Violations (Stage: Consolidate )==============================================');
        LOG.curr_txt = sprintf('\n=================================== Violations (Stage: Consolidate )==============================================');
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
    
        [violation.mean_err, violation.min_err, violation.max_err, violation.err_percent] = calcViolations(reg.X_noref,reg.eq_cons_all,3);
        fprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);            
    
        [violation.mean_err, ~, violation.max_err, violation.err_percent] = calcViolations(reg.X_noref,reg.up_bounds,1);
        fprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);        
        LOG.curr_txt = sprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);        
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);

        [violation.mean_err, violation.min_err, ~, violation.err_percent] = calcViolations(reg.X_noref,reg.sdp_lo_bounds,2);
        fprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);

        fprintf('\n\n After refinement');
        LOG.curr_txt = sprintf('\n\n After refinement');
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
        
        [violation.mean_err, violation.min_err, violation.max_err, violation.err_percent] = calcViolations(reg.X_refine,reg.eq_cons_all,3);
        fprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.max_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
        
        [violation.mean_err, ~, violation.max_err, violation.err_percent] = calcViolations(reg.X_refine,reg.up_bounds,1);
        fprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', violation.mean_err, violation.max_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);

        [violation.mean_err, violation.min_err, ~, violation.err_percent] = calcViolations(reg.X_refine,reg.sdp_lo_bounds,2);
        fprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
        LOG.curr_txt = sprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', violation.mean_err, violation.min_err, violation.err_percent);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
    
        fprintf('\n'); 
%% prepare pdb file            
        fprintf('\n_________write to pdb file: registration result_________\n');  
        %pdb_file_name = 'grp_registration.pdb';
        pdb_file_name = sprintf('grp_registration_%s.pdb',protein_name);  pdb_file = strcat(OUTPUT.local_folder_pdb,filesep,pdb_file_name);
        reg.resi_index={writeToPDB(pdb_file,reg.include_index,reg.chk_coord_ref',Comp)};
        
        LOG.curr_txt = sprintf('\n %s written.\n',pdb_file);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);

 %% save the ramachandran plot and info about the resi in the various region
 
        plot_filename_rama = sprintf('%s_ramachandran_run%d.fig',protein_name,reg.count);
        plot_file_rama    = strcat(local_folder_plots_postreg,filesep,plot_filename_rama);
 
        %[rama.out,fig_handler] = ramachandranMe(pdb_file,'Regions',true,'Glycine',true);   % 19-09-20 no fig run
        [rama.out,~] = ramachandranMe(pdb_file,'Regions',true,'Glycine',true, 'Plot','None'); % 19-09-20 no fig run
         
        %fprintf('\n_________write the ramachandran plot_________\n'); % 19-09-20 no fig run
        %savefig(fig_handler,plot_file_rama);                          % 19-09-20 no fig run
        
        LOG.curr_txt = sprintf('\n %s written.\n',plot_file_rama);
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
 
        %resi in various area in the ramachandran map
        [rama.region,rama.part_tot_resi] = getRamachandranReigionDistri(rama.out.Angles);
        
 %% get percent disallowed region for ramachandran map
 
     rama.disallowed = length(rama.region(end).resi) / length([rama.region.resi]);
     
     LOG.curr_txt = sprintf('\n Fraction of residues in disallowed region: %f \n',rama.disallowed);  
     LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
     
     
     if rama.disallowed >= 0.4 %0.5  BUG REPORT: DONE AFTER 1TFB
        reg.count = reg.count+1;
        fprintf('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Reflecting the rotation matrices for ramachandran correction !!!!!!!!!!!!!!!!!!!!!!!\n');
        LOG.curr_txt = sprintf('\n!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Reflecting the rotation matrices for ramachandran correction !!!!!!!!!!!!!!!!!!!!!!!\n');
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
        reg.I = eye(3); reg.I(1,1) = -1; 
        reg.tmp = size(reg.X_refine,2);
        reg.X_noref = reg.I * [reg.rot{1:end}] * reg.B * reg.L_inv;
        reg.X_noref = reg.X_noref(:,1:reg.tmp);        
     else
        break
     end
   end       
time.time_consolidate= toc(time_consolidate);

LOG.curr_txt = sprintf('Module: Conquer: time taken (includes pdb, ramachandran and log file write): %d',time.time_consolidate);
LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
   
clearvars fig_handler frag_tmp_ 
LOG.prev_txt = mexDolog(LOG.prev_txt,'\n___________END of global registration module__________\n',1,LOG.div_conquer);
