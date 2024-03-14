%% Post registration
%    -get regions for which all back bone atom coords are complete
%     * check if side chain (excluding H) exist for these regions 
%     * if not fill them (scwrl)
%    -prepare seq file for gap fill (modeller)
%     * generate the seq file
%     * generate the .ali file
%%     get the regions for which back bone is complete
clearvars fill ztmp

LOG.div_conquer = strcat(LOG.div_n_conquer,filesep,'div_n_conquer_',protein_name,'.txt');
LOG.prev_txt = '';

fill.resi_noH_details = getResiAtomnamesNoH(A);

[fill.bk_bone_tmp_ind, fill.bk_bone_atm_map] = onlyTheBackBone(reg.chk_coord_ref',reg.include_index,Comp);

fill.reg_resi_bkbn = unique(Comp.residue(fill.bk_bone_atm_map));

fill.tmp = Comp.residue(reg.include_index);
fill.reg_bkbn_atom_map=[];
for i=1:length(fill.reg_resi_bkbn)
    fill.indx = find(fill.tmp == fill.reg_resi_bkbn(i));
    fill.reg_bkbn_atom_map = [fill.reg_bkbn_atom_map, reg.include_index(fill.indx)'];
end
%%      get the regions which have complete side chains (excluding H) after prev step
fill.resi_full_sidec = getResiFullSidec(fill.reg_bkbn_atom_map, Comp, fill.resi_noH_details);

%%
ztmp.resi_present = fill.resi_full_sidec;
ztmp.count = 0;
 for i=1:length(ztmp.resi_present)
      ztmp.Comp_num_seq_indx =    find(Comp.num_seq == ztmp.resi_present(i));
      ztmp.tmp_resi = Comp.seq(ztmp.Comp_num_seq_indx);
      ztmp.count = ztmp.count + length(fill.resi_noH_details{ztmp.tmp_resi,2});
 end

%%      write the residues with complete backbone and sidechain (this are anchors)
  fill.Z_ref = nan(length(fill.reg_bkbn_atom_map),3);
  for i=1:length(fill.reg_bkbn_atom_map)
     fill.Z_ref(i,:) = reg.chk_coord_ref(find(reg.include_index==fill.reg_bkbn_atom_map(i)),:);
  end
  
  fill.bckbnfile = sprintf('grp_%s.pdb',protein_name);
  fill.filename = strcat(OUTPUT.local_folder_pdb, filesep, fill.bckbnfile);
  writeToPDB(fill.filename,fill.reg_bkbn_atom_map,fill.Z_ref',Comp);        % pdb read by modeller
  
  LOG.curr_txt = sprintf('\n Written only residues with complete backbone in grp_%s.pdb\n',protein_name);
  LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);  
%% -------------- create .ali file for modeller -----------------------

   fill.count = 0;
   
   %fill.seq = repmat('-',1,max_res);
   %fill.seq = repmat('-',1,max(org_num));
   fill.seq = repmat('-',1,length(org_num));   % 16-feb-20
         
   for i=1:length(ztmp.resi_present)
       ztmp.Comp_num_seq_indx =    find(Comp.num_seq == ztmp.resi_present(i));
       ztmp.tmp_resi = Comp.seq(ztmp.Comp_num_seq_indx);
       %fill.seq(ztmp.resi_present(i)) = A(ztmp.tmp_resi).code;  16-feb-20
       fill_seq_indx=find(org_num==ztmp.resi_present(i));        %16-feb-20
       fill.seq(fill_seq_indx) = A(ztmp.tmp_resi).code;          %16-feb-20
   end
   
   fill.ali_script_filename = strcat(OUTPUT.modeller_scripts, filesep, 'alignment_',protein_name,'.ali');
   fill.ali_script_fid      = fopen(fill.ali_script_filename,'w');
   
   if fill.ali_script_fid ~= -1
        fprintf(fill.ali_script_fid,'>P1;grp_%s\n',protein_name);
        %fprintf(fill.ali_script_fid,'structureX:grp_%s:%4d : :+%-5d: :::-1.00:-1.00\n',protein_name,fill.reg_resi_bkbn(1),max_res);
        fprintf(fill.ali_script_fid,'structureX:grp_%s:%4d : :+%-5d: :::-1.00:-1.00\n',protein_name,fill.reg_resi_bkbn(1),length(Comp.num_seq));
        fill.count=0;   
      %  for i=1:length(fill.seq)
      %  for i=1:length(num)
        for i=1:length(org_num)
            fill.count = fill.count+1;
            if fill.count > 75
                fill.count = 1;
                fprintf(fill.ali_script_fid,'\n');
            end
            %fprintf(fill.ali_script_fid,'%s',fill.seq(i));
            %fprintf(fill.ali_script_fid,'%s',fill.seq(num(i)));
            %fprintf(fill.ali_script_fid,'%s',fill.seq(org_num(i))); %16-feb-20
            z_index=i;%find(org_num==org_num(i));                    %16-feb-20
            fprintf(fill.ali_script_fid,'%s',fill.seq(z_index));     %16-feb-20
        end    
        fprintf(fill.ali_script_fid,'*\n');      
   
        fprintf(fill.ali_script_fid,'\n>P1;grp_%s_fill\nsequence:::::::::\n',protein_name);
        [seq_org, ~] = seq_reader(seq_file);    
        fill.count=0;
        for i=1:length(seq_org)
                fill.count = fill.count+1;
                if fill.count > 75
                    fill.count = 1;
                    fprintf(fill.ali_script_fid,'\n');
                end
                fprintf(fill.ali_script_fid,'%s',A(seq_org(i)).code);
        end
      fprintf(fill.ali_script_fid,'*');
      fclose(fill.ali_script_fid);
      
      LOG.curr_txt = sprintf('\n Written alignment_%s.ali file.\n',protein_name);
      LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
   else
        LOG.curr_txt = sprintf('Module: prepForFillGap : alignment_%s.ali file could not be created.',protein_name); 
        LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
        error('Module: prepForFillGap : alignment_%s.ali file could not be created.',protein_name); 
   end
%% ---------- create do_loop_dntmove.py for modeller ----------------------

  %fill.resi_abscent = setdiff(1:max_res, fill.reg_resi_bkbn);
  %fill.resi_abscent = setdiff(Comp.num_seq, fill.reg_resi_bkbn);
  fill.resi_abscent = setdiff(org_num, fill.reg_resi_bkbn);
  fill.tmp_ = fill.resi_abscent(2:end) - fill.resi_abscent(1:end-1);
  
  fill.tmp  = find(fill.tmp_ > 1);  fill.curr = 1;
  
  fill.py_script_filename = strcat(OUTPUT.modeller_scripts, filesep, 'do_loop_dntmove_',protein_name,'.py');
  fill.py_script_fid      = fopen(fill.py_script_filename,'w');
  
  if fill.py_script_fid ~= -1
       fprintf(fill.py_script_fid,'\nfrom modeller import *\nfrom modeller.automodel import *    # Load the automodel class');
       fprintf(fill.py_script_fid,'\nlog.verbose()\nenv = environ()\n');
       fprintf(fill.py_script_fid,'\n# directories for input atom files\nenv.io.atom_files_directory = [''.'', ''../atom_files'']\n');
       
       %fprintf(fill.py_script_fid,'class MyModel(automodel):\n\tdef select_atoms(self):\n\t\treturn selection(');
       %if num(1) == 1
       if org_num(1) == 1
         %fprintf(fill.py_script_fid,'class MyModel(automodel):\n\tdef select_atoms(self):\n\t\treturn selection(');
         fprintf(fill.py_script_fid,'class MyModel(automodel):\n\tdef special_patches(self, aln):\n');
         %fprintf(fill.py_script_fid,'\t\tself.rename_segments(segment_ids='''', renumber_residues=%d)\n',Comp.num_seq(1));
         fprintf(fill.py_script_fid,'\t\tself.rename_segments(segment_ids='''', renumber_residues=%d)\n',org_num(1));
         fprintf(fill.py_script_fid,'\tdef select_atoms(self):\n\t\treturn selection(');
       else
         fprintf(fill.py_script_fid,'class MyModel(automodel):\n\tdef special_patches(self, aln):\n');
         %fprintf(fill.py_script_fid,'\t\tself.rename_segments(segment_ids='''', renumber_residues=%d)\n',Comp.num_seq(1));
         fprintf(fill.py_script_fid,'\t\tself.rename_segments(segment_ids='''', renumber_residues=%d)\n',org_num(1));
         fprintf(fill.py_script_fid,'\tdef select_atoms(self):\n\t\treturn selection(');
       end
       for i=1:length(fill.tmp)
          if i==1
            fprintf(fill.py_script_fid,'self.residue_range(''%d'',''%d''),\n',fill.resi_abscent(fill.curr),fill.resi_abscent(fill.tmp(i)));   
          else
            fprintf(fill.py_script_fid,'\t\t\t\tself.residue_range(''%d'',''%d''),\n',fill.resi_abscent(fill.curr),fill.resi_abscent(fill.tmp(i)));
          end
          fill.curr = fill.tmp(i)+1;
       end
       fprintf(fill.py_script_fid,'\t\t\t\tself.residue_range(''%d'',''%d''))\n',fill.resi_abscent(fill.curr),fill.resi_abscent(end));
  
       fprintf(fill.py_script_fid,'\n\na = MyModel(env, alnfile = ''alignment_%s.ali'',\nknowns = ''grp_%s'', sequence = ''grp_%s_fill'')',protein_name,protein_name,protein_name);
       fprintf(fill.py_script_fid,'\na.starting_model= %d\na.ending_model  = %d',1,5);
       fprintf(fill.py_script_fid,'\n\na.make()\n');
       
       fclose(fill.py_script_fid);
       LOG.curr_txt = sprintf('\n Written do_loop_dntmove_%s.py file.\n',protein_name);
       LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
  else
      LOG.curr_txt = sprintf('Module: prepForFillGap: do_loop_dntmove_%s.py file creating failed.',protein_name);
      LOG.prev_txt = mexDolog(LOG.prev_txt,LOG.curr_txt,0,LOG.div_conquer);
      error('Module: prepForFillGap: do_loop_dntmove_%s.py file creating failed.',protein_name);
  end
  
LOG.prev_txt = mexDolog(LOG.prev_txt,'\n END of localization of groups',1,LOG.div_conquer);       
