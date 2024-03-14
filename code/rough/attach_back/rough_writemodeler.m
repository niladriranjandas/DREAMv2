[part2.back_bone_tmp_index, part2.back_bone_atom_map] = onlyTheBackBone(reg2.chk_coord_ref',reg2.include_index,Comp);
part2.reg_resi_bkbn = unique(Comp.residue(part2.back_bone_atom_map));

part2.anchors_resi = unique(Comp.residue(part2.back_bone_atom_map));

part2.tmp = Comp.residue(reg2.include_index);
part2.reg_bkbn_atom_map=[];
for i=1:length(part2.reg_resi_bkbn)
    part2.indx = find(part2.tmp == part2.reg_resi_bkbn(i));
    part2.reg_bkbn_atom_map = [part2.reg_bkbn_atom_map, reg2.include_index(part2.indx)'];
end

%%      write the residues with complete backbone and sidechain (this are anchors)
  part2.Z_ref = nan(length(part2.reg_bkbn_atom_map),3);
  for i=1:length(part2.reg_bkbn_atom_map)
     part2.Z_ref(i,:) = reg2.chk_coord_ref(find(reg2.include_index==part2.reg_bkbn_atom_map(i)),:);
  end
  
  part2.bckbnfile = sprintf('grp_%s.pdb',protein_name);
  part2.filename = part2.bckbnfile;
  %part2.filename = strcat(OUTPUT.local_folder_pdb, filesep, part2.bckbnfile);
  writeToPDB(part2.filename,part2.reg_bkbn_atom_map,part2.Z_ref',Comp);        % pdb read by modeller

%% prepare alignment.ali
part2.fid = fopen('reg2_alignment.ali','w');
if part2.fid < 0
   error('error in openening the .ali file.');
   fclose(part2.fid);
   exit -1;
end

part2.len     = length(part2.anchors_resi(1) : part2.anchors_resi(end));
part2.part2txt = repmat('-',1,part2.len);

fprintf(part2.fid,'>P1;grp_%s\n',protein_name);
fprintf(part2.fid,'structureX:grp_%s:%4d : :+%-5d: :::-1.00:-1.00\n',protein_name,part2.anchors_resi(1),length(part2.anchors_resi));
count=0;
for i=1:length(part2.anchors_resi)   
    ztmp.Comp_num_seq_indx = find(Comp.num_seq == part2.reg_resi_bkbn(i));
    ztmp.tmp_resi = Comp.seq(ztmp.Comp_num_seq_indx);
    part2.part2txt(part2.anchors_resi(i)-part2.anchors_resi(1)+1) = A(ztmp.tmp_resi).code;
end

    for i=1:part2.len
        count=count+1;
        if count > 75
            count = 1;
            fprintf(part2.fid,'\n');
        end
        fprintf(part2.fid,'%s',part2.part2txt(org_num(i)));
    end
    fprintf(part2.fid,'*\n');
    
    fprintf(part2.fid,'\n>P1;grp_%s_part2\nsequence:::::::::\n',protein_name);
    [seq_org, ~] = seq_reader('/home/niladri/Documents/Disco_etc_all_in_1/our_algo/protein/2kul/2kul.seq');    
    count=0;
    for i=part2.anchors_resi(1):part2.anchors_resi(end)
        count = count+1;
        if count > 75
             count = 1;
             fprintf(part2.fid,'\n');
         end
         fprintf(part2.fid,'%s',A(seq_org(i)).code);
    end
    fprintf(part2.fid,'*');        

fclose(part2.fid);    
%fprintf(part2.fid,'\n');    
%% prepare the do_loop_dntmove.py
part2.fid = fopen('reg2_do_loop_dntmove.py','w');
if part2.fid < 0
   error('error in openening the .ali file.');
   fclose(part2.fid);
   exit -1;
end

part2.resi_abscent = setdiff(part2.anchors_resi(1) : part2.anchors_resi(end), part2.anchors_resi);
  part2.tmp_ = part2.resi_abscent(2:end) - part2.resi_abscent(1:end-1);
  
  part2.tmp  = find(part2.tmp_ > 1);  part2.curr = 1;
   
       fprintf(part2.fid,'\nfrom modeller import *\nfrom modeller.automodel import *    # Load the automodel class');
       fprintf(part2.fid,'\nlog.verbose()\nenv = environ()\n');
       fprintf(part2.fid,'\n# directories for input atom files\nenv.io.atom_files_directory = [''.'', ''../atom_files'']\n');
       
       %fprintf(part2.fid,'class MyModel(automodel):\n\tdef select_atoms(self):\n\t\treturn selection(');
       %if num(1) == 1
       if part2.anchors_resi(1) == 1
         fprintf(part2.fid,'class MyModel(automodel):\n\tdef select_atoms(self):\n\t\treturn selection(');
       else
         fprintf(part2.fid,'class MyModel(automodel):\n\tdef special_patches(self, aln):\n');
         %fprintf(part2.fid,'\t\tself.rename_segments(segment_ids='''', renumber_residues=%d)\n',Comp.num_seq(1));
         fprintf(part2.fid,'\t\tself.rename_segments(segment_ids='''', renumber_residues=%d)\n',part2.anchors_resi(1));
         fprintf(part2.fid,'\tdef select_atoms(self):\n\t\treturn selection(');
       end
       for i=1:length(part2.tmp)
          if i==1
            fprintf(part2.fid,'self.residue_range(''%d'',''%d''),\n',part2.resi_abscent(part2.curr),part2.resi_abscent(part2.tmp(i)));   
          else
            fprintf(part2.fid,'\t\t\t\tself.residue_range(''%d'',''%d''),\n',part2.resi_abscent(part2.curr),part2.resi_abscent(part2.tmp(i)));
          end
          part2.curr = part2.tmp(i)+1;
       end
       fprintf(part2.fid,'\t\t\t\tself.residue_range(''%d'',''%d''))\n',part2.resi_abscent(part2.curr),part2.resi_abscent(end));
  
       fprintf(part2.fid,'\n\na = MyModel(env, alnfile = ''reg2_alignment.ali'',\nknowns = ''grp_%s'', sequence = ''grp_%s_part2'')',protein_name,protein_name);
       fprintf(part2.fid,'\na.starting_model= %d\na.ending_model  = %d',1,5);
       fprintf(part2.fid,'\n\na.make()\n');
       
fclose(part2.fid);           