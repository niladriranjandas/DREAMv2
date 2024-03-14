function [seqMat] = modifyCI_v2(cI,Comp)

%% dense subgraph algo creates clusters, such that points inside are in scattered clumps
%  This function separates those clumps
%%
 
 jump=100;
 
   count = 0;
   for i=1:length(cI)       
       %%detect peaks
       consec_diff = cI{i}(2:end) - cI{i}(1:end-1);
       [~,col,~] = find(consec_diff>jump);
       
       atom_no_start = cI{i}(1);
       seq_mat_row = [];
       for j=1:length(col)
         count = count+1;
         
         atom_no_end = cI{i}(col(j));
       %% find the residue no. corr to the atom no.
         seq_start = find(Comp.residue_bias <= atom_no_start);
         seq_start = seq_start(end);

         % seq_end   = find(Comp.residue_bias >= atom_no_end);
         % seq_end   = seq_end(1)
         seq_end   = find(Comp.residue_bias <= atom_no_end);
         seq_end   = seq_end(end);
         
         atom_no_start = cI{i}(col(j)+1);
       
       %% create the seqMat matrix
         seq_mat_row = [seq_mat_row;
                        seq_start,seq_end];
       
       end
      %%assign the last seq in the group
       atom_no_end = cI{i}(end);
       seq_start = find(Comp.residue_bias <= atom_no_start);
       seq_start = seq_start(end);
       
       % seq_end   = find(Comp.residue_bias >= atom_no_end);
       % seq_end   = seq_end(1)
       seq_end   = find(Comp.residue_bias <= atom_no_end);
       seq_end   = seq_end(end);
       
       seq_mat_row = [seq_mat_row;
                      seq_start,seq_end];
       
      seqMat{i} = seq_mat_row;
   end

end