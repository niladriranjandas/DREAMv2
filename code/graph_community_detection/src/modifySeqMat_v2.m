function modseqMat = modifySeqMat_v2(seqMat,adj_mat,residue_bias,tot_resi)
%% modify seqMat to change i,i pairs in seqMat to (i-1),i or i,(i+1)
%%
    if nargin ~=4
        error('Module modifySeqMat: Not enough inputs');
    end
    
    if ~issymmetric(adj_mat) && ~isuppertriangular(adj_mat)
        error('Module modifySeqMat: Adjacency matrix is incorrect');
    end 
    
%%  gather the atom seq belonging to (i,j) pair in other_atom and 
%   those belonging to (i,i) pair in atom
    modseqMat = seqMat;
 
    for i=1:length(seqMat)
      other_atom = [];
      resi=[];
      resi_pos = [];
      count = 0;
       for j=1:length(seqMat{i})
          if seqMat{i}(j,1) ~= seqMat{i}(j,2)
             st_atm = residue_bias(seqMat{i}(j,1))+1;
             if seqMat{i}(j,1) ~= tot_resi
                nd_atm = residue_bias(seqMat{i}(j,2)+1);                
             end
             other_atom = [other_atom,st_atm:nd_atm];
          else
             count = count +1;
             resi_pos(count) = j;
             resi(count) = seqMat{i}(j,1);
          end
       end
       
    %%select whether to include seqMat{i}(j-1) or seqMat{i}(j+1) in
    %%seqMat{i}(j,:)
       for j=1:length(resi)
     display('kjkjkj')
          %test for atoms belonging to resi(j) -1
          tmp_resi1 = resi(j) -1;
          st_atm1 = residue_bias(tmp_resi1)+1;
          if tmp_resi1 ~= tot_resi
            nd_atm1 = residue_bias(tmp_resi1+1);
          end
          atm_seq1 = st_atm1:nd_atm1;
          neigh_minus_1 = sum(sum(adj_mat(atm_seq1,other_atom),2))
          
          %test for atoms belonging to resi(j) +1
          tmp_resi2 = resi(j) +1;
          st_atm2 = residue_bias(tmp_resi2)+1;
          if tmp_resi2 ~= tot_resi
            nd_atm2 = residue_bias(tmp_resi2+1);
          end
          atm_seq2 = st_atm2:nd_atm2;
          neigh_plus_1 = sum(sum(adj_mat(atm_seq2,other_atom),2))
          
          %decide which one to include
          if neigh_minus_1 >= neigh_plus_1
             modseqMat{i}(resi_pos(j),:) = [modseqMat{i}(resi_pos(j),1)-1,modseqMat{i}(resi_pos(j),1)];
          else
             modseqMat{i}(resi_pos(j),:) = [modseqMat{i}(resi_pos(j),1),modseqMat{i}(resi_pos(j),1)+1];
          end
       end
       
    end
    
    
end    
        