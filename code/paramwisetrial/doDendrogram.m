function doDendrogram(matfiles, opfile, flagfilestage1)
%% load the matfiles and doDendogram
% 
%%
    prefix = "../protein/";

    count = 0;
    files = split(matfiles, ',');    
    for j = 1:length(files)
        clearvars -except matfiles opfiles params count files j opfile prefix flagfilestage1
        load(strcat(prefix,filesep,files{j}));
        count = count +1;
        % --------------------------------%
            index =  find([test.dup_flag]==1);
            regtmp.g=graph; sparse(regtmp.g);
            regtmp.patch_graph = (test.patch_adj(index,index)>3).*(test.patch_adj(index,index));  % consider only patches/fragments with >3 common points 
            set_matrix(regtmp.g,(regtmp.patch_graph+regtmp.patch_graph')>0);
            label(regtmp.g,cellstr(num2str(index(:))));
            regtmp.conc_compo  = components(regtmp.g);
            regtmp.parts_compo_tmp = parts(regtmp.conc_compo);                                       % connected components in the patch/fragment graph
            regtmp.parts_compo     = cell(1,length(regtmp.parts_compo_tmp));

            for i=1:length(regtmp.parts_compo_tmp)
                regtmp.parts_compo(i) = {index(regtmp.parts_compo_tmp{i})};
            end

            [~,regtmp.I] = max(cellfun(@length, regtmp.parts_compo));
            regtmp.grp_nums = regtmp.parts_compo{regtmp.I};
        
        
            index = regtmp.grp_nums;% find([test.dup_flag]==1);
            z_tmp = (test.patch_adj(index,index)>3).*(test.patch_adj(index,index));
            z     = (z_tmp>0);
            
            [L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = graphProperties(z+z');            
            MC               = maximalCliques(z+z');
            %NODES_degs       = sum(z+z');
            Edges            = nnz(z);
            Verts            = size(z,1);
            Circuit_rank     = Edges - (Verts - 1);
            Distance_allpair = allspath(z+z');
            Dia              = max(max(Distance_allpair));
            Mean_path_length = mean(Distance_allpair(:));
            Coverage         = 1;%double(length(regtmp.resi_index{1}))/double(max_res - min_res +1);
            fragment         = unique(cell2mat(MC));
            fragment_indexed = index(unique(cell2mat(MC)));
            % --- avg degree --- % 
            scaled_adj          = (z_tmp + z_tmp') .* (100/length(test.adj_mat));
            deg                 = sum(scaled_adj(fragment,fragment),1);
            [mean_deg, var_deg] = normfit(deg);
            avg_deg = mean(deg);
            % --- avg size fragments --- %
            size_frags = zeros(1,length(fragment_indexed));
            for i = 1:length(fragment_indexed)
                size_frags(i) = 1;%length(localize(fragment_indexed(i)).atoms)*100/length(regtmp.X_refine);
            end
            [mean_fragsize, var_fragsize] = normfit(size_frags);
            avg_fragsize = mean(size_frags);               
        % --------------------------------%
        j
        params.name(j,:) = protein_name;
        params.graphprops(j, :) = [L,EGlob,CClosed,ELocClosed,COpen,ELocOpen,Circuit_rank,Dia,Mean_path_length,Coverage,1,full(mean_deg),full(var_deg),full(avg_deg),mean_fragsize,var_fragsize,avg_fragsize,length(MC)];       
        if exist('structure_done','var') == 1
           if structure_done == 1
               params.solved(j) = 1;
           else
               params.solved(j) = -1;
           end
        else
            params.solved(j) = 0;
        end
    end
    
    % ------------------------------------------
    [z_r, z_c] = find(isnan(params.graphprops));
    
    no_nan_index = setdiff([1:size(params.graphprops,2)],unique(z_c));
    % -------------------------------------------
    [params.distmat, params. z_linkage, params.outperm] = findNearestParam(params, protein_name, no_nan_index,1);
%%     
    indx_solved = find(params.solved);
    
    fid_flagfilestage1=fopen(flagfilestage1,'w');
    if ~isempty(indx_solved)
        fprintf(fid_flagfilestage1, '%d',length(indx_solved));
    else
        indx_solved = -1; % so that indx_solved doesnot give index err
    end
    fclose(fid_flagfilestage1)
    
    dendrogram_op_indx = circshift(params.outperm, indx_solved(1)+1);
        
    fid = fopen(opfile,'w');
    for i = 1:length(dendrogram_op_indx)
        if params.solved(params.outperm(i)) == 0
           fprintf(fid,'%s\n',params.name(params.outperm(i),:)); 
        end
    end       
    fclose(fid);
end
