%load('svm/predict_svmrfe_new_12.mat')
load('svm/predictmodel96.mat')

index = find([test.dup_flag]==1);
z_tmp = (test.patch_adj(index,index)>3).*(test.patch_adj(index,index));
z = (z_tmp>0);
[L,EGlob,CClosed,ELocClosed,COpen,ELocOpen] = graphProperties(z+z');
MC = maximalCliques(z+z');
NODES_degs = sum(z+z');
Edges = nnz(z);
Verts = size(z,1);
Circuit_rank = Edges - (Verts - 1);
Distance_allpair = allspath(z+z');
Dia = max(max(Distance_allpair));
Mean_path_length = mean(Distance_allpair(:));
Coverage = 1;%double(length(reg.resi_index{1}))/double(max_res - min_res +1);
fragment         = unique(cell2mat(MC));
fragment_indexed = index(unique(cell2mat(MC)));
	% --- avg degree --- %
scaled_adj = (z_tmp + z_tmp') .* (100/length(test.adj_mat));
deg = sum(scaled_adj(fragment,fragment),1);
[mean_deg, var_deg] = normfit(deg);
avg_deg = mean(deg);
	% --- avg size fragments --- %
size_frags = zeros(1,length(fragment_indexed));
for i = 1:length(fragment_indexed)
size_frags(i) =1;% length(localize(fragment_indexed(i)).atoms)*100/length(test.adj_mat); % correct this
end
[mean_fragsize, var_fragsize] = normfit(size_frags);
avg_fragsize = mean(size_frags);
%x_predict=[L,EGlob,CClosed,ELocClosed,COpen,ELocOpen,Circuit_rank,Dia,Mean_path_length,Coverage,1,full(mean_deg),full(var_deg),full(avg_deg),mean_fragsize,var_fragsize,avg_fragsize,length(MC)];
x_predict=[L,EGlob,CClosed,ELocClosed,COpen,ELocOpen,Circuit_rank,Dia,Mean_path_length,full(mean_deg),full(var_deg),full(avg_deg),mean_fragsize,var_fragsize,avg_fragsize,length(MC)];
x_predict_tmp_lin =x_predict(:,ftRank_linear_indexed(1:stop_index));
x_predict_tmp_rbf =x_predict(:,ftRank_rbf_indexed(1:stop_index));
   % ------ kernel svms -------%

[predict1.predict_row_sigmoid  ,predict1.predict_col_sigmoid] = predict(c1_sigmoid_chosenx, x_predict_tmp_rbf);
[predict1.predict_row_rbf      ,predict1.predict_col_rbf]     = predict(c1_rbf_chosenx, x_predict_tmp_rbf);
[predict1.predict_row_poly     ,predict1.predict_col_poly]    = predict(c1_poly_chosenx, x_predict_tmp_rbf);
[predict1.predict_row_linear   ,predict1.predict_col_linear]  = predict(c1_linear_chosenx, x_predict_tmp_lin);

consolidate_flag = predict1.predict_row_sigmoid | predict1.predict_row_rbf | predict1.predict_row_poly | predict1.predict_row_linear;



