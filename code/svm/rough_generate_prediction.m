% see do_prediction.m for generating param structure

x = [param.L;param.EGlob;param.CClosed;param.ELocClosed;param.COpen;param.ELocOpen;param.Circuit_rank;param.Dia;param.Mean_path_length;param.mean_deg;param.var_deg;param.avg_deg;param.mean_fragsize;param.var_fragsize;param.avg_fragsize;param.Maximalcliques]';
y_ = [param.label]';
cutoff = 0.31;
y  = double(y_ <= cutoff);
%% 
no_nan_index=[1:4,7:16];
new_x = x(:,no_nan_index);
%% get rank from svmrfe
PARAM.kerType=0; [ftRank_linear,ftScore_linear] = ftSel_SVMRFECBR(x(:,no_nan_index),y,PARAM);
PARAM.kerType=2; [ftRank_rbf,ftScore_rbf]       = ftSel_SVMRFECBR(x(:,no_nan_index),y,PARAM);
ftRank_linear_indexed = no_nan_index(ftRank_linear);
ftRank_rbf_indexed    = no_nan_index(ftRank_rbf);
%% train svm
stop_index = 10;
chosen_x_linear = x(:,ftRank_linear_indexed(1:stop_index));
chosen_x_rbf    = x(:,ftRank_rbf_indexed(1:stop_index));
c1_sigmoid_chosenx = fitcsvm(chosen_x_rbf,y,'KernelFunction','mysigmoid','Standardize',true);
c1_rbf_chosenx = fitcsvm(chosen_x_rbf,y,'KernelFunction','rbf','BoxConstraint',Inf,'ClassNames',[0,1]);
c1_linear_chosenx = fitcsvm(chosen_x_linear,y,'KernelFunction','linear','BoxConstraint',Inf,'ClassNames',[0,1]);
c1_poly_chosenx = fitcsvm(chosen_x_rbf,y,'KernelFunction','polynomial','BoxConstraint',Inf,'ClassNames',[0,1]);

