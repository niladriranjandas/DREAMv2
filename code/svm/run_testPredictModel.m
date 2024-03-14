cd ..
addpath(genpath(pwd));
cd code;

%%
pred153 = load('predictmodel153.mat');

test_x = [pred153.param.L;pred153.param.EGlob;pred153.param.CClosed;pred153.param.ELocClosed;pred153.param.COpen;pred153.param.ELocOpen;pred153.param.Circuit_rank;pred153.param.Dia;pred153.param.Mean_path_length;pred153.param.mean_deg;pred153.param.var_deg;pred153.param.avg_deg;pred153.param.mean_fragsize;pred153.param.var_fragsize;pred153.param.avg_fragsize;pred153.param.Maximalcliques]';

test_y_ = [pred153.param.label]';
cutoff  = 0.31
test_y  = double(test_y_ <= cutoff);

test_chosen_x_linear = test_x(:,ftRank_linear_indexed(1:stop_index));
test_chosen_x_rbf    = test_x(:,ftRank_rbf_indexed(1:stop_index));


test_accuracy_rbf = testPredictModel(test_chosen_x_rbf, test_y, c1_rbf_chosenx, 'random', 50, 'itter', 50);
mean(test_accuracy_rbf)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

>> test_accuracy_rbf = testPredictModel(test_chosen_x_rbf, test_y, c1_rbf_chosenx, 'random', 50, 'itter', 50);
mean(test_accuracy_rbf)

ans =

    0.1728

>> test_accuracy_rbf = testPredictModel(test_chosen_x_rbf, test_y, c1_rbf_chosenx, 'random', 50, 'itter', 1000);
>> mean(test_accuracy_rbf)

ans =

    0.1807

>> test_accuracy_rbf = testPredictModel(test_chosen_x_rbf, test_y, c1_rbf_chosenx, 'random', 50, 'itter', 100000);
>> mean(test_accuracy_rbf)

ans =

    0.1798
    
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:10
	test_accuracy_rbf = testPredictModel(test_chosen_x_rbf, test_y, c1_rbf_chosenx, 'random', 35, 'itter', 100000); 
	mean(test_accuracy_rbf)    
end

    0.1796
    0.1796
    0.1795
    0.1794
    0.1795
    0.1793
    0.1797
    0.1795
    0.1798
    0.1797



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%	Training set %%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
pred96 = load('predictmodel96.mat');

test_x = [pred96.param.L;pred96.param.EGlob;pred96.param.CClosed;pred96.param.ELocClosed;pred96.param.COpen;pred96.param.ELocOpen;pred96.param.Circuit_rank;pred96.param.Dia;pred96.param.Mean_path_length;pred96.param.mean_deg;pred96.param.var_deg;pred96.param.avg_deg;pred96.param.mean_fragsize;pred96.param.var_fragsize;pred96.param.avg_fragsize;pred96.param.Maximalcliques]';

test_y_ = [pred96.param.label]';
cutoff  = 0.31
test_y  = double(test_y_ <= cutoff);

test_chosen_x_linear = test_x(:,ftRank_linear_indexed(1:stop_index));
test_chosen_x_rbf    = test_x(:,ftRank_rbf_indexed(1:stop_index));

%%

for i=1:10
	test_accuracy_rbf = testPredictModel(test_chosen_x_rbf, test_y, c1_rbf_chosenx, 'random', 35, 'itter', 100000); 
	mean(test_accuracy_rbf)    
end


    0.0895
    0.0895
    0.0895
    0.0892
    0.0894
    0.0895
    0.0892
    0.0895
    0.0895
    0.0892



