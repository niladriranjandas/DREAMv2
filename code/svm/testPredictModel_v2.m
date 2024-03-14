function [accuracy, confuse_mat] = testPredictModel_v2(Xtest, Ytest, mdl, varargin)
%%
%	accuray: real no. [0,1]: number of failed an accepted cases
%	
%	Xtest: rows: no. of data points
%	       cols: no. of variables
%	Ytest: rows: no. of data points (= no. of rows in Xtest)	
%	       cols: 1 (labels) 
%	mdl:	classifier learned from ftcsvm in MATLAB
%	NameValue pairs:
%		'random_pts': integer (choose pts from Xtest,Ytest to test predictions)
%		'itter'	    : no. of iterations to choose 'random'
%%

	[r_Xtest, c_Xtest] = size(Xtest);
	[r_Ytest, c_Ytest] = size(Ytest);
	
	if r_Xtest ~= r_Ytest
		error('No. of rows (data points) in data (Xtest) and labels (Ytest) dont match');
	end	
	if c_Ytest ~= 1
		error('No. of rows labels (Ytest) should be 1');
	end		
	
	for i = 1:2:length(varargin)
		if ischar(varargin{i})
			params.(varargin{i}) = varargin{i+1};
		end
	end
	
	accuracy = zeros(1,params.itter);
	confuse_mat = zeros(2,2,params.itter);
%%
	for i=1:params.itter
		rand_rows  = randperm(r_Xtest, params.random);
		rand_Xtest = Xtest(rand_rows,:);
		rand_Ytest = Ytest(rand_rows);

		% find no. of 1's (positive) and 0's (negative)
		num_1 = sum(rand_Ytest)
		num_0 = params.random - num_1
		
		pred_Ytest = predict(mdl, rand_Xtest);

		% calculate the individual cases 00,01,10,11
		for k=1:length(pred_Ytest)
           if (rand_Ytest(k)== 0) && (pred_Ytest(k) == 0)
               confuse_mat(1,1,i) = confuse_mat(1,1,i) + 1;
           elseif (rand_Ytest(k)== 0) && (pred_Ytest(k) == 1)  
               confuse_mat(1,2,i) = confuse_mat(1,2,i) + 1;
           elseif (rand_Ytest(k)==1) && (pred_Ytest(k) == 0)
               confuse_mat(2,1,i) = confuse_mat(2,1,i) + 1;
           elseif (rand_Ytest(k)==1) && (pred_Ytest(k) == 1)
               confuse_mat(2,2,i) = confuse_mat(2,2,i) + 1;
           end
		end
		
		confuse_mat(1,1,i) = confuse_mat(1,1,i)/num_0; confuse_mat(1,2,i) = confuse_mat(1,2,i)/num_0; 
		confuse_mat(2,1,i) = confuse_mat(2,1,i)/num_1; confuse_mat(2,2,i) = confuse_mat(2,2,i)/num_1; 
		accuracy(i) = sum(abs(pred_Ytest - rand_Ytest)) / params.random;
	end
	
	
	
