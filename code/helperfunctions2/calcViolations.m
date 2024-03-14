function [mean_err, min_err, max_err, err_percent] = calcViolations(X_coord, bounds, ind)
%% function to calculate the violations
%  Input:    X_coord: 3 * n
%             bounds: (i,j,d_ij) sometimes 4th col indicates weights
%                ind: 1 - upper bound
%                     2 - lower bound
%                     3 - equality bound
% (Note that indices i and j in bounds array are indexed as per X_coord
%  i.e. i and j refers to i^th and j^th extries in X_coord)
%  Output:         mean: mean of the violations               
%               min_err: min violation
%               max_err: max violation
%           err_percent: % of entries violated
%%
     err_tol = 1e-4;
    
     if nargin ~= 3
         error('module: calcViolations: incorrect no. of inputs.');
     end
     
     if size(X_coord,1) ~= 3
         error('module: calcViolations: X_coord is 3 times n.');
     end
     
     [num_bounds, cols] = size(bounds);
     if cols < 3
         error('module: calcViolations: bounds matrix size error.');         
     end
     
     min_ind = min( min(bounds(:,1)), min(bounds(:,2)) );
     max_ind = max( max(bounds(:,1)), max(bounds(:,2)) );
     
     if max_ind > size(X_coord,2)
         error('module: calcViolations: index not proper in the bounds matrix.');
     end
     
     if min_ind ~= 1
         warning('module: calcViolations: lowest index of bounds matrix is not 1.');
     end
     
     if ~ismember(ind, [1,2,3])
         error('module: calcViolations: ind must be 1,2 or 3.');
     end
%%
  err_sum = 0;
  min_err   = 1e3;   max_err = 0;
  count_err = 0;
      for i=1:num_bounds       % proceeding each entry wise instead of whole array to save space
        dist_ij = norm(X_coord(:,bounds(i,1)) - X_coord(:,bounds(i,2)));
        diff    = bounds(i,3) - dist_ij;
        switch(ind)
            case 1
                if diff < 0                     
                    err_sum = err_sum + abs(diff);
                    %count_err = count_err +1;
                    count_err = count_err + ( abs(diff) < err_tol); %err count only if error < err_tol
                end
            case 2
                if diff > 0
                    err_sum = err_sum + diff;
                    %count_err = count_err +1;
                    count_err = count_err + (diff < err_tol); %err count only if error < err_tol
                end
            case 3
                if diff
                    err_sum = err_sum + abs(diff);
                    %count_err = count_err +1;
                    count_err = count_err + ( abs(diff) < err_tol); %err count only if error < err_tol
                end
        end
        if abs(diff) > max_err
           max_err = abs(diff); 
        end
        if abs(diff) < min_err
           min_err = abs(diff); 
        end
      end
  
      mean_err = err_sum / num_bounds;
      err_percent = (count_err*100)/num_bounds;
      
end