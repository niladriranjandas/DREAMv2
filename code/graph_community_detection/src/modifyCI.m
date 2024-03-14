function [mod_cI,range_cI] = modifyCI(cI)

%% dense subgraph algo creates clusters, such that points inside are in scattered clumps
%  This function separates those clumps
%%
 
 jump=100;
 
   count = 0;
   for i=1:length(cI)       
       %%detect peaks
       consec_diff = cI{i}(2:end) - cI{i}(1:end-1);
       [~,col,~] = find(consec_diff>jump);
       
       count = count+1;
       mod_cI{count} = cI{i}(1:col(1));
       range_cI{count} = cI{i}(1):cI{i}(col(1));
       %cI{i}(1:col(1))
       for j=2:length(col)
           count=count+1;
           mod_cI{count} = cI{i}(col(j-1)+1:col(j));
           range_cI{count} = cI{i}(col(j-1)+1):cI{i}(col(j));
           %cI{i}(col(j-1)+1:col(j))
       end
       count = count+1;
       mod_cI{count} = cI{i}(col(end)+1:end);
       range_cI{count} = cI{i}(col(end)+1):cI{i}(end);
   end

end