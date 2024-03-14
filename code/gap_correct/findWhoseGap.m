function [modelgaps, gaps_fill_struct] = findWhoseGap(gapmat, distmat,gaps)
%%
%
%%

	[r1,c1] = size(gapmat);
	[r2,c2] = size(distmat);
	
	if (r1 ~= r2) || (c1 ~= c2)
		error('size mismatch between gapmat and distmat')
	end
	
	modelgaps = zeros(c1,r1);
%%
        for i = 1:c1
        	tmp_col = gapmat(:,i);
        	pos1 = find(tmp_col == 1);
        	if isempty(pos1)
        		pos0 = find(tmp_col == 0);
        		if ~isempty(pos0)
        		   dists_to_boundary = distmat(pos0,i);
        		   [m,I] = max(dists_to_boundary); %check if this can be sorted in decreasing order
        		   modelgaps(i,:) = pos0(I);
        		else
        		   pos_minus1 = find(tmp_col == -1);        		   
        		   dists_to_allowed = distmat(pos_minus1,i);
        		   [m,I] = min(dists_to_allowed);  % check if this can be sorted in increasing order
        		   modelgaps(i,:) = pos_minus1(I);
                end
            else
			    modelgaps(i,1:length(pos1)) = pos1;	
            end
        end	
%% 
  models = unique(modelgaps(:,1));
  for i=1:length(models)
     pos = find(modelgaps(:,1) == models(i));
     disp([models(i),-1,pos']);
  end
  disp('---');
  for i=1:length(models)
     gaps_fill_struct(i).models = models(i);
     pos = find(modelgaps(:,1) == models(i));
     gaps_fill_struct(i).resi   = [gaps(pos).resi]; 
  end
  
  
  
  
  
  
  
  