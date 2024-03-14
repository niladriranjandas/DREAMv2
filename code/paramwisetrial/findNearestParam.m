function [distmat, z_linkage, outperm] = findNearestParam(param, titlename, chosen_index, dispflag)
   %%
   %    req_param = string 
   %
   %         param.name  : (n*1) list of names of params
   %         param.label : (n*1) list of (1/0) labels for each param
   %    param.graphprops : graph params n * d
   %
   %        chosen_index : which param to consider , all (default)
   %%
        if nargin == 2
           chosen_index =  [1:size(param.graphprops,2)];
           dispflag = 1;           
        elseif nargin == 3
            dispflag = 1;
        end
        
%%
    work_param = param.graphprops(:,chosen_index);
    distmat = pdist2(work_param, work_param, 'euclidean');
       
    % build the dendogram
    z_linkage = linkage(work_param); % by default linkage = 1 and eucliden dist used
    
    if dispflag == 1
        figure;
        %dendrogram(z_linkage);        
        cutoff = median([z_linkage(end-2,3) z_linkage(end-1,3)]);  % from https://in.mathworks.com/help/stats/linkage.html
        [H,T,outperm] = dendrogram(z_linkage,'ColorThreshold',cutoff);
        dendrogram(z_linkage,'ColorThreshold',cutoff)
        
        set(gca,'XTickLabel',param.name);
        xtickangle(45)
        title(titlename);
    end
end