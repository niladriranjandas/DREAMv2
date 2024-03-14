function cI = preBreakByDenseGraph(adj_mat,min_compo_size,eta_lo)
%% "breakByDenseGraph" algorithm fails if the graph is disconnected. Hence this module
%  finds the connected components in the graph before calling
%  "breakByDenseGraph" module
%    (requires matgraph package)
%  Input:       adj_mat: adjacenecy matrix of the graph
%        min_compo_size: min component size to call the breakByDenseGraph
%                        module
%                eta_lo: break graph with density from eta_lo - eta_hi (=1)
% Output:            cI: cells each element corr to a dense subgraph
%%
   if ~issymmetric(adj_mat)
%       if ~isuppertriangular(adj_mat)
%           error('Module: preBreakByDenseGraph: adjacency matrix neighter symmetric nor uppertriangular.');
%       end
%   else
%       adj_mat = triu(adj_mat,1);
      error('Module: preBreakByDenseGraph: adjacency matrix not symmetric.');
   end
   
   if min_compo_size <3
       error('Module: preBreakByDenseGraph: min component size must atleast be 3.');
   end
   
%%

   g=graph; sparse(g);
   set_matrix(g,adj_mat>0);
   
    if ~isconnected(g)
       warning('Module: preBreakByDenseGraph: adjacency graph not connected.');
       fprintf('\n Choosing the components having size more than %d only\n',min_compo_size);
       
       conc_compo  = components(g);
       parts_compo = parts(conc_compo);
    else
       cI = breakByDenseGraph(adj_mat,eta_lo);
       return;
    end
    
%% choose only the parts with size >= min_compo_size
  cI = [];
  for i=1:length(parts_compo)
     if(length(parts_compo{i}) >= min_compo_size )         
         adj_mat_tmp = adj_mat(parts_compo{i},parts_compo{i});
         cI_tmp = breakByDenseGraph(adj_mat_tmp,eta_lo);
         if ~isempty(cI_tmp)
           % map  the atoms the original index
            cI_org_index = mapToIndex(cI_tmp,parts_compo{i});
            cI = [cI,cI_org_index];
         end
     end
  end
end

function cI_mapped = mapToIndex(cI_tmp,index)
%% function to map elements of each block of the cel cI_tmp as per index
%%
   num_cells = length(cI_tmp);
   cI_mapped = cell(1,num_cells);
     for i=1:num_cells
       num_elem = length(cI_tmp{i});
       tmp= zeros(1,num_elem);
        for j=1:num_elem
            tmp(j) = index(cI_tmp{i}(j));
        end
       cI_mapped(i) = {tmp};
     end
end
   