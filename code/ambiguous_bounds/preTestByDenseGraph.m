function [cI_test, ret_densities] = preTestByDenseGraph(adj_mat,min_compo_size,eta_lo)
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
       [cI_test, ret_densities] = testByDenseGraph(adj_mat,eta_lo);
       return;
    end
    
%% choose only the parts with size >= min_compo_size
  cI_test_ = [];
  ret_densities_ = [];
  ret_densities  = [];
  for i=1:length(parts_compo)
     cI_tmp_indxed = [] ; % changed 16-07-22  while testing for 1pbu_ambi5r1_noe 
     ret_densities_tmp = []; % changed 16-07-22  while testing for 1pbu_ambi5r1_noe 
     if(length(parts_compo{i}) >= min_compo_size )         
         adj_mat_tmp = adj_mat(parts_compo{i},parts_compo{i});
         %cI_tmp = testByDenseGraph(adj_mat_tmp,eta_lo);
         [cI_tmp, ret_densities_tmp] = testByDenseGraph(adj_mat_tmp,eta_lo);  % changed 16-07-22  while testing for 1pbu_ambi5r1_noe   
         for k=1:length(cI_tmp)
            if ~isempty(cI_tmp{k})
            % map  the atoms the original index
                %[cI_tmp_indxed{k}, ret_densities_tmp{k}] = mapToIndex(cI_tmp{k},parts_compo{i});  % changed 16-07-22  while testing for 1pbu_ambi5r1_noe 
                cI_tmp_indxed{k} = mapToIndex(cI_tmp{k},parts_compo{i});   % changed 16-07-22  while testing for 1pbu_ambi5r1_noe  
            else
                cI_tmp_indxed{k} = cI_tmp{k};
            end
         end
     end
     cI_test_{i} = cI_tmp_indxed;
     ret_densities_{i} = ret_densities_tmp;
  end
  
  % concatenate cI_test and ret_densities into 1, otherwise
  % getSubgraphDensity creates error
  cI_test = cat(1, cI_test_{:});
  for k=1:length(ret_densities_)
     ret_densities = [ret_densities, ret_densities_{k}]; 
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
   