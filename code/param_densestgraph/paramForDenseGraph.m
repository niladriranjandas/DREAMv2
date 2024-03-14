function [param_test_global,param_test_local] = paramForDenseGraph(num_runs,start_num_vertices,end_num_vertices)

%% This module creates a graphs of size varying from start_num_vertices to  end_num_vertices.
%  For each of these graph generated, spanning tree is found and then one
%  by one edges are added and tested for global rigidity (using gotler
%  rigidity test). The value of (num_edges/num_edges in a complete graph)
%  is returned for the first instant at which the global_rigit param is 1 
%  in the gotler rigidity test.
%
%  Input: num_runs : no. of times experiment is repeated for a graph with
%                    certain no. of vertices
%         start_num_vertices : size of graph to begin exp with.
%         end_num_vertices   : size of graph to end exp with.
%
%  Output:param_test: 2 * (start_num_vertices - end_num_vertices +1) matrix
%         num_vertices | param value at which graph with num_vertices
%                        is globaly rigid (avg of exp repeated num_runs of
%                        times.
%                      
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% check the inputs
if nargin ~= 3
    error('Module paramGorDenseGraph: Not enough inputs');
end

if num_runs < 1
    error('Module paramGorDenseGraph: num_runs has to be atleas 1');
end

if (start_num_vertices < 4 || end_num_vertices < 4)
    error('Module paramGorDenseGraph: the number of vertices must be atleast 4');
end

if (start_num_vertices > end_num_vertices)
    error('Module paramGorDenseGraph: the graph size wrongly mentioned');
end

%% call the required modules
param_test_global = zeros(end_num_vertices - start_num_vertices + 1,2);
param_test_local  = zeros(end_num_vertices - start_num_vertices + 1,2);
count = 1;
for i = start_num_vertices:end_num_vertices
 fprintf('\n-------Vertex : %d-----',i)    
     tmp_param_g = zeros(1,num_runs);
     tmp_param_l = zeros(1,num_runs);
     
     for j=1:num_runs
         tmp_param_g(j) = testGlobalRigityAddingEdges(i);
         tmp_param_l(j) = testLocalRigityAddingEdges(i);
     end
 %%excluding zeros (if present) from mean calculation
   [~,ii_g,val_g] = find(tmp_param_g);
   avg_res_g = sum(val_g) / length(ii_g);
   
   [~,ii_l,val_l] = find(tmp_param_l);
   avg_res_l = sum(val_l) / length(ii_l);
   
  param_test_global(count,:) = [i,avg_res_g];
  param_test_local(count,:) = [i,avg_res_l];
  count = count + 1;
end