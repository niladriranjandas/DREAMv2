function writeTogexf(adj_mat,bgr_val, file_name)
%%
%
% z_bgr = round(255 .* z_bgr);
% writeTogexf(test.resi_adj_mat, z_bgr, 'resi_graph.gexf');
%%
   if nargin ~= 3
      error('Module writeTogexf: not enough inputs.'); 
   end
   
   if size(adj_mat,1) ~= size(adj_mat,2)
      error('Module writeTogexf: adj_mat should be a square matrix.'); 
   end
   if ~isuppertriangular(adj_mat)
      adj_mat = triu(adj_mat); 
   end
   
   n_vert = length(adj_mat);
   if size(bgr_val,2) ~= 3
       error('Module writeTogexf: bgr_val matrix should be (n_vertices x 3).');
   end
   if size(bgr_val,1) ~= n_vert
       error('Module writeTogexf: bgr_val matrix size incompatible with adj_mat.');
   end
%%

   format_nodes = '\t\t\t<node id="%d" label="%d">';
   format_viz   = '\t\t\t\t<viz:color b="%d" g="%d" r="%d"/>';
   format_edge  = '\t\t\t<edge id="%d" source="%d" target="%d"/>';

%%
prefix1 = '<?xml version="1.0" encoding="UTF-8"?>';
prefix2 = '<gexf xmlns="http://www.gexf.net/1.1draft"';
prefix3 = '\txmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"';
prefix4 = '\txsi:schemaLocation="http://www.gexf.net/1.1draft http://www.gexf.net/1.1draft/gexf.xsd"';
prefix5 = '\txmlns:viz="http://www.gephi.org/gexf/viz"';
prefix6 = '\tversion="1.1">';
      
%%
  fid = fopen(file_name,'w');
    if fid < 0 
       error('Module writeTogexf: couldn''t open the file');
    end
  
    % write prefix
    fprintf(fid,strcat(prefix1,'\n'));
    fprintf(fid,strcat(prefix2,'\n'));
    fprintf(fid,strcat(prefix3,'\n'));
    fprintf(fid,strcat(prefix4,'\n'));
    fprintf(fid,strcat(prefix5,'\n'));
    fprintf(fid,strcat(prefix6,'\n'));
    
    % write to file
  fprintf(fid,'\t<graph mode="static" defaultedgetype="undirected">\n');
    % write the vertices first
    fprintf(fid,'\t\t<nodes>\n');
    for i=1:n_vert
        fprintf(fid, strcat(format_nodes,'\n'),i,i);
        fprintf(fid, strcat(format_viz,'\n'),bgr_val(i,1),bgr_val(i,2),bgr_val(i,3));
        fprintf(fid,'\t\t\t</node>\n');
    end
    fprintf(fid,'\t\t</nodes>\n');
    
    % write the edges 
    fprintf(fid,'\t\t<edges>\n');
    count = 0;
    for i=1:n_vert
       i_edge = find(adj_mat(i,:));
       for  j=1:length(i_edge)
           count = count+1;
           fprintf(fid,strcat(format_edge,'\n'),count,i,i_edge(j));
       end
    end
    fprintf(fid,'\t\t</edges>\n');
    
  fprintf(fid,'\t</graph>\n');
  fprintf(fid,'</gexf>'); 
  
  fclose(fid);
end