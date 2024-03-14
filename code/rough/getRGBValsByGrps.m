% z_bgr = getRGBValsByGrps(length(test.resi_adj_mat), test.cI_resi, find([localize.method]==1))
function bgr_val = getRGBValsByGrps(n_vert, div_cells,grp_num)

bgr_val = zeros(n_vert,3);
color_grp = colormap(lines);

  for i=1:length(grp_num)
     vert_i = div_cells{grp_num(i)};
     color_no = mod(i,length(color_grp) );
     for j=1:length(vert_i)
         bgr_val(vert_i(j),:) = color_grp(color_no,:);
     end
  end

end