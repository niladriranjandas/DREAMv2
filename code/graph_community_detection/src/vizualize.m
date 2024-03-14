%%visualize the clusters created
bin_size = 50; % bin size for the histogram
%map = brewermap(length(test.cI),'Set1'); 
%map = cool(length(test.cI)); % using colormap to generate colors
map = colormap(jet(length(test.cI)));
  figure
  for i=1:length(test.cI)
      hold on
      histf(test.cI{i},ceil(length(test.adj_mat)/bin_size),'facecolor',map(i,:),'facealpha',.5,'edgecolor','none')
  end
  hold on
  %axis tight
  legend;%('H1','H2','H3','H4','location','northwest')  
  legend boxoff
  hold off
  
  
  tmp=[];
 for i=1:length(test.cI)
     tmp=union(tmp,test.cI{i});
 end
 
 test.outlier = setdiff(1:length(test.adj_mat),tmp);