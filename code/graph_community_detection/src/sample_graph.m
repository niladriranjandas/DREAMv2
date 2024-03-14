adj_mat = zeros(18,18);

adj_mat(1,2)=1;
adj_mat(1,6)=1;
adj_mat(1,7)=1;
adj_mat(1,5)=1;
adj_mat(1,3)=1;

adj_mat(2,3)=1;
adj_mat(2, 4 )=1;
adj_mat(2, 7 )=1;
adj_mat(2,10 )=1;
             
adj_mat(3, 4 )=1;
adj_mat(3, 5 )=1;
adj_mat(3, 7 )=1;
adj_mat(3,10 )=1;
adj_mat(3,13 )=1;

adj_mat(4, 5 )=1;
adj_mat(4, 7 )=1;
            
adj_mat(5, 6 )=1;
adj_mat(5, 7 )=1;
adj_mat(5,16 )=1;
            
adj_mat(6, 7 )=1;
adj_mat(6,17 )=1;
            
adj_mat(8,10 )=1;
            
adj_mat(9,10 )=1;
            
adj_mat(11,13)=1;
            
adj_mat(12,13)=1;
            
adj_mat(14,15)=1;
adj_mat(14,16)=1;
            
adj_mat(15,16)=1;

adj_mat(17,18)=1;

label_vertex = [1:size(adj_mat,1)]';

[row,col] = find(adj_mat);