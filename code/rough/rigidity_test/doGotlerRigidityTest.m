function gotler_rigidity = doGotlerRigidityTest(localize)

indx = find([localize.method]==1);

is_globally_rigid = zeros(1,length(indx));
is_locally_rigid  = zeros(1,length(indx));

for i=1:length(indx)
    i
    num_nodes = length(localize(indx(i)).atoms);    
    
    eq_cons = localize(indx(i)).eq_cons_all;
    up_cons = localize(indx(i)).up_bounds;
    lo_cons = localize(indx(i)).sdp_lo_bounds;
    
    adj_mat_ = bsxfun(@or, sparse(eq_cons(:,1), eq_cons(:,2),1,num_nodes,num_nodes),...
                          sparse(up_cons(:,1), up_cons(:,2),1,num_nodes,num_nodes));
    adj_mat  = bsxfun(@or, adj_mat_,...
                          sparse(lo_cons(:,1), lo_cons(:,2),1,num_nodes,num_nodes));    
                      
     adj_mat  = (adj_mat + adj_mat')>0;
     diag(adj_mat) = 0;
        
     d=3;
     runs=10;
     [ is_globally_rigid(i), is_locally_rigid(i), GL_dof, WW ] = ... 
                            SVD_glob_rigidity_many_times_G(adj_mat,d,runs);                          
end

  gotler_rigidity = [is_globally_rigid;is_locally_rigid];
end