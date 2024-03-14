function transform=findPairwiseTransform(targetStruct, mobileStruct)
%% function to find pairwise transformation
% Input:   targetStruct.x: dim * n1  targetStruct.indx: 1*n1 (global index)
%          mobileStruct.x: dim * n2  mobileStruct.indx: 1*n2 (global index)
% Output:  transform.X: dim * n   transform.indx: 1*n (global index)
%%
    X     = targetStruct.X;
    Xindx = targetStruct.indx;
    
    Y     = mobileStruct.X;
    Yindx = mobileStruct.indx;
    
%% validate
    if size(X,2) ~= length(Xindx)
        error('Target no. of points and index missmatch')
    end
    if size(Y,2) ~= length(Yindx)
        error('Mobile no. of points and index missmatch')
    end
%%
    ComIndx = intersect(Xindx,Yindx);
    XlocalComIndx = findLocal(ComIndx, Xindx);
    YlocalComIndx = findLocal(ComIndx, Yindx);
    
    target = X(:,XlocalComIndx);
    mobile = Y(:,YlocalComIndx);
    [~,Z,rigid] = procrustes(target',mobile','scaling',false);
    %[R,t] = procrustes(target,mobile);
    
    t = rigid.c(1,:)';
    Y_new = rigid.T * Y + repmat(t,1,length(Y));
    
    
%% populate transform.X in increasing order of indices of union(Xindx,Yindx)
    I_XnotY = setdiff(Xindx,Yindx);
    local_inX = findLocal(I_XnotY, Xindx);
    
    transform.indx = sort(union(I_XnotY,Yindx));
    
    localIndx_intransform_X = findLocal(I_XnotY, transform.indx);
    localIndx_intransform_Y = findLocal(Yindx, transform.indx);
    
    transform.X = nan(size(X,1),length(transform.indx));
    transform.X(:,localIndx_intransform_X) = X(:,local_inX);
    transform.X(:,localIndx_intransform_Y) = Y_new;
end

function localIndx = findLocal(findindx, fromIndx)
    localIndx = nan(1,length(findindx));
    for i=1:length(findindx)
       localIndx(i) = find(fromIndx==findindx(i));
    end
end