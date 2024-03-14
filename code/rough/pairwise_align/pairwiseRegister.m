function regXnIndx = pairwiseRegiser(grpnum_seq, x_n_indx)
%%
%   Input : grpnum_seq : 1*k
%           x_n_indx: (same as doRegForGroup)
%   Output: regXnIndx.X: dim * n     regXnIndx.indx: 1*n
%%
   if( length(grpnum_seq) ~= length(x_n_indx))
       error('seq of aligments and x_n_indx struct donot match');
   end
%%    
    target.X    = x_n_indx(grpnum_seq(1)).x;
    target.indx = x_n_indx(grpnum_seq(1)).ind;
    for i=2:length(grpnum_seq)
        mobile.X    = x_n_indx(grpnum_seq(i)).x;
        mobile.indx = x_n_indx(grpnum_seq(i)).ind;
        
        target = findPairwiseTransform(target, mobile);
    end        
    
    regXnIndx = target;
end