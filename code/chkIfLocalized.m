%% 
% small code which checks for successful localizations
%%

index_chk = find([localize.method]==1);

for chk_i = 1:length(index_chk)
   [chk_r,~] = size(localize(index_chk(chk_i)).info.obj );
   if chk_r == 1
       localize(index_chk(chk_i)).method = -3;
   end
end