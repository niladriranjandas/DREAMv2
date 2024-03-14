function [anclocrefarr, anclocref] = doMultiAnclocRefine(wh_eq_cons_all_org, wh_up_bounds_org, wh_lo_bounds_org, fillQ, wh_Comp,path)
%%
% alias matlab='LD_PRELOAD=/usr/lib64/libstdc++.so.6.0.22 /usr/local/bin/matlab'
% OR
% setenv('LD_LIBRARY_PATH', '/usr/lib64/libstdc++.so.6.0.22') before
% calling the function
%% extract the  bounds
 wh_eq_cons_all = extractCons(wh_eq_cons_all_org,fillQ.atom_map);
 wh_lo_bounds   = extractCons(wh_lo_bounds_org, fillQ.atom_map);
 wh_up_bounds   = extractCons(wh_up_bounds_org, fillQ.atom_map);
%%
 flag=1;
 i=1;
 
 tol=-1*(1e-4);
 maxitter = 10;%5;
 
 name = fillQ.filename(1:strfind(fillQ.filename,'.pdb')-1);
 %path = '/home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/rough/fragmentwisereport/density_diff/1xxe_rough';
%%% path = '/home/niladri/Documents/Disco_etc_all_in_1/our_algo/code/multi_refine';
 anclocrefarr(i).X_refine = fillQ.x_ref_chiral;
 [anclocrefarr(i).eq_mean_err, anclocrefarr(i).eq_min_err, anclocrefarr(i).max_err, anclocrefarr(i).eq_err_percent] = calcViolations(fillQ.x_ref_chiral,wh_eq_cons_all,3);        
 [anclocrefarr(i).up_mean_err, ~, anclocrefarr(i).up_max_err, anclocrefarr(i).up_err_percent] = calcViolations(fillQ.x_ref_chiral,wh_up_bounds,1);
 [anclocrefarr(i).lo_mean_err, anclocrefarr(i).lo_min_err, ~, anclocrefarr(i).lo_err_percent] = calcViolations(fillQ.x_ref_chiral,wh_lo_bounds,2);
 anclocrefarr(i).name_i = fillQ.file_full;
 cmd = strcat('./percentSecondary.sh',{' '},anclocrefarr(i).name_i);
 [anclocrefarr(i).shellstatus, shellout] = system(cmd{1});
 anclocrefarr(i).shellout = splitlines(shellout);
 
 while flag
      i = i+1
      [anclocrefarr(i).X_refine,anclocrefarr(i).info] =  postProcessingMe(anclocrefarr(i-1).X_refine,...
                                                    wh_eq_cons_all, wh_lo_bounds, full(wh_up_bounds),...
                                                   [1e2,1e2,1e2,-0.001],[10,10,10,10,10]);

      [anclocrefarr(i).eq_mean_err, anclocrefarr(i).eq_min_err, anclocrefarr(i).eq_max_err, anclocrefarr(i).eq_err_percent] = calcViolations(anclocrefarr(i).X_refine,wh_eq_cons_all ,3);
      fprintf('\n Equality bounds: Mean error: %d \t Min error: %d \t Max error: %d \t Violation percent: %f', anclocrefarr(i).eq_mean_err, anclocrefarr(i).eq_min_err, anclocrefarr(i).eq_max_err, anclocrefarr(i).eq_err_percent);
        
      [anclocrefarr(i).up_mean_err, ~, anclocrefarr(i).up_max_err, anclocrefarr(i).up_err_percent] = calcViolations(anclocrefarr(i).X_refine,wh_up_bounds,1);
      fprintf('\n Upper bounds:    Mean error: %d \t Max error: %d \t Violation percent: %f', anclocrefarr(i).up_mean_err, anclocrefarr(i).up_max_err, anclocrefarr(i).up_err_percent);

      [anclocrefarr(i).lo_mean_err, anclocrefarr(i).lo_min_err, ~, anclocrefarr(i).lo_err_percent] = calcViolations(anclocrefarr(i).X_refine,wh_lo_bounds,2);
      fprintf('\n Lower bounds:    Mean error: %d \t Min error: %d \t Violation percent: %f', anclocrefarr(i).lo_mean_err, anclocrefarr(i).lo_min_err, anclocrefarr(i).lo_err_percent);
      
      eq_mean_diff = anclocrefarr(i).eq_mean_err - anclocrefarr(i-1).eq_mean_err;
      eq_max_err   = anclocrefarr(i).eq_max_err  - anclocrefarr(i-1).eq_max_err;
      
      up_mean_diff = anclocrefarr(i).up_mean_err - anclocrefarr(i-1).up_mean_err;
      up_max_err   = anclocrefarr(i).up_max_err  - anclocrefarr(i-1).up_max_err;
      
      anclocrefarr(i).name_i = strcat(path,filesep,name,'_',string(i),'.pdb');
      {writeToPDB(anclocrefarr(i).name_i,fillQ.atom_map, anclocrefarr(i).X_refine,wh_Comp)};
      
      cmd = strcat('./percentSecondary.sh',{' '},anclocrefarr(i).name_i);
      [anclocrefarr(i).shellstatus, shellout] = system(cmd{1});
      anclocrefarr(i).shellout = splitlines(shellout);
      % {writeToPDB('1xxe_ancloc_0001_refine.pdb',fillQ.atom_map, anclocref.X_refine,wh_Comp)}
      
     if (eq_max_err < tol) flag1=1; else flag1=0; end
     if (up_max_err < tol) flag2=1; else flag2=0; end
      
      flag = (flag1 || flag2) && (i<maxitter);            
      %flag = (i<maxitter);
 end
 
 anclocref.X_refine = anclocrefarr(i-1).X_refine;
 anclocref.info     = anclocrefarr(i-1).info;
 
anclocref.eq_mean_err    = anclocrefarr(i).eq_mean_err;
anclocref.eq_min_err     = anclocrefarr(i).eq_min_err ;
anclocref.eq_max_err     = anclocrefarr(i).eq_max_err ;
anclocref.eq_err_percent = anclocrefarr(i).eq_err_percent;

anclocref.up_mean_err    = anclocrefarr(i).up_mean_err;
anclocref.up_max_err     = anclocrefarr(i).up_max_err ;
anclocref.up_err_percent = anclocrefarr(i).up_err_percent;

anclocref.lo_mean_err    = anclocrefarr(i).lo_mean_err ;
anclocref.lo_min_err     = anclocrefarr(i).lo_min_err ;
anclocref.lo_err_percent = anclocrefarr(i).lo_err_percent;

end
  

