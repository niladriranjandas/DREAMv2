omit.percent = 10;
omit.lo_indx = randi(length(up_bounds),1, round((omit.percent/100)*length(up_bounds)));

up_bounds_bkp = up_bounds;
omit.include_indx = setdiff(1:length(up_bounds),omit.lo_indx);
up_bounds = up_bounds(omit.include_indx,:);