#!/bin/csh

foreach i (*.tab)
   set a    = (`echo $i | sed -e 's/\// /g' | sed -e 's/\.tab//g'`)
   set root = $a[$#a]

   echo $root
   cat /home/shenyang/projects/SPARTA2/src/ANN_tab/hbond_s2_EF/tab/errorSurface/header $i > a
   mv a $i
end

