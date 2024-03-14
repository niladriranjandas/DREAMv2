#!/bin/csh

foreach i (*.tab)
   set a    = (`echo $i | sed -e 's/\// /g' | sed -e 's/\.tab//g'`)
   set root = $a[$#a]

   echo $root
   cat $i | sed s/"8.5f"/"6.3f"/g > temp

   FMT_GDB -in temp -out $i

end

