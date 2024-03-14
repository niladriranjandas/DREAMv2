#!/bin/csh -f
#$Id: run_fromdisman_cgen.csh,v 4.10 2010/11/22 15:41:35 roberto Exp roberto $
#
# *** Last modified in July 1997
#	Now (JUL97) it is able to generate random seeds for the dynamics.
#	Though this is not really needed as we start from different
#	coordinates from a previous calc (we HAVE PDB coordinates)
# 	I liked the idea of having more freedom for the DYNAMICS using
#	also different seed number for the random generator.
#	* Also added the path for running in TORTOLA.CABM.RUTGERS.EDU
#	it is a R10000 SGI processor.
#
#  (JUL1995)
#  Now it takes three arguments, 
#  1. The name of the file where the PDB coordinates are stored as a 
#     multiple model 
#  2. The number of the starting structure 
#  3. The number of final structure
#
#
# changes to do in each case
#----------
@ salida = 0
unsetenv TOPOL
setenv TOPOL `pwd`/
#
# Things to modify 
#
set usage = "  Syntax:  $0 FILE.pdb initial_str final_str "
set example = "  Example: $0 zdomain.pdb 1 16 "
if (($1 == "") || ($2 == "") || ($3 == "") || \
    ($1 == " ") || ($2 == " ") || ($3 == " ")) then
 echo " "
 echo "  $usage "
 echo "  $example "
 echo " "
 exit
endif
#
set file_pdb = $1			#name of file with multiple PDB coords
set work_file = "MOLECULAUSUARIO_sa_fromdisman" #name for work file
set dir_base = "MOLECULAUSUARIO"		#base DIRECTORY to store data
set dir_base0 = ${dir_base}0		# i.e. ZDOMCHK01 , ZDOMCHK02 , ...
set chain_id = "MOLE_ID"		#one letter ID for the chain
#
set  from = $2         # starting structure
set  to = $3	       # ending structures
@ nice_number = 19
#
# for timing purposes
#
set inicio = `date`
set inicio_dias = `date | awk -F" " '{print $3}'`
set inicio_hors = `date | awk -F" " '{print $4}' |awk -F":" '{print $1}'`
set inicio_mins = `date | awk -F" " '{print $4}' |awk -F":" '{print $2}'`
set inicio_secs = `date | awk -F" " '{print $4}' |awk -F":" '{print $3}'`

#
# Preliminary checking
#
echo " "
if !(-f $file_pdb) then
	echo " Family file: '$file_pdb' doesn't exist. Please check\! "
 	@ salida = 1
endif
if !(-f ${work_file}_0.inp ) then
	echo " Work file: '${work_file}_0.inp' doesn't exist. Please check\! "
 	@ salida = 1
endif
set identi = `fgrep -i "atom" $file_pdb | fgrep  " $chain_id  " | wc -l`
if !( $identi > 0) then
	echo " "
	echo " ** It seems as if the identification for the chain "
	echo " ** doesn't agree with the id in the PDB file. Please check\!"
        echo " "
	@ salida = 1
endif
if ( $salida > 0 ) then
	echo " "
	echo " ** Something went wrong. Please check\!"
	echo " "
	exit (-1)
else
	echo " ** Checks passed. Continuing now ... "
	echo " "
endif
#
# Some checking
#
if ($from > $to) then
  @ this_script = $from - $to
  set direccion = "backward"
else
  @ this_script = $to - $from
  set direccion = "forward"
endif
@ this_script+=1
@ number_str = `fgrep -i model\ \ \   $file_pdb | wc -l`
#
 set maquina = "`cat /etc/sys_id | sed -e 's/\..*//'`"
 if ($maquina == "bioiris") then
 	set exe_cgen_path = "/usr/local/congen/bin"
	set exe_path = "/usr/people/tejero/bin"
 else if ($maquina == "IRIS") then
 	set exe_cgen_path = "/usr/local/congen/bin"
	set exe_path = "/usr/people/tejero/bin"
 else if ($maquina == "silicon") then
 	set exe_cgen_path = "/usr/people/tejero/bin"
	set exe_path = "/usr/people/tejero/bin"
 else if ($maquina == "orca") then
 	set exe_cgen_path = "/usr/people/tejero/bin"
	set exe_path = "/usr/people/tejero/bin"
 else if ($maquina == "puerto-rico") then
 	set exe_cgen_path = "/tmp2d/nmrapps/congen"
	set exe_path = "/nmrusr/people/tejero/bin"
 else if ($maquina == "tiger") then
 	set exe_cgen_path = "/tmp2d/nmrapps/congen"
	set exe_path = "/nmrusr/people/tejero/bin"
 else if ($maquina == "nmrlab") then
 	set exe_cgen_path = "/tmp2d/nmrapps/congen"
	set exe_path = "/nmrusr/people/tejero/bin"
 else if ($maquina == "lion") then
 	set exe_cgen_path = "/usr/local/congen/sgim4i6/v2/bin"
 else if ($maquina == "tortola") then
 	set exe_cgen_path = "/keck/congen/sgir10i6.4/v2/bin"
	set exe_path = "/nmrusr/local/bin"
 else if ($maquina == "biopower") then
 	set exe_cgen_path = "/usr/local/congen/sgim4i6/v2/bin"
	set exe_path = "/biothetas/local/bin"
 else if ($maquina == "biolin") then
 	set exe_cgen_path = "/usr/local/congen/sgim2i6/v2/bin"
	set exe_path = "/biothetas/local/bin"
 endif
echo " =========================================================== "
echo "   **** Working CONGEN from previuos calc. for "
echo "       PDB file: $file_pdb               "
echo "     input file: ${work_file}_0.inp      "
echo "    directories: ${dir_base}${from} to ${dir_base}${to}"
echo " ----------------------------------------------------------- "
echo " "
echo " Log : "
#exit
#
 if ($number_str != $this_script) then
	echo "  "
	echo " ---------------------  IMPORTANT ---------------------- "
	echo " The number of structures in file:  $file_pdb don't  "
	echo " agree with the number you gave to me in this script "
	echo " file, ** ----------- Please, check ----------------- ** "
	echo " "
	echo " Data in: $file_pdb --> $number_str "
	echo " Data in this script --> $this_script  (from $from to $to) "
	echo "  "
	echo " I assume that you are working in a SUBSET of the "
	echo " total number and continue working ... "
	echo "  "
 else
	echo "  "
	echo " Working in family: $dir_base from file: $file_pdb "
	echo " from structure: $from to structure: $to "
	echo "  "
 endif
#
#
# nuevo is used to set a new calculation seed number
#
@ nuevo=`date | awk -F" " '{print $4}' | awk -F: '{print $2$1$3}'`
#
if ($direccion == "backward" ) then
   while ($from >= $to)
else
   while ($from <= $to)
endif
	#
	# This has been added to allow also for different seeds in running from
	# previous calcs. Before we always used the same SEED as we start from
	# a previous *DIFFERENT* structure but I decided to use also different 
	# SEEDS to allow more freedom to move in the DYNAMIC stage.
	#
	# RTT. March 1997. Valencia University.
	#
	# --- 
	# ---   Random generation of seed
	# ---
    if ( -e semillas_usadas.stored ) then
	set used = \
	( `/usr/bin/awk -F" " '{print $8}' semillas_usadas.stored` )
    else
	@ used = 0 
    endif
    set semilla = ${nuevo}
    set todas = ${#used}
    newranges:
	if ( -e semillas_usadas.stored ) then
		@ newrangess++
	else
		@ newrangess = 0
        endif
	if ($newrangess > $todas ) then
	    echo " Process killed. UNABLE TO OBTAIN A GOOD NEW SEED "
	    exit 1
	endif
	@ indice = 1
	if ( -e semillas_usadas.stored ) then
	    echo "`wc semillas_usadas.stored`" > trash
	else
	    echo " 1 " > trash
	endif
	/bin/sleep 1
	if ( -e trash ) then
	    set nuevo = `/usr/bin/nawk '{srand();print int(rand*1000000)}' trash`
	else
	    echo " file trash doesn't exist, using seed in file: $nuevo "
	endif
	    # ---
	    # ---  Check if the new seed has been used before
	    # ---
	while ($indice <= $todas ) 
	    if ($nuevo >= $used[$indice]) then
		if ($nuevo == $used[$indice]) then
		    set semilla = $indice
		    goto newranges
		endif
	    endif
	    @ indice++
	end
	set oldseed = $nuevo
	if ( $oldseed == "" ) goto newranges
	echo "`date` ---> $nuevo" >> semillas_usadas.stored
	sort +7 -n semillas_usadas.stored > temporas
	mv temporas semillas_usadas.stored
# END OF THE SEED STUFF
#
    if ($from<10) set base="$dir_base0"
#
    if ($from>9) set base="$dir_base"
#        
	echo "$file_pdb" > tmp.$$_$from
	echo "$chain_id" >> tmp.$$_$from
	if ($2 != $3) then
		echo "$from" >> tmp.$$_$from
	endif
	echo "${work_file}_${from}.pdb">>tmp.$$_$from
	echo "n" >> tmp.$$_$from
	${exe_path}/ortocg < tmp.$$_$from >& log_build

        rm tmp.$$_$from
#
	mkdir $base$from
        mv ${work_file}_${from}.pdb $base$from
	sed -e /objetivo/s/objetivo/${work_file}_${from}.pdb/g \
	    -e /semillatmp/s/semillatmp/${nuevo}/g \
               ${work_file}_0.inp \
             > temporal 
        sed -e /_0/s/_0/_${from}/g \
	       temporal > $base$from/${work_file}_$from.inp
        rm temporal
	cd $base$from
	/bin/nice -${nice_number} ${exe_cgen_path}/congen \
               < ${work_file}_$from.inp \
                 > ${work_file}_$from.out
	cd ..
        echo " ... Done with $base$from ... "	
#
	if ($direccion == "backward") then
		@ from-=1	
        else
		@ from+=1
	endif
#
end
rm log_build

set final = `date`
set final_dias = `date | awk -F" " '{print $3}'`
set final_hors = `date | awk -F" " '{print $4}' |awk -F":" '{print $1}'`
set final_mins = `date | awk -F" " '{print $4}' |awk -F":" '{print $2}'`
set final_secs = `date | awk -F" " '{print $4}' |awk -F":" '{print $3}'`
echo " "
echo " Starting: $inicio "
echo " Ending  : $final "
echo " "
set tdias = `echo "${final_dias}-${inicio_dias}" | bc`
set thors = `echo "${final_hors}-${inicio_hors}" | bc`
set tmins = `echo "${final_mins}-${inicio_mins}" | bc`
set tsecs = `echo "${final_secs}-${inicio_secs}" | bc`
if ($tsecs < 0) then
	 @ tmins--
	 set tsecs = `echo "60+${tsecs}"|bc`
endif
if ($tmins < 0) then
	@ thors--
	 set tmins = `echo "60+${tmins}"|bc`
endif
if ($thors < 0) then
	@ tdias--
	 set thors = `echo "60+${thors}"|bc`
endif
#
set ttotal = \
 `echo "scale=2; (${tdias}*24*60+${thors}*60+${tmins})*60+${tsecs}"|bc`
echo " Tiempo usado: ${tdias} d. ${thors} h. ${tmins} m. ${tsecs} s. "
echo " Tiempo total: $ttotal secs. "
echo " "
echo "  .... All done buddy\!  .... "
echo "    "
exit
