#!/bin/csh -f
#$Id: run_auto_cgen.csh,v 4.10 2010/11/22 15:41:35 roberto Exp roberto $
# 
#-------------------------------------------------------------------------
# Input .csh file for running CONGEN in an automatic fashion. 
# It takes samples of the output at regular intervals checking for sanity 
# and good results. In case of bad results the program is stopped and the 
# SEED changed for a new calculation. If we have good results the program 
# is allowed to continue.
#
# **** THE VALUES OF SLEEP TIMES DEPEND ON THE TYPE OF PROCESSOR, SO FIRST
# **** IT IS WORTH TO MAKE SOME TESTS.
#
# --- IMPORTANT : some (I mean a lot) improvements are needed
#
# Roberto Tejero. April 1993. CABM-Rutgers Univ.      (RTT, Apr, 1993).
#-------------------------------------------------------------------------
#
# These are the different machines to which I have access to run CONGEN
# they are quite different in speed so we take this into account in the
# sleep times.
# Current machines are:
#
#	name		type			processor
#	--------       ----------------        --------------------------
#	bioiris		old indigo,		IP12  33MHz R2000A/3000
#	silicon		GT 2D/240, 		IP7   33MHz R2000A/3000
#	iris		Indigo2			IP22 150MHz R4400
#	biopower	Ingido2/POWER		IP26  75MHz R8000
#
#	nmrlab		Power Series		IP5   16MHz R2000A/3000
#	puero-rico	Indigo2			IP22 100MHz R4000
#	lion		Ingido2/POWER		IP21  75MHz R8000
#	tiger		GT 4D/240		IP7   33MHz R2000A/3000 
#	orca		Indigo2			IP20 150MHz R4400
#
# Roberto Tejero.
#
# -- Tweaks to this file
# -- The user is supposed to tweak with these values to have the
# -- script deciding with runs are good enough to follow
# -- These values depend greately on the system being studied
# -- how many atoms it has, how many residues, how many NOEs
# -- The only way to decide which values are the best is run the
# -- script with HIGH values, like 900 or 9000 for all the cutoffs
# -- to see how they change alogn the run.
# -- Some values are given as template in the comment part. They were
# -- used for proteins like EGF,CspA, between 60-80 residues and 
# -- with about 700-900 NOEs. 
# -- Also the CHARMM potential was used (as PARAM5 in CONGEN).
# -- If one is using AMBER potential in CONGEN, the values should be
# -- higher. 
#
 @ nice_number =     19     # the bigger the less priority
 @ vdw_cutoff =    1199     # usually here is about 150 - 175 !!
 @ vdw_cutoff_2 =   599     # usually here is about 10 - 40
 @ enoe_cutoff =   9999     # usually here is under 500 - 800
 @ enoe_cutoff_2 = 2999     # usually here is under 400
#
# =================================================================
#
# NOTE: Below this line the user SHOULD not tweak with variables
# unless he/she knows what is doing.
#
# =================================================================
#
 unsetenv TOPOL
 setenv TOPOL `pwd`/        # TOPOL variable is used by my iputs to CONGEN
 @ nuevo=`date | awk -F" " '{print $4}' | awk -F: '{print $2$1$3}'`
 set dire=`pwd`/MOLECULAUSUARIO
 set file=MOLECULAUSUARIO     
#
# Determine type of system and machine and set the paths
# accordingly
#
# set maquina="`cat /etc/sys_id | sed -e 's/\..*//'`"
 set maquina=`uname -a | awk -F" " '{print $2}'`
 set sistema=`uname -a | awk -F" " '{print $1}'`
#
 switch ($sistema) 
 case IRI* ) 
    set ps_opt="-ef"
    set NICE="/bin/nice"
    breaksw
 case Linux )
 #
 # It has been updated for work under RedHat 6.0
 #
    set maquina=`uname -a | awk -F" " '{print $2}' | awk -F"." '{print $1}'`
   # set ps_opt=" axuj"
    set ps_opt="-ef"
    set NICE="/bin/nice"
    breaksw
 case * )
    echo " New system, you should modify this file "
    exit
    breaksw
 endsw
#
 if ($maquina == "bioiris") then
 	set exe_path=/usr/local/congen/bin
 else if ($maquina == "IRIS") then
 	set exe_path=/usr/local/congen/bin
 else if ($maquina == "iris") then
 	set exe_path=/usr/local/congen/bin
 else if ($maquina == "silicon") then
 	set exe_path=/usr/people/tejero/bin
 else if ($maquina == "orca") then
 	set exe_path=/usr/people/tejero/bin
 else if ($maquina == "puerto-rico") then
 	set exe_path=/nmrapps/congen
 else if ($maquina == "tiger") then
 	set exe_path=/nmrapps/congen
 else if ($maquina == "nmrlab") then
 	set exe_path=/nmrapps/congen
 else if ($maquina == "lion") then
 	set exe_path=/usr/local/congen/sgim4i6/v2/bin
 else if ($maquina == "tortola") then
 	set exe_path=/keck/congen/sgir10i6.4/v2/bin
 else if ($maquina == "biopower") then
 	set exe_path=/usr/local/congen/sgim4i6/v2/bin
 else if ($sistema == "Linux") then
        set exe_path=/usr/local/congenLinux/linux/v2/bin
 endif
#
# Checking 
#
set usage = "  Syntax:  run_auto_MOLECULAUSUARIO.csh initial_str final_str "
set example = "  Example: run_auto_MOLECULAUSUARIO.csh 1 16 "
if (($1 == "") || ($2 == "") || \
    ($1 == " ") || ($2 == " ")) then
 echo " "
 echo "  $usage "
 echo "  $example "
 echo " "
 exit
endif

# 
# -- Some Definitions needed
#
 @ i = $1                   # number of the first structure
 @ max_structures = $2      # For generating TEN structures 
# ---
 set filetemp=${file}_0.inp
 set another_attempt="n"
 set last = "n"
# ---
 @ count = 0
# ---
# banner "good luck"
 echo " "
 echo " ----------------------------------------------------- "
 echo " Working on machine: $maquina "
 echo " ----------------------------------------------------- "
 echo "                                        "
 echo \
 " ===================================================================== "
 echo "         FILE: ${dire}/${file}          "
 echo "         DATE : `date`                  "
 echo "                                        "
# ---
# --- Check if this is a new attempt
# ---
 inicio:
 if ( $another_attempt == "n" ) then
    @ attempts = 0
    @ count = 0
    if ( ${i} > 9) then
	set directorio = ${dire}${i}
    else
	set directorio = ${dire}0${i}
    endif
    if ( -e ${directorio} ) then
       echo " Directory ${directorio} exists"
       echo " Making copy to ${directorio}.OLD "
       mv ${directorio} ${directorio}.OLD
    endif
    mkdir ${directorio}
    echo " "
    echo " ------------------- File: ${file}_${i}.out ----------------- "
 endif
# ---
# --- initialize a new attempt. IT'S RIGHT!!
# ---
 set another_attempt = "n"   
 set last = "n"
 @ veces = 0
 @ dynas = 0
 @ tmpdynas = 444
 @ old_dato = 111111
 @ newrangess = 0
# ---
# --- Copy template file to directory
# ---
 cp ${filetemp} ${directorio}/${file}_${i}.inp
 cd ${directorio} 
# --- 
# ---   Random generation of seed
# ---
 if ( -e ../semillas_usadas.stored) then
    set used = \
    ( `awk -F" " '{print $8}' ../semillas_usadas.stored` )
 else
    @ used = 0
 endif
 set semilla = ${nuevo}
 set todas = ${#used}
 newranges:
 @ newrangess++
 if ($newrangess > $todas ) then
    echo " Process killed. UNABLE TO OBTAIN A GOOD NEW SEED "
    exit 1
 endif
 @ indice = 1
 if ( -e ../semillas_usadas.stored) then
     echo "`wc ../semillas_usadas.stored`" > trash
 else
     echo " 1 " > trash
 endif
 sleep 1
 if ( -e trash ) then
    # set nuevo = `nawk '{srand();print int(rand()*1000000)}' trash`
    set nuevo = `awk '{srand();print int(rand()*1000000)}' trash`
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
 echo "`date` ---> $nuevo" >> ../semillas_usadas.stored
 sort +7 -n ../semillas_usadas.stored > ../temporas
 mv ../temporas ../semillas_usadas.stored
 @ attempts++
 sed -e /_0/s/_0/_${i}/g  \
          -e /objetivo/s/objetivo/\ ${oldseed}/ \
             ${file}_${i}.inp > tmp2
 mv tmp2 ${file}_${i}.inp
 echo \
 "DATE: `date` File: ${file}_${i} ISEED: $oldseed  Attempt: $attempts" \
     >> ../iseeds.used_all
 @ dyna_max = `fgrep -i 'dyna -' ${file}_${i}.inp | wc -l `
 if ($dyna_max == "") then
    @ dyna_max = 10
    echo " - ** -- WARNING -- **: something wrong with dyna_max "
    echo " - ** -- WARNING -- **: using dyna_max = $dyna_max "
 endif
    echo -n " - Process:  "
 endif   
 ${NICE} -${nice_number} ${exe_path}/congen \
     < ${file}_${i}.inp > ${file}_${i}.out &
#
 if (($maquina == "lion") || ($maquina == "biopower")) then
  sleep 20
 else if (($maquina == "biotech") || ($maquina == "leioa")) then
  sleep 20
 else if (($maquina == "orca") || ($maquina == "iris") || \
          ($maquina == "hyper")) then
  sleep 45
 else if (($maquina == "bioiris") || ($maquina == "silicon")) then
  sleep 60
 else
  sleep 40
 endif
#
 set identi =\
 `ps $ps_opt | grep -v 'ps $ps_opt' |grep -v 'grep'|grep congen|grep $$|awk '{print $2}'`
 echo " - Date : `date` "
 echo  \
" - Attempt: $attempts, dynas: $dyna_max PID: $identi, SEED: $oldseed "
# ---
# --- Main loop for calculations the ifs are routed to here
# ---
 chequea:
 if ($last == "y") then
    @ count = 0
    goto again
 endif		
 @ veces++
#
 if (($maquina == "lion") || ($maquina == "biopower")) then
  sleep 150
 else if (($maquina == "biotech") || ($maquina == "leioa")) then
  sleep 30
 else if (($maquina == "orca") || ($maquina == "iris")) then
  sleep 200
 else if (($maquina == "bioiris") || ($maquina == "silicon")) then
  sleep 300
 else
  sleep 50
 endif
# ---
# --- obtain line numbers for the word 'DYNA -' in output file
# ---
 fgrep -in 'dyna -' ${file}_${i}.out > tmp0.$$
 set dynas = `cat tmp0.$$ | wc -l`
 if  (($dynas == $tmpdynas)  && \
      ($veces < 700) ) then
    goto chequea
 else
    set tmpdynas = $dynas
 endif
 set datocheck = `tail -1 tmp0.$$|awk -F: '{print $1}'`
# ---
# --- Looking for VdW data in the output file
# --- 
 if ( $old_dato != $datocheck ) then
    @ old_dato = $datocheck
    set dato = \
    `tail +$old_dato ${file}_${i}.out|fgrep -in 'total e'|head -1|awk -F: '{print $1}'`
    @ dato+=3
    @ dato+=$old_dato
    tail +$dato ${file}_${i}.out|head -2 > tmp.$$
    set vdw = `tail -1 tmp.$$|awk -F" " '{printf "%d", $2}'`
    set enoe = \
       `head -1 tmp.$$|awk -F" " '{printf "%d", $6}'`
    set nuevo = \
  `head -1 tmp.$$|awk -F" " '{printf "%d", ($2+$3)-($4+$5)-$6}'`
   # ---
   # --- Check VDW and if greater than cutoff stop, if not sleep a little bit 
   # --- more
   # ---
    if ( $vdw > $vdw_cutoff ) then
       set another_attempt="y"
    endif
   # ---
   # -- Check Enoe values and stop program if they are greater than cutoff
   # ---
    if (($dynas > 12) && (($enoe > $enoe_cutoff) || ($vdw > $vdw_cutoff_2))) then
          set another_attempt="y"
    endif
   # ---
   # --- Check if Enoe is going down and stop program if it is greater than
   # --- cutoff minus certain value (it's kind of arbitrary), please take
   # --- a look to the definitions of different numbers
   # ---
    if ( $dynas > 23 ) then
       if ( $enoe > $enoe_cutoff_2 ) then
          set another_attempt="y"
       endif
    endif
    echo \
 "    Line: $dato  VdW: $vdw  Enoe: $enoe  Dynamics: $dynas  Loops: $veces"
    if ( $another_attempt == "y" ) then
       kill -9 $identi
       goto continua
    else
   # ---
       if ( $veces > 700 ) then
	  @ count = 0
	  goto again
       endif 

       if ( $dynas < $dyna_max ) goto chequea
       @ count = 0
   again:
       set fin=`tail -1 ${file}_${i}.out|awk -F: '{print $1}'`
       if ( "${fin}" != "CPU TIME" ) then
	 if (($maquina == "lion") || ($maquina == "biopower")) then
	  sleep  50
	 else if (($maquina == "orca") || ($maquina == "iris")) then
	  sleep 100
	 else if (($maquina == "bioiris") || ($maquina == "silicon")) then
	  sleep 120
	 else
	  sleep 50
	 endif
       else                       
          set another_attempt="n"
          goto continua
       endif                      
       @ count+=1
       if ( $count > 999 ) then
          echo "      *** Process killed for security. WOW, Buddy! *** "
          echo "      This has been too much time working one structure! "
          echo "      *** May be it is worth to change sleep times! *** "
          kill -9 $identi
       else                       
          goto again
       endif                      
     endif                      # if another_attempt #
   else                         # if old_dato  #  
   # ---
   # --- More security, if time is bigger than 83 hours, last steps
   # ---
     if ( $veces > 999 ) then
        set last = "y"
     endif
     goto chequea
 endif                          #   if old_dato  #
# ---
# --- Check if we have more than the wanted number of structures
# --- this number is usually TEN structures if you didn't modify
# --- the value of max_structures
# ---
  continua:
# --- new attempt (the structure has failed to finish)
 if ( $another_attempt == "y" ) then
    cd ..
    goto inicio
 else
 # --- new structure (the work has been done)
    if ( $i < $max_structures ) then 
     echo " -> A good ${file}0${i} at attempt: $attempts, SEED: $oldseed "
     echo " -> A good ${file}0${i} at attempt: $attempts, SEED: $oldseed " >>\
                ../good_seeds_stored.out
     @ i+=1
     rm tmp* trash                   # cleaning up the directory
     cd ..
     goto inicio
    endif
 endif
 exit 0      # An finally THE END of the procedure and/or script
#
#   ---   G   O   O   D       L   U   C   K   !   ---    (RTT, 1994)

  ####    ####    ####   #####           #       #    #   ####   #    #
 #    #  #    #  #    #  #    #          #       #    #  #    #  #   #
 #       #    #  #    #  #    #          #       #    #  #       ####
 #  ###  #    #  #    #  #    #          #       #    #  #       #  #
 #    #  #    #  #    #  #    #          #       #    #  #    #  #   #
  ####    ####    ####   #####           ######   ####    ####   #    #


