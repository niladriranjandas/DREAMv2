remark   file parallhdg_new_db.pro
remark   NOTE: these parameters are DEPRECATED. Please instead use protein.par.
remark   geometric energy function parameters for distance geometry and
remark   simulated annealing.  Modified by G.M.C. to have sigma's of Param19x and parafloat

set message off echo off end

! The values of the force constants are somewhat arbitrary. The aim
! is to make ALL of them stiff compared to the weight on the NOE 
! term. The force constant for the angles is adjusted such that the
! 1-3 distance is approximately as well maintained as the 1-2 distance.
! 
! All angles around tetrahedral carbons (type CT) have been set
! to the ideal methane value, 109.5 degrees, by JK.
!

!
remark corrected proline angles to Engh & Huber, Acta Cryst., A47, 392 (1991) values.
remark Needed to define a new atom type, CP, which is equivalent to their 
remark CH2P atom type. JJK 9/22/95
!

!
remark changed to allow the tau3 bond angle (N..Ca..C) and 
remark the omega torsion angle (Ca..C..N..Ca) to be turned off
remark by introducing kdbang and kdbimpr
remark This is necessary for the angle database stuff to work.
remark JJK 8/5/96
!

                                             {* set energy constants *} 
evaluate ($kbon = 1000)  ! kcal / mol-A^2
evaluate ($kang =  500)  ! kcal / mol-rad^2
evaluate ($kchi =  500)  ! kcal / mol-rad^2
evaluate ($kback = 500)
evaluate ($kssbon = 1000)
evaluate ($kssang = 500) 
evaluate ($kpla =  500)  ! kcal / mol-rad^2
evaluate ($kdih =    0)  ! kcal / mol-rad^2

! BONDS

 bonds     H    NA                 $kbon     0.98
 bonds     H    NB                 $kbon     0.98
 bond      H    NH2                $kbon     0.98
 bond      H    NH1                $kbon     0.98
 bond      H    OH                 $kbon     0.96
 bond      H    S                  $kbon     0.96
 bond      HA   CT                 $kbon     1.08
 bond      HA   CP                 $kbon     1.08 ! CP added
 bond      HA   C                  $kbon     1.08
 bond      HC   NC2                $kbon     1.00
 bond      HC   NH1                $kbon     0.98
 bond      HC   NH3                $kbon     1.04
 bond      C    C                  $kbon     1.38
 bond      C    CT                 $kbon     1.53
 bond      C    N                  $kbon     1.305
 bond      C    NP                 $kbon     1.305
 bond      C    NR                 $kbon     1.305
 bond      C    NH1                $kbon     1.305
 bond      C    NH2                $kbon     1.305
 bond      C    NC2                $kbon     1.305
 bond      C    O                  $kbon     1.215
 bond      C    OC                 $kbon     1.22
 bond      C    OH                 $kbon     1.38
 bond      C    CA                 $kbon     1.400
 bond      CA   CA                 $kbon     1.400
 bond      CA   CB                 $kbon     1.404
 bond      CA   CN                 $kbon     1.400
 bond      CA   CT                 $kbon     1.510
 bond      CA   HA                 $kbon     1.080
 bond      CB   CN                 $kbon     1.419
 bond      CC   CT                 $kbon     1.504
 bond      CC   CV                 $kbon     1.375
 bond      CC   NA                 $kbon     1.385
 bond      CN   NA                 $kbon     1.380
 bond      CR   HA                 $kbon     1.080
 bond      CR   NA                 $kbon     1.343
 bond      CR   NB                 $kbon     1.3350
 bond      CT   CT                 $kbon     1.53
 bond      CP   CP                 $kbon     1.53
 bond      CP   CT                 $kbon     1.492 ! CP hydrogen bonds are guesses
 bond      CT   N                  $kbon     1.49
 bond      CP   N                  $kbon     1.473 ! CP
 bond      CT   NH1                $kbon     1.49
 bond      CT   NH2                $kbon     1.49
 bond      CT   NH3                $kbon     1.49
 bond      CT   NC2                $kbon     1.49
 bond      CT   OH                 $kbon     1.42
 bond      CT   S                  $kbon     1.81
 bond      CV   HA                 $kbon     1.080
 bond      CV   NB                 $kbon     1.394
 bond      CX   CB                 $kbon     1.459
 bond      CX   CT                 $kbon     1.495
 bond      CX   CW                 $kbon     1.352
 bond      CW   HA                 $kbon     1.080
 bond      CW   NA                 $kbon     1.381
 bond      S    S                  $kssbon   2.02

! ANGLES

 angle     H    NH1  H             $kang     107.5
 angle     H    NH1  C             $kang     120.0
 angle     H    NH1  CT            $kang     120.0
 angle     H    NH2  H             $kang     120.0
 angle     H    NH2  C             $kang     120.0
 angle     H    NH2  CT            $kang     120.0
 angle     H    OH   CT            $kang     108.0
 angle     H    S    CT            $kang     108.0
 angle     H    OH   C             $kang     108.0
 angle     HC   NH3  HC            $kang     109.5
 angle     HC   NH3  CT            $kang     109.5
 angle     HC   NH3  CP            $kang     109.5 ! CP, in N-terminal pro
 angle     HC   NC2  HC            $kang     120.0
 angle     HC   NC2  C             $kang     120.0
 angle     HC   NC2  CT            $kang     107.5
 angle     HA   C    C             $kang     120.0
 angle     HA   C    NH1           $kang     120.0
 angle     HA   C    NH2           $kang     120.0
 angle     HA   C    NR            $kang     120.0
 angle     HA   C    O             $kang     120.0
 angle     HA   CR   NA            $kang     124.2     !corrected by JK
 angle     HA   CR   NB            $kang     124.2     !corrected by JK
 angle     HA   CT   HA            $kang     109.5
 angle     HA   CT   C             $kang     109.5
 angle     HA   CT   CT            $kang     109.5
 angle     HA   CT   N             $kang     109.5
 angle     HA   CT   NH1           $kang     109.5
 angle     HA   CT   NH3           $kang     109.5
 angle     HA   CT   NC2           $kang     109.5
 angle     HA   CT   OH            $kang     109.5 !corrected by JK
 angle     HA   CT   S             $kang     109.5
 angle     HA   CV   NB            $kang     125.05    !corrected by JK
 angle     HA   CW   NA            $kang     125.65    !corrected JK
 angle     HA   CT   CP            $kang     109.5 ! CP hydrogen angles are guesses
 angle     HA   CP   CT            $kang     109.5 ! CP
 angle     HA   CP   HA            $kang     109.5 ! CP
 angle     HA   CP   CP            $kang     109.5 ! CP
 angle     HA   CP   N             $kang     109.5 ! CP
 angle     C    C    C             $kang     126.5
 angle     C    C    CT            $kang     122.3
 angle     C    C    NH1           $kang     108.6
 angle     C    C    NP            $kang     109.5
 angle     C    C    NR            $kang     109.5
 angle     C    C    OH            $kang     122.3
 angle     C    CA   CA            $kang     120.000
 angle     C    CA   HA            $kang     120.0
 angle     C    CT   CT            $kang     109.5
 angle     C    CT   N             $kbbang   111.8 ! updated with CP
 angle     C    CT   NH1           $kbbang   109.5 ! changed to support angle database JJK 
 angle     C    CT   NH2           $kang     109.5
 angle     C    CT   NH3           $kang     109.5
 angle     C    N    CT            $kang     122.6 ! updated with CP
 angle     C    N    CP            $kang     125.0 ! CP
 angle     C    NH1  C             $kang     119.1
 angle     C    NH1  CT            $kang     120.0  !corrected by JK unq. ARG
 angle     C    NC2  CT            $kang     107.5
 angle     C    NP   C             $kang     106.0
 angle     C    NR   C             $kang     106.0
 angle     CA   C    CA            $kang     120.0
 angle     CA   C    OH            $kang     120.0
 angle     CA   CA   CA            $kang     120.0
 angle     CA   CA   CB            $kang     120.0
 angle     CA   CA   CN            $kang     120.0
 angle     CA   CA   CT            $kang     120.0
 angle     CA   CA   HA            $kang     120.000 
 angle     CA   CB   CN            $kang     118.0
 angle     CA   CN   CB            $kang     122.0
 angle     CA   CT   CT            $kang     109.5
 angle     CA   CT   HA            $kang     109.5000
 angle     CB   CA   HA            $kang     120.0
 angle     CB   CX   CT            $kang     128.6
 angle     CB   CX   CW            $kang     106.4 
 angle     CB   CN   NA            $kang     104.40
 angle     CC   CT   CT            $kang     109.5
 angle     CC   CT   HA            $kang     109.5000
 angle     CC   CV   HA            $kang     125.05    !corrected by JK
 angle     CC   CV   NB            $kang     109.9000
 angle     CC   NA   CR            $kang     107.3000
 angle     CC   NA   H             $kang     126.300
 angle     CN   CA   HA            $kang     120.00
 angle     CN   NA   CW            $kang     111.6000
 angle     CN   NA   H             $kang     123.100
 angle     CR   NA   H             $kang     126.3000
 angle     CR   NB   H             $kang     126.3000
 angle     CV   NB   H             $kang     126.3000
 angle     CR   NB   CV            $kang     105.3000
 angle     CT   C    N             $kang     116.9 ! updated with CP
 angle     CT   C    NH1           $kang     117.5
 angle     CT   C    NH2           $kang     120.0
 angle     CT   C    O             $kang     121.25  !corrected by JK (pept planarity)
 angle     CT   C    OH            $kang     117.5
 angle     CT   C    OC            $kang     118.0
 angle     CT   CC   CV            $kang     131.900
 angle     CT   CC   NA            $kang     122.2000
 angle     CT   CT   CT            $kang     109.5 !JK
 angle     CT   CT   N             $kang     103.0 ! updated with CP
 angle     CP   CP   N             $kang     103.2 ! CP
 angle     CP   CP   CT            $kang     106.1 ! CP
 angle     CP   CT   CT            $kang     104.5 ! CP
 angle     CT   CT   NH1           $kang     109.5
 angle     CT   CT   NH2           $kang     109.5
 angle     CT   CT   NH3           $kang     109.5
 angle     CT   CT   NC2           $kang     109.5
 angle     CT   CT   OH            $kang     109.5 !corrected by JK
 angle     CT   CT   S             $kang     109.5
 angle     CT   CX   CW            $kang     125.0000
 angle     CT   N    CT            $kang     119.5 !PRO
 angle     CT   N    CP            $kang     112.0 ! CP
 angle     CT   NH3  CT            $kang     109.5 !PRO N-ter
 angle     CT   S    CT            $kang     97.20
 angle     CT   S    S             $kssang   104.2
 angle     CV   CC   NA            $kang     105.900
 angle     CW   NA   H             $kang     125.300
 angle     CX   CB   CN            $kang     108.800
 angle     CX   CT   CT            $kang     109.5  
 angle     CX   CT   HA            $kang     109.5
 angle     CX   CW   HA            $kang     125.65   !corrected by JK
 angle     CX   CW   NA            $kang     108.7
 angle     N    C    O             $kang     122.0 ! updated with CP
 angle     NA   CR   NB            $kang     111.60
 angle     NH1  C    NC2           $kang     120.0  !corrected by JK
 angle     NH1  C    NR            $kang     117.0
 angle     NH1  C    O             $kang     121.25 !corrected by JK (pept planarity)
 angle     NH2  C    O             $kang     120.6
 angle     NC2  C    NC2           $kang     120.0  !corrected by JK
 angle     O    C    OH            $kang     124.5
 angle     OC   C    OC            $kang     124.0  !corrected by JK
 angle     CA   CN   NA            $kang     133.6  !trp addition
 angle     CA   CB   CX            $kang     133.2  !trp addition

! IMPROPERS

! For dihedrals and impropers, the following convention was adopted:
! All dihedral terms maintaining planarity (esp. omega) have been
! converted into impropers. The only dihedrals left are around
! rotatable bonds. 

 improper  H    X    X    C        $kpla    0    0.0
 improper  H    X    X    NH1      $kpla    0    0.0
 improper  H    X    X    NH2      $kpla    0    0.0
 improper  H    X    X    O        $kpla    0    0.0
 improper  HC   X    X    NH1      $kpla    0    0.0
 improper  HC   X    X    NC2      $kpla    0    0.0
 improper  HC   X    X    NH3      $kpla    0    0.0
 improper  HA   X    X    C        $kpla    0    0.0
 improper  C    X    X    C        $kpla    0    0.0
 improper  C    X    X    CT       $kpla    0    0.0
 improper  C    X    X    N        $kpla    0    0.0
 improper  C    X    X    NH1      $kpla    0    0.0
 improper  C    X    X    NH2      $kpla    0    0.0
 improper  C    X    X    NC2      $kpla    0    0.0
 improper  C    X    X    NR       $kpla    0    0.0
 improper  C    X    X    O        $kpla    0    0.0
 improper  C    X    X    OH       $kpla    0    0.0
 improper  C    X    X    OC       $kpla    0    0.0
 improper  CT   X    X    N        $kpla    0    0.0
 improper  CT   X    X    NH1      $kpla    0    0.0
 improper  CT   X    X    NH3      $kpla    0    0.0
 improper  CT   X    X    NC2      $kpla    0    0.0
 improper  CT   X    X    NR       $kpla    0    0.0
 improper  NH1  X    X    NH1      $kpla    0    0.0
 improper  NH1  X    X    NR       $kpla    0    0.0
 improper  NH1  CT   NC2  HC       $kpla    0    0.0    ! Arg
 improper  CT   C    NH2  H        $kpla    0    0.0    ! Asn and Gln

! note the way the omega planarity is defined which allows switch
! between trans and cis by a patch. 

 improper  CT   C    NH1  CT       $kbbimp  0  180.0    ! peptide omega--
                                                          ! changed to support angle db JJK 
 improper  CT   C    N    CT       $kbbimp  0  180.0    ! proline
 improper  CT   N    C    CP       $kback   0  180.0    ! CP peptide bond
 improper  CT   N    CP   C        $kback   0  180.0    ! CP cis peptide bond
 improper  O    C    NH1  CT       $kback   0    0.0    ! CO peptide planarity
 improper  O    C    N    CT       $kback   0    0.0    ! proline
 improper  H    NH1  C    CT       $kback   0    0.0    ! NH peptide planarity
 improper  C    CT   CT   NH1      $kback   0  180.0    ! trans peptide bond
 improper  C    CT   NH1  CT       $kback   0  180.0    ! cis   peptide bond
 improper  C    CT   CT   N        $kback   0  180.0    ! Pro trans peptide bond
 improper  C    CT   N    CT       $kback   0  180.0    ! Pro cis   peptide bond

! 65.977 degrees for the chirality restraints are from
! isolated CAs minimized with bonds and angles only --JK

 improper  HA   NH1  C    CT       $kchi    0   65.977  ! CA chirality
 improper  HA   N    C    CT       $kchi    0   65.977  ! Pro CA chirality
 improper  HA   NH3  C    CT       $kchi    0   65.977  ! N-terminal CA chirality
 improper  HA   C    NH1  CT       $kchi    0   65.977  ! D  CA chirality
 improper  HA   C    N    CT       $kchi    0   65.977  ! D  Pro CA chirality
 improper  HA   C    NH3  CT       $kchi    0   65.977  ! D  N-terminal CA chirality
 improper  HA   CT   OH   CT       $kchi    0   65.977  ! Thr CB chirality
 improper  HA   CT   CT   CT       $kchi    0   -65.977  ! val, leu, ile CB chirality

! stereospecific hydrogen impropers
 improper  HA   HA   NH1  C        $kchi   0   -70.874  ! stereo GLY CA
 improper  HA   HA   CT   HA       $kchi   0   -66.514  ! methyl
 improper  HA   HA   S    HA       $kchi   0   -66.0    ! met methyl
 improper  HC   HC   CT   HC       $kchi   0   -66.0    ! amine for lys & N terminus
 improper  HC   HC   CT   CT       $kchi   0   -70.874  ! PRO N-terminus
 improper  HC   HC   CT   CP       $kchi   0   -70.874  ! CP, N-terminal pro
 improper  HA   HA   CT   CT       $kchi   0   -70.874  ! methylene
! need to get better values for this
 improper  HA   HA   CT   CP       $kchi   0   -70.874  ! CP methylene
 improper  HA   HA   CP   N        $kchi   0   -70.874  ! CP methylene

 improper  HA   HA   CT   CA       $kchi   0   -70.874  ! methylene phe 
 improper  HA   HA   CT   CX       $kchi   0   -70.874  ! methylene trp
 improper  HA   HA   CT   CC       $kchi   0   -70.874  ! methylene his
 improper  HA   HA   CT   C        $kchi   0   -70.874  ! methylene glu, asp
 improper  HA   HA   C    NH1      $kchi   0   -70.874  ! methylene gly
 improper  HA   HA   CT   NH1      $kchi   0   -70.874  ! methylene arg
 improper  HA   HA   CT   NH3      $kchi   0   -70.874  ! methylene lys
 improper  HA   HA   CT   N        $kchi   0   -70.874  ! methylene pro
 improper  HA   HA   CT   S        $kchi   0   -73.230  ! methylene cys, met
 improper  HA   HA   CT   OH       $kchi   0   -71.884  ! methylene ser


! hold rings flat (from parallhgeo)
! trp
 improper   CB   CN   CA   CA      $kpla    0    0.0
 improper   CN   CA   CA   CA      $kpla    0    0.0
 improper   CA   CA   CA   CA      $kpla    0    0.0
 improper   CA   CA   CA   CB      $kpla    0    0.0
 improper   CA   CA   CB   CN      $kpla    0    0.0
 improper   CA   CB   CN   CA      $kpla    0    0.0
 improper   CW   NA   CN   CA      $kpla    0  180.0
 improper   CW   CX   CB   CA      $kpla    0  180.0
 improper   NA   CN   CA   CA      $kpla    0  180.0
 improper   NA   CN   CB   CA      $kpla    0  180.0
 improper   CX   CB   CA   CA      $kpla    0  180.0
 improper   CX   CB   CN   CA      $kpla    0  180.0
 improper   HA   CA   CA   CA      $kpla    0  180.0
 improper   H    NA   CN   CB      $kpla    0  180.0
 improper   HA   CW   NA   CN      $kpla    0  180.0
 improper   CT   CX   CB   CN      $kpla    0  180.0
!phe 
 improper   CT   CA   CA   CA      $kpla    0  180.0
!tyr
 improper   OH   C    CA   CA      $kpla    0  180.0
 improper   HA   CA   CA   C       $kpla    0  180.0
 improper   HA   CA   C    CA      $kpla    0  180.0
 
! improper   CA   C    CA   HA      $kpla    0    0.0
 improper   CA   CA   C    CA      $kpla    0    0.0
 improper   CA   CA   C    HA      $kpla    0    0.0
 improper   CA   CA   CA   C       $kpla    0    0.0
 improper   CA   CA   CB   HA      $kpla    0    0.0
 improper   CA   CA   CN   CA      $kpla    0    0.0
 improper   CA   CB   CA   CA      $kpla    0    0.0
 improper   CA   CN   CA   HA      $kpla    0    0.0
 improper   CB   CX   CN   CA      $kpla    0    0.0
 improper   CC   CT   NA   CV      $kpla    0    0.0
 improper   CC   NA   CR   NB      $kpla    0    0.0
 improper   CN   CB   CA   NA      $kpla    0    0.0
 improper   CN   CB   CX   CW      $kpla    0    0.0
 improper   CN   NA   CW   CX      $kpla    0    0.0
 improper   CR   NA   NB   HA      $kpla    0    0.0
 improper   CR   NB   CV   CC      $kpla    0    0.0
 improper   CV   CC   NA   CR      $kpla    0    0.0
 improper   CV   NB   CC   HA      $kpla    0    0.0
 improper   CW   CX   NA   HA      $kpla    0    0.0
 improper   CW   NA   CN   CB      $kpla    0    0.0
 improper   CX   CW   CB   CT      $kpla    0    0.0
 improper   NA   CN   CB   CX      $kpla    0    0.0
 improper   NA   CR   CC   H       $kpla    0    0.0
 improper   NB   CV   CR   H       $kpla    0    0.0
 improper   NA   CR   NB   CV      $kpla    0    0.0
 improper   NA   CW   CN   H       $kpla    0    0.0
 improper   NA   CW   CX   CB      $kpla    0    0.0
 improper   NB   CV   CC   NA      $kpla    0    0.0


! DIHEDRALS

 dihedral  CA   CA   CT   CT       $kdih    3    0.0
 dihedral  CW   CX   CT   CT       $kdih    3    0.0
 dihedral  NA   CC   CT   CT       $kdih    3    0.0
 dihedral  X    NH1  CT   X        $kdih    3    0.0  ! chi1 - chi4
 dihedral  X    CT   CT   X        $kdih    3    0.0  ! chi1 - chi4
 dihedral  X    C    CT   X        $kdih    3    0.0  ! chi1 - chi4
 dihedral  X    S    CT   X        $kdih    3    0.0  ! chi1 - chi4

!Radii as in CHARMM param19x
! the radius is sigma*2^(-5/6)
!  use repel of 0.80
!                   eps     sigma       eps(1:4) sigma(1:4)  
!                                                                               
 nonbonded  C       0.0903   3.2072      0.0903   3.2072   !         
 NONBonded  CA      0.120    3.2072      0.120    3.2072
 NONBonded  CB      0.145    3.2072      0.145    3.2072
 NONBonded  CC      0.145    3.2072      0.145    3.2072
 NONBonded  CN      0.145    3.2072      0.145    3.2072
 NONBonded  CR      0.1200   3.2072      0.1200   3.2072
 nonbonded  CT      0.0903   3.2072      0.0903   3.2072
 nonbonded  CP      0.0903   3.2072      0.0903   3.2072 ! CP
 NONBonded  CV      0.1200   3.2072      0.1200   3.2072
 NONBonded  CW      0.1200   3.2072      0.1200   3.2072
 NONBonded  CX      0.1450   3.2072      0.1450   3.2072
 nonbonded  H       0.0498   1.4254      0.0498   1.4254   !         
 nonbonded  HA      0.0045   2.6157      0.0045   2.6157
 nonbonded  HC      0.0498   1.4254      0.0498   1.4254
 nonbonded  N       0.1592   2.7618      0.1592   2.7618   !         
 NONBonded  NA      0.1592   2.7618      0.1592   2.7618
 NONBonded  NB      0.1592   2.7618      0.1592   2.7618
 nonbonded  NC2     0.1592   2.7618      0.1592   2.7618
 nonbonded  NH1     0.1592   2.7618      0.1592   2.7618
 nonbonded  NH2     0.1592   2.7618      0.1592   2.7618
 nonbonded  NH3     0.1592   2.7618      0.1592   2.7618
 nonbonded  O       0.2342   2.6406      0.2342   2.6406   !         
 nonbonded  OC      1.0244   2.6406      1.0244   2.6406
 nonbonded  OH      0.2342   2.6406      0.2342   2.6406
 nonbonded  S       0.0239   3.3854      0.0239   3.3854   !         

! the following nbfixes allow hydrogen bonding 
! the distance used is (2A/B)^(1/6)*repel                                    distances
!                         A    B   A1-4  B1-4
 nbfix      H     NB    44.2  1.0  44.2  1.0               !         2.111   1.900   1.689   1.583
 nbfix      H     O     44.2  1.0  44.2  1.0 
 nbfix      H     OC    44.2  1.0  44.2  1.0 
 nbfix      H     OH    44.2  1.0  44.2  1.0 
 nbfix      HC    NB    44.2  1.0  44.2  1.0 
 nbfix      HC    O     44.2  1.0  44.2  1.0 
 nbfix      HC    OC    44.2  1.0  44.2  1.0 
 nbfix      HC    OH    44.2  1.0  44.2  1.0 

set message on echo on end
