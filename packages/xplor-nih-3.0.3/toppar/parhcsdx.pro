REMARK Parameter file including bond and angle parameters
REMARK derived from Cambridge Data Base model structures
REMARK (R. A. Engh and R. Huber, Acta Cryst. Sect. A., 1991).
REMARK Dihedral, improper, and non-bonded parameters taken
REMARK from param19x (XPLOR--Axel T. Brunger, Yale University,
REMARK BRUNGER@YALEVMS) and assigned to new atom types
REMARK where appropriate.

!
! Please cite the following reference when using these parameters:
! Engh, R.A. and  Huber, R. (1991). Accurate Bond and 
!  Angle Parameters for X-ray Protein-Structure Refinement, 
!  Acta Cryst. A47, 392-400.
!
!

set echo=false end

REMARK These parameters use additional atom types (14 more than in
REMARK param19x) and requires the corresponding topology file for
REMARK psf file generation.

! The XPL19X comments are the original XPLOR param19x parameters.  Under
! these lines are the corresponding new parameters.  Where new atom
! types are defined, several bond or angle types correspond to a single
! XPL19X definition.  The numerical values in comments represent the
! CSD determined sample standard deviation from which the force constants
! were derived by assuming the standard deviation represents expected
! deviations from potential minima at a temperature of 293 K.  Although
! some forces have up to seven significant figures, they might be considered
! useful only to about one significant figure.

! Some values were not determined from the Cambridge Data Base (CSD) but were
! taken either from related geometry types or from X-PLOR parameters.  These
! are marked with an asterisk and include hydrogen related parameters and some
! other rarely occurring geometries.

! Only bond and angle parameters were taken from the CSD.  Dihedrals and
! improper dihedrals retain the original XPL19X target values and force
! constants.  One result is that the bonds and angles are weighted
! approximately 3 and 7 times as heavily, respectively, compared to
! XPL19X parameters.  This is not an error, but is a matter for consideration.
! In general, the relative weighting of different types of restraints
! is not standard, and depends primarily on the desired deviations from
! target values.  We think that the relative weighting of bonds and angles
! is correctly represented in these parameters.  The other parameters
! may be scaled according to any special criteria.

! IF any parameters are missing or inconsistent, contact me by E-MAIL
! at engh@nmrvex.biochem.mpg.de or engh@mpib-martinsried.mpg.dbp.de.  If there
! is no reply, please use physical, tangible, slow, reliable mail to R. Engh,
! Max-Planck-Institut fuer Biochemie, D8033 Martinsried bei Muenchen,
! Fed. Rep. of Germany.

  !XPL19X BOND C    C      450.0  1.38
bond C5W  CW   1827.161  1.433 ! 0.018
bond CW   CW   2048.443  1.409 ! 0.017

  !XPL19X BOND C    CH1E   405.0  1.52
bond C    CH1E 1342.404  1.525 ! 0.021

  !XPL19X BOND C    CH2E   405.0  1.52
bond C5   CH2E 3020.408  1.497 ! 0.014
bond C5W  CH2E  616.025  1.498 ! 0.031
bond CF   CH2E 1119.093  1.502 ! 0.023
bond CY   CH2E 1223.141  1.512 ! 0.022
bond C    CH2E  947.200  1.516 ! 0.025
bond CN   CH2E 1639.889  1.503 ! 0.019
bond C    CH2G 1827.161  1.516 ! 0.010

  !XPL19X BOND C    CH3E   405.0  1.52

  !XPL19X BOND C    CR1E   450.0  1.38
bond C5W  CR1E  947.200  1.365 ! 0.025
bond CW   CR1E 2312.500  1.398 ! 0.016
bond CW   CR1W 1342.404  1.394 ! 0.021
bond CF   CR1E 1342.404  1.384 ! 0.021
bond CY   CR1E 1342.404  1.389 ! 0.021
bond CY2  CR1E 1027.778  1.378 ! 0.024
bond C5   CR1H 4892.562  1.354 ! 0.011
bond C5   CR1E 4892.562  1.356 ! 0.011

  !XPL19X BOND C    CT     405.0  1.53

  !XPL19X BOND C    N      471.0  1.33
bond C    N    2312.500  1.341 ! 0.016

  !XPL19X BOND C    NC2    400.0  1.33
bond C    NC2  1827.161  1.326 ! 0.018

  !XPL19X BOND C    NH1    471.0  1.33
bond C5   NH1  4892.562  1.378 ! 0.011
bond CW   NH1  4892.562  1.370 ! 0.011
bond C    NH1  3020.408  1.329 ! 0.014

  !XPL19X BOND C    NH2    471.0  1.33
bond C    NH2  1342.404  1.328 ! 0.021

  !XPL19X BOND C    NP     471.0  1.33

  !XPL19X BOND C    NR     471.0  1.33
bond C5   NR   2048.443  1.371 ! 0.017

  !XPL19X BOND C    O      580.0  1.23
bond C    O    1480.000  1.231 ! 0.020
bond CN   O    1119.093  1.208 ! 0.023

  !XPL19X BOND C    OC     580.0  1.23
bond C    OC   1639.889  1.249 ! 0.019

  !XPL19X BOND C    OH1    450.0  1.38
bond CY2  OH1  1342.404  1.376 ! 0.021
bond C    OH1  1223.141  1.304 ! 0.022

  !XPL19X BOND CH1E CH1E   225.0  1.53
bond CH1E CH1E  812.071  1.540 ! 0.027

  !XPL19X BOND CH1E CH2E   225.0  1.52
bond CH1E CH2E 1480.000  1.530 ! 0.020

  !XPL19X BOND CH1E CH3E   225.0  1.52
bond CH1E CH3E  543.618  1.521 ! 0.033

  !XPL19X BOND CH1E N      422.0  1.45
bond CH1E N    2631.111  1.466 ! 0.015

  !XPL19X BOND CH1E NH1    422.0  1.45
bond CH1E NH1  1639.889  1.458 ! 0.019

  !XPL19X BOND CH1E NH2    422.0  1.45

  !XPL19X BOND CH1E NH3    422.0  1.45
bond CH1E NH3  1342.404  1.491 ! 0.021

  !XPL19X BOND CH1E OH1    400.0  1.42
bond CH1E OH1  2312.500  1.433 ! 0.016

  !XPL19X BOND CH2E CH2E   225.0  1.52
bond CH2E CH2E  657.778  1.520 ! 0.030
bond CH2P CH2E  236.800  1.492 ! 0.050
bond CH2P CH2P  512.111  1.503 ! 0.034

  !XPL19X BOND CH2E CH3E   225.0  1.54
bond CH2E CH3E  389.218  1.513 ! 0.039

  !XPL19X BOND CH2E CR1E   250.0  1.45

  !XPL19X BOND CH2E N      422.0  1.45
bond CH2P N    3020.408  1.473 ! 0.014
bond CH2P NH3  3020.408  1.473 ! 0.014  !added 5/19/93 ATB

  !XPL19X BOND CH2E NH1    422.0  1.45
bond CH2G NH1  2312.500  1.451 ! 0.016
bond CH2E NH1  1827.161  1.460 ! 0.018
bond CH3E NH1  1827.161  1.460 ! copied for special residue

  !XPL19X BOND CH2E NH2    422.0  1.45

  !XPL19X BOND CH2E NH3    422.0  1.45
bond CH2E NH3   657.778  1.489 ! 0.030
bond CH2G NH3   657.778  1.489 ! 0.030*

  !XPL19X BOND CH2E OH1    400.0  1.42
bond CH2E OH1  1480.000  1.417 ! 0.020

  !XPL19X BOND CH2E S      450.0  1.81
bond CH2E S    1480.000  1.822 ! 0.020
bond CH2E SM    512.111  1.803 ! 0.034

  !XPL19X BOND CH2E SH1E   450.0  1.81
bond CH2E SH1E  543.618  1.808 ! 0.033

  !XPL19X BOND CH3E NH1    422.0  1.49

  !XPL19X BOND CH3E NR     422.0  1.49

  !XPL19X BOND CH3E S      450.0  1.77
bond CH3E SM    170.066  1.791 ! 0.059

  !XPL19X BOND CM   OM    1115.0  1.128

  !XPL19X BOND CR1E CR1E   450.0  1.38
bond CR1E CR1E  657.778  1.382 ! 0.030
bond CR1E CR1W  947.200  1.400 ! 0.025
bond CR1W CR1W 1639.889  1.368 ! 0.019

  !XPL19X BOND CR1E NH1    450.0  1.305
bond CR1E NH1  1342.404  1.374 ! 0.021
bond CRH  NH1  1480.000  1.345 ! 0.020
bond CRHH NH1  5920.000  1.321 ! 0.010
bond CR1H NH1  4892.562  1.374 ! 0.011

  !XPL19X BOND CR1E NR     450.0  1.305
bond CR1E NR    657.778  1.382 ! 0.030

  !XPL19X BOND OC   S      400.0  1.43

  !XPL19X BOND OM   OM     600.0  1.23

  !XPL19X BOND S    S      500.0  2.02
bond S    S    2410.110  2.030 ! 0.008

  !XPL19X BOND CR1E NR     450.0  1.305
bond CRH  NR   3502.958  1.319 ! 0.013

bond H    NH1    405.0  0.98 !*
bond H    NH2    405.0  0.98 !*
bond H    OH1    450.0  0.96 !*
bond HA   C      350.0  1.08 !*
bond HA   CT     300.0  1.08 !*
bond HC   NC2    405.0  1.00 !*
bond HC   NH1    405.0  0.98 !*
bond HC   NH3    405.0  1.04 !*




  !XPL19X ANGLE C    C    C       70.0 106.5
angle C5W  CW   CW   1349.600  107.200 !1.2

  !XPL19X ANGLE C    C    CH2E    65.0 126.5
angle CW   C5W  CH2E  991.543  126.800 !1.4

  !XPL19X ANGLE C    C    CR1E    70.0 122.5
angle C5W  CW   CR1E 1943.424  133.900 !1.0
angle CW   CW   CR1E 1943.424  118.800 !1.0
angle CW   CW   CR1W 1943.424  122.400 !1.0
angle CW   C5W  CR1E  759.150  106.300 !1.6

  !XPL19X ANGLE C    C    NH1  65.0 109.0
angle CW   CW   NH1  1149.956  107.400 !1.3

  !XPL19X ANGLE CH1E C    N       20.0 117.5
angle CH1E C    N     863.744  116.900 !1.5

  !XPL19X ANGLE CH1E C    NH1     20.0 117.5
angle CH1E C    NH1   485.856  116.200 !2.0

  !XPL19X ANGLE CH1E C    O       85.0 121.5
angle CH1E C    O     672.465  120.800 !1.7

  !XPL19X ANGLE CH1E C    OC      85.0 117.5
angle CH1E C    OC    310.948  117.000 !2.5

  !XPL19X ANGLE CH2E C    CR1E    70.0 121.5
angle CH2E C5   CR1E 1149.956  129.100 !1.3
angle CH2E C5   CR1H 1149.956  131.200 !1.3
angle CH2E CF   CR1E  672.465  120.700 !1.7
angle CH2E C5W  CR1E  863.744  126.900 !1.5
angle CH2E CY   CR1E  863.744  120.800 !1.5

  !XPL19X ANGLE CH2E C    N       20.0 117.5
angle CH2E C    N     440.686  118.200 !2.1
angle CH2G C    N     440.686  118.200 !2.1 *

  !XPL19X ANGLE CH2E C    NH1     20.0 117.5
angle CH2E C5   NH1   863.744  122.700 !1.5
angle CH2E C    NH1   440.686  116.500 !2.1
angle CH2G C    NH1   440.686  116.400 !2.1 *

  !XPL19X ANGLE CH2E C    NH2     20.0 117.5
angle CH2E C    NH2   863.744  116.400 !1.5

  !XPL19X ANGLE CH2E C    NR      60.0 116.0
angle CH2E C5   NR    863.744  121.600 !1.5

  !XPL19X ANGLE CH2E C    O       85.0 121.6
angle CH2E C    O     485.856  120.800 !2.0
angle CH2G C    O     440.686  120.800 !2.1

  !XPL19X ANGLE CH2E C    OC      85.0 118.5
angle CH2E C    OC    367.377  118.400 !2.3
angle CH2G C    OC    367.377  118.400 !2.3 *

  !XPL19X ANGLE CR1E C    CR1E    65.0 120.5
angle CR1E CY2  CR1E  485.856  120.300 !2.0
angle CR1E CY   CR1E  863.744  118.100 !1.5
angle CR1E CF   CR1E  863.744  118.600 !1.5 

  !XPL19X ANGLE CR1E C    NH1     65.0 110.5
angle CR1W CW   NH1   863.744  130.100 !1.5
angle CR1E C5   NH1  1943.424  105.200 !1.0
angle CR1H C5   NH1  1943.424  106.100 !1.0

  !XPL19X ANGLE CR1E C    NP      65.0 122.5
  !XPL19X ANGLE CR1E C    NR      65.0 122.5
angle CR1E C5  NR    1943.424  109.200  !0.7 HISE, taken as 1.0

  !XPL19X ANGLE CR1E C    OH1     65.0 119.0
angle CR1E CY2  OH1   215.936  119.900 !3.0

  !XPL19X ANGLE N    C    O       85.0 121.0
angle N    C    O     991.543  122.000 !1.4

  !XPL19X ANGLE NC2  C    NC2     70.0 120.0
angle NC2  C    NC2   599.823  119.700 !1.8

  !XPL19X ANGLE NC2  C    NH1     70.0 120.0
angle NC2  C    NH1   538.345  120.000 !1.9

  !XPL19X ANGLE NH1  C    O       65.0 121.0
angle NH1  C    O     759.150  123.000 !1.6

  !XPL19X ANGLE NH2  C    O       65.0 121.0
angle NH2  C    O    1943.424  122.600 !1.0

  !XPL19X ANGLE OC   C    OC      85.0 122.5
angle OC   C    OC    337.400  122.900 !2.4

  !XPL19X ANGLE C    CH1E CH1E    70.0 110.0
angle C    CH1E CH1E  401.534  109.100 !2.2

  !XPL19X ANGLE C    CH1E CH2E    70.0 109.5
angle C    CH1E CH2E  538.345  110.100 !1.9

  !XPL19X ANGLE C    CH1E CH3E    70.0 106.5
angle C    CH1E CH3E  863.744  110.500 !1.5

  !XPL19X ANGLE C    CH1E N       45.0 111.6
angle C    CH1E N     310.948  111.800 !2.5

  !XPL19X ANGLE C    CH1E NH1     45.0 111.6
angle C    CH1E NH1   247.886  111.200 !2.8
angle C    CH1E NH3   247.886  111.200 !2.8 *


  !XPL19X ANGLE CH1E CH1E CH2E    45.0 112.5
angle CH1E CH1E CH2E  672.465  110.400 !1.7

  !XPL19X ANGLE CH1E CH1E CH3E    45.0 111.0
angle CH1E CH1E CH3E  672.465  110.500 !1.7

  !XPL19X ANGLE CH1E CH1E NH1     50.0 110.0
angle CH1E CH1E NH1   672.465  111.500 !1.7
angle CH1E CH1E NH3   672.465  111.500 !1.7 *

  !XPL19X ANGLE CH1E CH1E OH1     50.0 104.5
angle CH1E CH1E OH1   863.744  109.600 !1.5

  !XPL19X ANGLE CH2E CH1E CH3E    50.0 111.5
angle CH2E CH1E CH3E  215.936  110.700 !3.0

  !XPL19X ANGLE CH2E CH1E N       65.0 104.0
angle CH2E CH1E N    1606.136  103.000 !1.1

  !XPL19X ANGLE CH2E CH1E NH1     65.0 110.0
angle CH2E CH1E NH1   672.465  110.500 !1.7
angle CH2E CH1E NH3   672.465  110.500 !1.7 *

  !XPL19X ANGLE CH3E CH1E CH3E    50.0 111.0
angle CH3E CH1E CH3E  401.534  110.800 !2.2

  !XPL19X ANGLE CH3E CH1E NH1     65.0 108.5
angle CH3E CH1E NH1   863.744  110.400 !1.5

  !ANGLE CH3E CH1E NH2     65.0 109.5

  !ANGLE CH3E CH1E NH3     65.0 109.5
angle CH3E CH1E NH3   672.465  110.500 !1.7 *

  !XPL19X ANGLE CH3E CH1E OH1     60.0 110.5
angle CH3E CH1E OH1   485.856  109.300 !2.0

  !XPL19X ANGLE C    CH2E CH1E    70.0 112.5
angle C    CH2E CH1E 1943.424  112.600 !1.0
angle C5   CH2E CH1E 1943.424  113.800 !1.0
angle CF   CH2E CH1E 1943.424  113.800 !1.0
angle C5W  CH2E CH1E  538.345  113.600 !1.9
angle CY   CH2E CH1E  599.823  113.900 !1.8

  !XPL19X ANGLE C    CH2E CH2E    70.0 113.0
angle C    CH2E CH2E  672.465  112.600 !1.7

  !XPL19X ANGLE C    CH2E NH1     70.0 111.6
angle C    CH2G NH1   231.085  112.500 !2.9
angle C    CH2G NH3   231.085  112.500 !2.9 *

  !XPL19X ANGLE CH1E CH2E CH1E    45.0 117.0
angle CH1E CH2E CH1E  158.647  116.300 !3.5

  !XPL19X ANGLE CH1E CH2E CH2E    45.0 112.5
angle CH1E CH2E CH2P  538.345  104.500 !1.9
angle CH1E CH2E CH2E  485.856  114.100 !2.0

  !XPL19X ANGLE CH1E CH2E CH3E    45.0 113.0
angle CH1E CH2E CH3E  440.686  113.800 !2.1

  !XPL19X ANGLE CH1E CH2E OH1     45.0 111.0
angle CH1E CH2E OH1   485.856  111.100 !2.0

  !XPL19X ANGLE CH1E CH2E SH1E    50.0 112.5
angle CH1E CH2E S     367.377  114.400 !2.3
angle CH1E CH2E SH1E  367.377  114.400 !2.3 *

  !XPL19X ANGLE CH2E CH2E CH2E    45.0 110.0
angle CH2E CH2E CH2E  367.377  111.300 !2.3
angle CH2E CH2P CH2P  189.788  106.100 !3.2

  !XPL19X ANGLE CH2E CH2E N       65.0 105.0
angle CH2P CH2P N     863.744  103.200 !1.5
angle CH2P CH2P NH3    863.744  103.200 !1.5  !added ATB 5/19/93

  !XPL19X ANGLE CH2E CH2E NH1     65.0 111.0
angle CH2E CH2E NH1   401.534  112.000 !2.2

  !XPL19X ANGLE CH2E CH2E NH3     65.0 110.5
angle CH2E CH2E NH3   189.788  111.900 !3.2

  !XPL19X ANGLE CH2E CH2E S       50.0 112.5
angle CH2E CH2E SM    215.936  112.700 !3.0

  !XPL19X ANGLE C    CR1E CR1E    90.0 119.0
angle CY2  CR1E CR1E  599.823  119.600 !1.8
angle CW   CR1E CR1E 1149.956  118.600 !1.3
angle CW   CR1W CR1W 1149.956  117.500 !1.3
angle CF   CR1E CR1E  672.465  120.700 !1.7
angle CY   CR1E CR1E  863.744  121.200 !1.5

  !XPL19X ANGLE C    CR1E NH1     90.0 109.5
angle C5   CR1E NH1  1943.424  106.500 !1.0
angle C5   CR1H NH1  1943.424  107.200 !1.0
angle C5W  CR1E NH1  1149.956  110.200 !1.3

  !XPL19X ANGLE C    CR1E NR      90.0 106.5
angle C5   CR1E NR    367.377  109.500 !2.3

  !XPL19X ANGLE CR1E CR1E CR1E    90.0 120.5
angle CR1E CR1E CR1W 1149.956  121.100 !1.3
angle CR1W CR1W CR1E 1149.956  121.500 !1.3
angle CR1E CR1E CR1E  599.823  120.000 !1.8

  !XPL19X ANGLE NH1  CR1E NH1     70.0 109.0
angle NH1  CRHH NH1  1943.424  108.400 !1.0

  !XPL19X ANGLE NH1  CR1E NR      70.0 109.0
angle NH1  CRH NR   1149.956  111.700 !1.3  HisE, applied also to HisD


  !XPL19X ANGLE C    N    CH1E    80.0 120.0
angle C    N    CH1E   77.737  122.600 !5.0

  !XPL19X ANGLE C    N    CH2E    80.0 120.0
angle C    N    CH2P  115.611  125.000 !4.1
angle HC   NH3  CH2P  500.     125.000 !4.1  !added ATB 5/19/93

  !XPL19X ANGLE CH1E N    CH2E 60.0 110.0
angle CH1E N    CH2P  991.543  112.000 !1.4
angle CH1E NH3  CH2P  991.543  112.000 !1.4  ! added ATB 5/19/93

  !XPL19X ANGLE C    NH1  CH1E    77.5 120.0
angle C    NH1  CH1E  599.823  121.700 !1.8

  !XPL19X ANGLE C    NH1  CH2E    77.5 120.0
angle C    NH1  CH2G  672.465  120.600 !1.7
angle C    NH1  CH2E  863.744  124.200 !1.5
angle C    NH1  CH3E  672.465  120.600 !1.7 *

  !XPL19X ANGLE C    NH1  CR1E    60.0 108.0
angle C5   NH1  CRHH  672.465  109.300 !1.7
angle C5   NH1  CRH   672.465  109.000 !1.7 *
angle CW   NH1  CR1E  599.823  108.900 !1.8

  !XPL19X ANGLE CR1E NH1  CR1E    65.0 110.0
angle CRHH NH1  CR1H 1943.424  109.000 !1.0
angle CRH  NH1  CR1E 1149.956  106.900 !1.3

  !XPL19X ANGLE C    NR   CR1E    70.0 109.5
angle C5   NR   CRH 1943.424  105.600 !1.0 HisE

  !XPL19X ANGLE CR1E NR   CR1E    65.0 110.0
angle CR1E NR   CRH  215.936  107.000 !3.0, HisD, weak statistics

  !XPL19X ANGLE CH2E S    CH3E    50.0  99.5
angle CH2E SM   CH3E  401.534  100.900 !2.2

  !XPL19X ANGLE CH2E S    S       50.0 104.2
angle CH2E S    S     599.823  103.800 !1.8

  !XPL19X ANGLE C    C    HA      40.0 120.0

angle C    NC2  HC      35.0 120.0 !*
angle HC   NC2  HC      40.0 120.0 !*

angle C    NH1  H       30.0 120.0 !*
angle CW   NH1  H       30.0 120.0 !*
angle C5W  NH1  H       30.0 120.0 !*
angle C5   NH1  H       30.0 120.0 !*
angle CH1E NH1  H       35.0 120.0 !*
angle CRH  NH1  H       35.0 120.0 !*
angle CRHH NH1  H       35.0 120.0 !*
angle CR1H NH1  H       35.0 120.0 !*

angle CH2E NH1  H       35.0 120.0 !*
angle CH2G NH1  H       35.0 120.0 !*

angle CH3E NH1  H       35.0 120.0 !*
angle CR1E NH1  H       35.0 120.0 !*
angle C    NH2  H       30.0 120.0 !*
angle CH1E NH2  H       35.0 120.0 !*
angle CH2E NH2  H       35.0 120.0 !*
angle H    NH2  H       40.0 125.0 !*
angle CH1E NH3  HC      35.0 109.5 !*
angle CH2E NH3  HC      35.0 109.5 !*
angle CH2G NH3  HC      35.0 109.5 !*
angle HC   NH3  HC      40.0 109.5 !*
angle C    OH1  H       50.0 109.5 !*
angle CY2  OH1  H       50.0 109.5 !*
angle CH1E OH1  H       35.0 109.5 !*
angle CH2E OH1  H       35.0 109.5 !*

! SEE ABOVE REGARDING DIHEDRALS
{ Weis 5/11/92 scaled up 3x for consistency with bond and angle parameters.}
{ Note: didn't change commented lines.                                     }
dihe CH1E C    N    CH1E 1250.0    2 180.0! reduced -> to allow cis PRO
dihe CH2E C    N    CH1E 1250.0       2     180.0!  "
dihe CH2G C    N    CH1E 1250.0       2     180.0!  "
!dihe CR1E C    C    CR1E 100.0       2     180.0 ! increased
dihe CR1E CW   CW   CR1E 300.0       2     180.0 ! increased
dihe CR1E CW   CW   CR1W 300.0       2     180.0 ! increased
!dihe CR1E C    C    C   100.0       2     180.0 !
dihe CR1W CW   CW   C5W  300.0       2     180.0 !
!dihe CR1E C    C    NH1 100.0       2     180.0 !
dihe CR1E CW   CW   NH1  300.0       2     180.0 !
dihe X    C    CH1E X    0.0       3       0.0! FROM GELIN THESIS AMIDES
dihe X    C    CH2E X    0.0       3       0.0! USING A SINGLE
dihe X    C    CH2G X    0.0       3       0.0! USING A SINGLE
dihe X    C5   CH2E X    0.0       3       0.0! USING A SINGLE
dihe X    C5W  CH2E X    0.0       3       0.0! USING A SINGLE
dihe X    CF   CH2E X    0.0       3       0.0! USING A SINGLE
dihe X    CY   CH2E X    0.0       3       0.0! USING A SINGLE
dihe X    C    CR1E X   30.0       2     180.0! DIHEDRAL PER BOND RATHER
!dihe X    C    CT   X    0.0       3       0.0! THAN MULTIPLE TORSIONS.
dihe X    C    N    X   24.6       2     180.0! ALKANE TORSION REDUCED TO
dihe X    C    NC2  X   24.6       2     180.0! 1.6 FROM 1.8 TO COINCIDE WITH
dihe X    C    NH1  X  1250.0       1       0.0    ! always trans
dihe X    C    NH2  X   24.6       2     180.0
!dihe X    C    OH1  X    1.8       2     180.0
dihe X    CY2  OH1  X    5.4       2     180.0
dihe X    CH1E CH1E X    4.8       3       0.0
dihe X    CH1E CH2E X    4.8       3       0.0
dihe X    CH1E N    X    0.9       3       0.0! FROM HAGLER ET AL TABULATION OF
dihe X    CH1E NH1  X    0.9       3       0.0! EXP. DATA AND 6 31G CALC.
dihe X    CH1E NH2  X    5.4       3       0.0! PROTONATED SECONDARY AMINE
dihe X    CH1E NH3  X    1.8       3       0.0! 1/PROTON SO 3 FOR THE BOND
dihe X    CH1E OH1  X    1.5       3       0.0! CHANGED TO ROUGHLY MEOH
dihe X    CH2E CH2E X    4.8       3       0.0
dihe X    CH2E CH2P X    4.8       3       0.0
dihe X    CH2P CH2P X    4.8       3       0.0
dihe X    CH2E N    X    0.9       3       0.0! SEE CH1E COMMENTS
dihe X    CH2P N    X    0.9       3       0.0! SEE CH1E COMMENTS
dihe X    CH2P NH3  X    0.9       3       0.0! added ATB 5/19/93
dihe X    CH2E NH1  X    0.9       3       0.0
dihe X    CH2G NH1  X    0.9       3       0.0
dihe X    CH2E NH2  X    1.8       3       0.0
dihe X    CH2E NH3  X    1.8       3       0.0
dihe X    CH2G NH3  X    1.8       3       0.0
dihe X    CH2E OH1  X    1.5       3       0.0
dihe X    CH2E S    X    3.6       2       0.0
dihe X    CH2E SM   X    3.6       2       0.0
!dihe X    CT   CT   X    1.6       3       0.0
!dihe X    CT   N    X    0.3       3       0.0! SEE CH1E COMMENTS
!dihe X    CT   NC2  X    0.3       3       0.0
!dihe X    CT   NH1  X    0.3       3       0.0
!dihe X    CT   NH2  X    0.6       3       0.0
!dihe X    CT   NH3  X    0.6       3       0.0
!dihe X    CT   OH1  X    0.5       3       0.0
!dihe X    CT   S    X    1.2       2       0.0
!dihe X    FE   NR   X    0.05      4       0.0
!dihe X    FE   CM   X    0.05      4       0.0
!dihe X    FE   OM   X    0.00      4       0.0
dihe X    S    S    X   12.0       2       0.0! FROM EXP.R BARRI
dihe CH2E NH1  C  NC2    300.0      2      180.    ! special dihedral for ARG


! SEE ABOVE REGARDING IMPROPERS
{ Weis 5/11/92 scaled up 3x for consistency with bond and angle parameters.}
{ Note: didn't change commented lines.                                     }
impr C    C    CR1E CH2E 750.0    0   0.0!!
impr C    CR1E C    CH2E 750.0    0   0.0!!
impr C5W  CR1E CW   CH2E 750.0    0   0.0!!
impr C    CR1E CR1E CH2E 750.0    0   0.0!! increased, ring torsions
impr CF   CR1E CR1E CH2E 750.0    0   0.0!! increased, ring torsions
impr CY   CR1E CR1E CH2E 750.0    0   0.0!! increased, ring torsions
impr C    CR1E NH1  CH2E 750.0    0   0.0!!
impr C5   NH1  CR1E CH2E 750.0    0   0.0!!
impr C5   NH1  CR1H CH2E 750.0    0   0.0!!
impr C5   NH1  CRHH CH2E 750.0    0   0.0!!
impr C    CR1E CR1E OH1  750.0    0   0.0!!
impr CY2  CR1E CR1E OH1  750.0    0   0.0!!
impr C    H    H    NH2  135.0    0   0.0! PRIMARY AMIDES (ASN AND GLN) OOP
impr C    OC   OC   CH1E 300.0    0   0.0! CARBOXYL OUT OF PLANE.
impr C    OC   OC   CH2E 300.0    0   0.0!
impr C    X    X    C     75.0    0   0.0! FROM BENZENE NORMAL MODE ANALYSIS
impr CY   X    X    CY2   75.0    0   0.0! FROM BENZENE NORMAL MODE ANALYSIS
impr C5W  X    X    CW    75.0    0   0.0! FROM BENZENE NORMAL MODE ANALYSIS
impr C    X    X    CH2E 270.0    0   0.0! FROM TOLUENE METHYL OOP. 217 CM 1
impr C5   X    X    CH2E 270.0    0   0.0! FROM TOLUENE METHYL OOP. 217 CM 1
impr C    X    X    CH3E 270.0    0   0.0
impr C   X    X    CR1E  750.0    0   0.0  !increased !!! ring torsions
impr C5  X    X    CRH   750.0    0   0.0  !increased !!! ring torsions
impr CF   X    X    CR1E 750.0    0   0.0
impr CW   X    X    CR1E 750.0    0   0.0
impr CW   X    X    CR1W 750.0    0   0.0
impr C5   X    X    CRHH 750.0    0   0.0
impr C    X    X    H    225.0    0   0.0! FROM BENZENE NORMAL MODE ANALYSIS
impr C    X    X    HA   225.0    0   0.0!
impr C    X    X    NH1  300.0    0   0.0! AMIDES FIT TO N METHYL ACETAMIDE.
impr CW   X    X    NH1  300.0    0   0.0! AMIDES FIT TO N METHYL ACETAMIDE.
impr C5W  X    X    NH1  300.0    0   0.0! AMIDES FIT TO N METHYL ACETAMIDE.
impr C5   X    X    NH1  300.0    0   0.0! AMIDES FIT TO N METHYL ACETAMIDE.
impr C    X    X    O    300.0    0   0.0
impr C    X    X    OC   300.0    0   0.0
impr C    X    X    OH1  450.0    0   0.0! USED FOR TYR HYDROXYL OOP
impr CH1E X    X    CH1E 300.0    0 35.26439 !
impr CH1E X    X    CH2E 300.0    0 35.26439 ! INCREASED  ! chirality
impr CH1E X    X    CH3E 300.0    0 35.26439 !
impr CR1E X    X    CR1E 750.0    0   0.0!  ! increased !! ring torsions
impr CR1E X    X    CRH  750.0    0   0.0!  ! increased !! ring torsions
impr CR1W X    X    CR1E 750.0    0   0.0! EXTENDED ATOM VERSION OF BENZENE
impr CR1H X    X    CRHH 750.0    0   0.0! EXTENDED ATOM VERSION OF BENZENE
impr CR1H X    X    NH1  750.0    0   0.0! SAME AS ABOVE FOR LACK OF VALUES
impr CR1E X    X    NH1  750.0    0   0.0!  ! increased  !! ring torsions
!impr FE   X    X    NP    20.0    0   0.0! FROM PARMFIX9
impr H    X    X    O    135.0    0   0.0
!impr N    CH1E CH2E C     45.0    0   0.0! PROLINE NITROGENS
impr N    CH1E CH2P C    135.0    0   0.0! PROLINE NITROGENS
impr N    X    X    CH2E 135.0    0   0.0
impr N    X    X    CT   135.0    0   0.0
impr NC2  X    X    CT   135.0    0   0.0
impr NC2  X    X    HC   135.0    0   0.0
impr NH1  X    X    CH1E 135.0    0   0.0
impr NH1  X    X    CH2E 135.0    0   0.0
impr NH1  X    X    CH3E 135.0    0   0.0
impr NH1  X    X    CT   135.0    0   0.0
impr NH1  X    X    H    135.0    0   0.0! AMIDES PROTON OOP
impr NH1  X    X    NH1   75.0    0   0.0!
impr NH1  X    X    NR    75.0    0   0.0
impr NH2  X    X    H    135.0    0   0.0
impr NR   X    X    C5    75.0    0   0.0
impr NR   X    X    CR1E  750.0    0   0.0   ! increased !! ring torsions
impr NR   X    X    CT    75.0    0   0.0
impr NR   X    X    CH3E  75.0    0   0.0 ! FOR NETROPSIN



{* nonbonding parameter section *}
{* ============================ *}
!! for use with:
!! NBXMOD=5  ATOM CDIEL SHIFT vswitch
!!    CUTNB=8.0  CTOFNB=7.5  CTONNB=6.5  EPS=1.0  E14FAC=0.4  WMIN=1.5
!!
{* nonbonding parameter section *}
{* ============================ *}
   nbonds                                    { This statement specifies the   }
      atom cdie shift eps=1.0  e14fac=0.4    { nonbonded interaction energy   }
      cutnb=7.5 ctonnb=6.0 ctofnb=6.5        { options.  Note the reduced     }
      nbxmod=5 vswitch  tolerance=0.5        { nonbonding cutoff to save      }
   end                                       { some CPU time                  }


 !                  eps     sigma       eps(1:4) sigma(1:4)
 !                  (kcal/mol) (A)
 !                  ---------------------------------------
 NONBonded  H       0.0498   1.4254      0.0498   1.4254
 NONBonded  HA      0.0450   2.6157      0.0450   2.6157 !- charged group.
 NONBonded  HC      0.0498   1.0691      0.0498   1.0691 !   Reduced vdw radius
 !
 NONBonded  C       0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  C5      0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  C5W     0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  CF      0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  CW      0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  CY      0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  CY2     0.1200   3.7418      0.1000   3.3854 ! carbonyl carbon
 NONBonded  CH1E    0.0486   4.2140      0.1000   3.3854 ! \
 NONBonded  CH2E    0.1142   3.9823      0.1000   3.3854 !  extended carbons
 NONBonded  CH2G    0.1142   3.9823      0.1000   3.3854 !  extended carbons
 NONBonded  CH2P    0.1142   3.9823      0.1000   3.3854 !  extended carbons
 NONBonded  CH3E    0.1811   3.8576      0.1000   3.3854 ! /
!! NONBonded  CM      0.0262   4.4367      0.1000   3.3854
 NONBonded  CR1E    0.1200   3.7418      0.1000   3.3854 !  ring carbons
 NONBonded  CR1H    0.1200   3.7418      0.1000   3.3854 !  ring carbons
 NONBonded  CR1W    0.1200   3.7418      0.1000   3.3854 !  ring carbons
 NONBonded  CRHH    0.1200   3.7418      0.1000   3.3854 !  ring carbons
 NONBonded  CRH     0.1200   3.7418      0.1000   3.3854 !  ring carbons
!! NONBonded  CT      0.0262   4.4367      0.1000   3.3854
 !
 NONBonded  N       0.2384   2.8509      0.2384   2.8509
 NONBonded  NC2     0.2384   2.8509      0.2384   2.8509
 NONBonded  NH1     0.2384   2.8509      0.2384   2.8509
 NONBonded  NH2     0.2384   2.8509      0.2384   2.8509
 NONBonded  NH3     0.2384   2.8509      0.2384   2.8509
 NONBonded  NP      0.2384   2.8509      0.2384   2.8509
 NONBonded  NR      0.2384   2.8509      0.2384   2.8509
 !
 NONBonded  O       0.1591   2.8509      0.1591   2.8509
 NONBonded  OC      0.6469   2.8509      0.6469   2.8509
 NONBonded  OH1     0.1591   2.8509      0.1591   2.8509
!! NONBonded  OM      0.1591   2.8509      0.1591   2.8509
 !
 NONBonded  S       0.0430   3.3676      0.0430   3.3676
 NONBonded  SM      0.0430   3.3676      0.0430   3.3676
 NONBonded  SH1E    0.0430   3.3676      0.0430   3.3676
 !
!! NONBONDED FE        0.0000    1.1582      0.0000 1.1582

set echo=true end
