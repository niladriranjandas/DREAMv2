remark - parameter file PROLSQ.PRO -
remark  Trying to match PROLSQ energy functions - ak 

set echo=false end

!! PROLSQ uses a standard deviation for energy parameters
!!  in general the form is (1 / sigma ^ 2) 
!!  so k-factors are calculated to match these standard deviations:

!!  Bond Length                           sigma = .02 A  
!!  Angle-related distance                sigma = .03 A
!!  Planar groups			  sigma = .02 A
!!  Chiral centers 			  sigma = .15 A
!!  Planar peptide w torsion angle	  sigma =   3 degrees

!! For the non-bonded energy term, the function is of the form
!!    (1 / sigma ^ 4)
!!  Nonbonded interactions                sigma = .5  A

evaluate ( $kbond       =  2500.0 ) !dimension: ( kcal / mole-A^2   )
evaluate ( $kangle      =   584.0 ) !dimension: ( kcal / mole-rad^2 )
evaluate ( $omega      =  1700.0 ) !dimension: ( kcal / mole       )
evaluate ( $pro_omega  =   100.0 ) !dimension: ( kcal / mole       )
evaluate ( $plane      =   360.0 ) !dimension: ( kcal / mole-rad^2 )
evaluate ( $chiral     =   796.0 ) !dimension: ( kcal / mole-rad^2 )

!! All bond paramaters are set to 2500.0 to match a standard deviation of
!!  sigma = .02 in PROLSQ

bond C    C      $kbond    1.38
bond C    CH1E   $kbond    1.52
bond C    CH2E   $kbond    1.52
bond C    CH3E   $kbond    1.52
bond C    CR1E   $kbond    1.38
bond C    N      $kbond    1.33
bond C    NC2    $kbond    1.33
bond C    NH1    $kbond    1.33
bond C    NH2    $kbond    1.33
bond C    NP     $kbond    1.33
bond C    NR     $kbond    1.33
bond C    O      $kbond    1.23
bond C    OC     $kbond    1.23
bond C    OH1    $kbond    1.38
bond C    OS     $kbond    1.43
bond CH1E CH1E   $kbond    1.53
bond CH1E CH2E   $kbond    1.52
bond CH1E CH3E   $kbond    1.52
bond CH1E N      $kbond    1.45
bond CH1E NH1    $kbond    1.45
bond CH1E NH2    $kbond    1.45
bond CH1E NH3    $kbond    1.45
bond CH1E OH1    $kbond    1.42
bond CH2E CH2E   $kbond    1.52
bond CH2E CH3E   $kbond    1.54
bond CH2E CR1E   $kbond    1.45
bond CH2E N      $kbond    1.45
bond CH2E NH1    $kbond    1.45
bond CH2E NH2    $kbond    1.45
bond CH2E NH3    $kbond    1.45
bond CH2E OH1    $kbond    1.42
bond CH2E S      $kbond    1.81
bond CH2E SH1E   $kbond    1.81
bond CH3E NH1    $kbond    1.49
bond CH3E NR     $kbond    1.49
bond CH3E S      $kbond    1.77
bond CH3E OS     $kbond    1.38
bond CM   OM     $kbond    1.128
bond CR1E CR1E   $kbond    1.38
bond CR1E NH1    $kbond    1.305
bond CR1E NR     $kbond    1.305
bond H    NH1    $kbond    0.98
bond H    NH2    $kbond    0.98
bond H    OH1    $kbond    0.96
bond HA   C      $kbond    1.08
bond HC   NC2    $kbond    1.00
bond HC   NH1    $kbond    0.98
bond HC   NH3    $kbond    1.04
bond OC   S      $kbond    1.43
bond OM   OM     $kbond    1.23
bond S    S      $kbond    2.02
      
!!  The energy parameter for angles was determinded by taking an ideal 
!!   situation of two bonds of length 1.45 A and an angle of 120 degrees
!!   and letting it deviate by 1.5 degrees in either direction and 
!!   taking the average of the two directions.

angle C    C    C      $kangle  106.5
angle C    C    CH2E   $kangle  126.5
angle C    C    CH3E   $kangle  126.5
angle C    C    CR1E   $kangle  122.5
angle C    C    HA     $kangle  120.0
angle C    C    NH1    $kangle  109.0
angle C    C    NP     $kangle  112.5
angle C    C    NR     $kangle  112.5
angle C    C    OH1    $kangle  119.0
angle C    C    O      $kangle  119.0
angle CH1E C    N      $kangle  117.5
angle CH1E C    NH1    $kangle  117.5
angle CH1E C    O      $kangle  121.5
angle CH1E C    OC     $kangle  117.5
angle CH1E C    OH1    $kangle  120.0
angle CH2E C    CR1E   $kangle  121.5
angle CH2E C    N      $kangle  117.5
angle CH2E C    NH1    $kangle  117.5
angle CH2E C    NH2    $kangle  117.5
angle CH2E C    NC2    $kangle  117.5
angle CH2E C    NR     $kangle  116.0
angle CH2E C    O      $kangle  121.6
angle CH2E C    OC     $kangle  118.5
angle CH2E C    OH1    $kangle  120.0
angle CH3E C    N      $kangle  117.5
angle CH3E C    NH1    $kangle  117.5
angle CH3E C    O      $kangle  121.5
angle CR1E C    CR1E   $kangle  120.5
angle CR1E C    NH1    $kangle  110.5
angle CR1E C    NP     $kangle  122.5
angle CR1E C    NR     $kangle  122.5
angle CR1E C    OH1    $kangle  119.0
angle N    C    O      $kangle  121.0
angle NC2  C    NC2    $kangle  120.0
angle NC2  C    NH1    $kangle  120.0
angle NH1  C    NR     $kangle  120.0
angle NH1  C    O      $kangle  121.0
angle NH2  C    O      $kangle  121.0
angle O    C    OH1    $kangle  120.0
angle OC   C    OC     $kangle  122.5
angle OS   C    CH1E   $kangle  125.3
angle OS   C    CH2E   $kangle  125.3
angle OS   C    O      $kangle  120.0
angle C    CH1E CH1E   $kangle  110.0
angle C    CH1E CH2E   $kangle  109.5
angle C    CH1E CH3E   $kangle  106.5
angle C    CH1E N      $kangle  111.6
angle C    CH1E NH1    $kangle  111.6
angle C    CH1E NH2    $kangle  111.6
angle C    CH1E NH3    $kangle  111.6
angle CH1E CH1E CH2E   $kangle  112.5
angle CH1E CH1E CH3E   $kangle  111.0
angle CH1E CH1E NH1    $kangle  110.0
angle CH1E CH1E NH2    $kangle  109.5
angle CH1E CH1E NH3    $kangle  107.5
angle CH1E CH1E OH1    $kangle  104.5
angle CH2E CH1E CH3E   $kangle  111.5
angle CH2E CH1E N      $kangle  104.0
angle CH2E CH1E NH1    $kangle  110.0
angle CH2E CH1E NH2    $kangle  110.0
angle CH2E CH1E NH3    $kangle  110.0
angle CH3E CH1E CH3E   $kangle  111.0
angle CH3E CH1E NH1    $kangle  108.5
angle CH3E CH1E NH2    $kangle  109.5
angle CH3E CH1E NH3    $kangle  109.5
angle CH3E CH1E OH1    $kangle  110.5
angle C    CH2E CH1E   $kangle  112.5
angle C    CH2E CH2E   $kangle  113.0
angle C    CH2E NH1    $kangle  111.6
angle C    CH2E NH2    $kangle  111.6
angle C    CH2E NH3    $kangle  111.6
angle CH1E CH2E CH1E   $kangle  117.0
angle CH1E CH2E CH2E   $kangle  112.5
angle CH1E CH2E CH3E   $kangle  113.0
angle CH1E CH2E OH1    $kangle  111.0
angle CH3E CH2E OH1    $kangle  111.0
angle CH1E CH2E S      $kangle  112.5
angle CH1E CH2E SH1E   $kangle  112.5
angle CH2E CH2E CH2E   $kangle  110.0
angle CH2E CH2E CH3E   $kangle  111.0
angle CH2E CH2E N      $kangle  105.0
angle CH2E CH2E NH1    $kangle  111.0
angle CH2E CH2E NH2    $kangle  109.5
angle CH2E CH2E NH3    $kangle  110.5
angle CH2E CH2E S      $kangle  112.5
angle C    CR1E C      $kangle  126.5
angle C    CR1E CH2E   $kangle  122.0
angle C    CR1E CR1E   $kangle  119.0
angle C    CR1E NH1    $kangle  109.5
angle C    CR1E NR     $kangle  106.5
angle CR1E CR1E CR1E   $kangle  120.5
angle NH1  CR1E NH1    $kangle  109.0
angle NH1  CR1E NR     $kangle  109.0
angle C    N    CH1E   $kangle  120.0
angle C    N    CH2E   $kangle  120.0
angle C    N    CT     $kangle  120.0
angle CH1E N    CH2E   $kangle  110.0
angle CH1E N    CH3E   $kangle  110.0
angle CH2E N    CH3E   $kangle  109.5
angle CT   N    CT     $kangle  110.0
angle C    NC2  CT     $kangle  120.0
angle C    NC2  HC     $kangle  120.0
angle HC   NC2  HC     $kangle  120.0
angle C    NH1  C      $kangle  102.5
angle C    NH1  CH1E   $kangle  120.0
angle C    NH1  CH2E   $kangle  120.0
angle C    NH1  CH3E   $kangle  120.0
angle C    NH1  CR1E   $kangle  108.0
angle C    NH1  H      $kangle  120.0
angle CH1E NH1  CH3E   $kangle  120.0
angle CH1E NH1  H      $kangle  120.0
angle CH2E NH1  CH3E   $kangle  120.0
angle CH2E NH1  H      $kangle  120.0
angle CH3E NH1  H      $kangle  120.0
angle CR1E NH1  CR1E   $kangle  110.0
angle CR1E NH1  H      $kangle  120.0
angle C    NH2  H      $kangle  120.0
angle CH1E NH2  CH2E   $kangle  120.0
angle CH1E NH2  H      $kangle  120.0
angle CH2E NH2  H      $kangle  120.0
angle H    NH2  H      $kangle  125.0
angle C    NP   C      $kangle  102.5
angle C    NR   C      $kangle  102.5
angle C    NR   CR1E   $kangle  109.5
angle CH3E NR   CR1E   $kangle  109.5 
angle CH3E NR   C      $kangle  109.5
angle CR1E NR   CR1E   $kangle  110.0
angle CH1E NH3  HC     $kangle  109.5
angle CH1E NH3  CH2E   $kangle  109.5
angle CH2E NH3  HC     $kangle  109.5
angle HC   NH3  HC     $kangle  109.5
angle C    OH1  H      $kangle  109.5
angle CH1E OH1  H      $kangle  109.5
angle CH2E OH1  H      $kangle  109.5
angle C    OS   CH3E   $kangle  120.5
angle CH2E S    CH3E   $kangle   99.5
angle CH2E S    S      $kangle  104.2
angle OC   S    OC     $kangle  109.5

!! The only dihedral term that is not being set to zero is the term that
!!  preserves planarity of the peptide bond.  The constant came from a 
!!  comparison of the strength of the two different energy functions under 
!!  small deviations.

dihe X    C    NH1  X    $omega     1       0.0    
dihe CH1E C    N    CH1E $pro_omega 2     180.0
dihe CH2E C    N    CH1E $pro_omega 2     180.0

dihe CR1E C    C    CR1E $pro_omega 2     180.0  !fixed, ATB
dihe CR1E C    C    C    $pro_omega 2     180.0
dihe CR1E C    C    NH1  $pro_omega 2     180.0
dihe X    C    CH1E X    0.0       3       0.0
dihe X    C    CH2E X    0.0       3       0.0
dihe X    C    CR1E X    0.0       2     180.0
dihe X    C    CT   X    0.0       3       0.0
dihe X    C    N    X    0.0       2     180.0
dihe X    C    NC2  X    0.0       2     180.0
dihe X    C    NH2  X    0.0       2     180.0
dihe X    C    OH1  X    0.0       2     180.0
dihe X    C    OS   X    0.0       2     180.0
dihe X    CH1E CH1E X    0.0       3       0.0
dihe X    CH1E CH2E X    0.0       3       0.0
dihe X    CH1E N    X    0.0       3       0.0
dihe X    CH1E NH1  X    0.0       3       0.0
dihe X    CH1E NH2  X    0.0       3       0.0
dihe X    CH1E NH3  X    0.0       3       0.0
dihe X    CH1E OH1  X    0.0       3       0.0
dihe X    CH2E CH2E X    0.0       3       0.0
dihe X    CH2E N    X    0.0       3       0.0
dihe X    CH2E NH1  X    0.0       3       0.0
dihe X    CH2E NH2  X    0.0       3       0.0
dihe X    CH2E NH3  X    0.0       3       0.0
dihe X    CH2E OH1  X    0.0       3       0.0
dihe X    CH2E S    X    0.0       2       0.0
dihe X    CT   CT   X    0.0       3       0.0
dihe X    CT   N    X    0.0       3       0.0
dihe X    CT   NC2  X    0.0       3       0.0
dihe X    CT   NH1  X    0.0       3       0.0
dihe X    CT   NH2  X    0.0       3       0.0
dihe X    CT   NH3  X    0.0       3       0.0
dihe X    CT   OH1  X    0.0       3       0.0
dihe X    CT   S    X    0.0       2       0.0
dihe X    FE   NR   X    0.0       4       0.0
dihe X    FE   CM   X    0.0       4       0.0
dihe X    FE   OM   X    0.0       4       0.0
dihe X    S    S    X    0.0       2       0.0



!! There are only two kinds of improper terms that are not set to zero:
!!  the terms that preserve the planarity of the rings and the terms that
!!   preserve the chirality of tetrahedral atoms.  The chiral constant comes
!!   from comparing the relative strength of the two functions when the 
!!   chirality is completely reversed.  The planarity constant comes from
!!   a comparison of the two functions under small deviaions from a perfect
!!   plane.
impr C    C    CR1E CH2E  $plane    0   0.0
impr C    CR1E C    CH2E  $plane    0   0.0
impr C    CR1E CR1E CH2E  $plane    0   0.0
impr C    CR1E NH1  CH2E  $plane    0   0.0
impr C    NH1  CR1E CH2E  $plane    0   0.0
impr C    CR1E CR1E OH1   $plane    0   0.0
impr C    X    X    C     $plane    0   0.0
impr C    X    X    CH2E  $plane    0   0.0
impr C    X    X    CH3E  $plane    0   0.0
impr C    X    X    CR1E  $plane    0   0.0
impr CR1E X    X    CR1E  $plane    0   0.0
impr CR1E X    X    NH1   $plane    0   0.0
impr NR   X    X    CR1E  $plane    0   0.0
impr CH1E X    X    CH1E  $chiral   0  35.26439 
impr CH1E X    X    CH2E  $chiral   0  35.26439 
impr CH1E X    X    CH3E  $chiral   0  35.26439 

impr C    H    H    NH2   0.0    0   0.0
impr C    OC   OC   CH1E  0.0    0   0.0
impr C    OC   OC   CH2E  0.0    0   0.0
impr C    X    X    H     0.0    0   0.0
impr C    X    X    HA    0.0    0   0.0
impr C    X    X    NH1   0.0    0   0.0
impr C    X    X    O     0.0    0   0.0
impr C    X    X    OC    0.0    0   0.0
impr C    X    X    OH1   0.0    0   0.0
impr H    X    X    O     0.0    0   0.0
impr N    CH1E CH2E C     0.0    0   0.0
impr N    X    X    CH2E  0.0    0   0.0
impr N    X    X    CT    0.0    0   0.0
impr NC2  X    X    CT    0.0    0   0.0
impr NC2  X    X    HC    0.0    0   0.0
impr NH1  X    X    CH1E  0.0    0   0.0
impr NH1  X    X    CH2E  0.0    0   0.0
impr NH1  X    X    CH3E  0.0    0   0.0
impr NH1  X    X    CT    0.0    0   0.0
impr NH1  X    X    H     0.0    0   0.0
impr NH1  X    X    NH1   0.0    0   0.0 
impr NH1  X    X    NR    0.0    0   0.0
impr NH2  X    X    H     0.0    0   0.0
impr NR   X    X    C     0.0    0   0.0
impr NR   X    X    CT    0.0    0   0.0
impr NR   X    X    CH3E  0.0    0   0.0 


{* nonbonding parameter section *}
{* ============================ *}
!!
!  This uses a new form of the REPEL function:
!    fVDW(R) =  RCON *( Rmin ^ IREX - R  ^ IREX ) ^ REXP
!
!  PROLSQ uses a function of the form:
!    fVDW(R) =  (1 / 0.5) ^ 4 * ( Rmin ^ 4 - R  ^ 4 ) 
!
!  The epsilon values are arbitrary since the repel function does not depend
!   on epsilon.  The sigma values come from converting the Van der Waals 
!   radii of the PROLSQ program into sigma values using the formula:
!     Rmin = sigma * 2 ^ (1/6)
!   Note:  Prolsq increments Van der Waals radii for non-bonded contacts
!           that involve torsion angles (1:4 contacts) by -.30 A.
!          This accounts for the difference between the sigma values for
!           regular contacts and (1:4) contacts.

 NBONds
  CUTNB=5.0   EPS=1.0  E14FAC=0.4  WMIN=1.5
  REPEl = 1.0    !            turns on the repel function
  REXPonent = 4
  IREXponent = 1
  RCONst = 16.0
  TOLErance = 0.5      NBXMOD = 5
END
 !                  eps     sigma       eps(1:4) sigma(1:4)
 !                  (kcal/mol) (A)
 !                  ---------------------------------------
 NONBonded  H       0.1000   2.1381      0.1000   1.8709
 NONBonded  HC      0.1000   2.1381      0.1000   1.8709 
 !
 NONBonded  C       0.1000   3.0291      0.1000   2.7618
 NONBonded  CH1E    0.1000   3.2963      0.1000   3.0291
 NONBonded  CH2E    0.1000   3.2963      0.1000   3.0291
 NONBonded  CH3E    0.1000   3.2963      0.1000   3.0291
 NONBonded  CR1E    0.1000   3.2963      0.1000   3.0291
 !
 NONBonded  N       0.1000   2.7618      0.1000   2.4945
 NONBonded  NC2     0.1000   2.7618      0.1000   2.4945
 NONBonded  NH1     0.1000   2.7618      0.1000   2.4945
 NONBonded  NH2     0.1000   2.7618      0.1000   2.4945
 NONBonded  NH3     0.1000   2.7618      0.1000   2.4945
 NONBonded  NP      0.1000   2.7618      0.1000   2.4945
 NONBonded  NR      0.1000   2.7618      0.1000   2.4945
 !
 NONBonded  O       0.1000   2.6727      0.1000   2.4054
 NONBonded  OC      0.1000   2.6727      0.1000   2.4054
 NONBonded  OH1     0.1000   2.6727      0.1000   2.4054
 !
 NONBonded  S       0.1000   3.2072      0.1000   2.9400
 NONBonded  SH1E    0.1000   3.2072      0.1000   2.9400
 !


set echo=true end

