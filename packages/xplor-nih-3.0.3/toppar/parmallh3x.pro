remarks   PARAMETER FILE FOR ALL EXPLICIT HYDROGENS BASED ON CFF PARAMETERS
remarks   IN BRUCE GELIN'S THESIS (PARMFIX7 VALUES HAVE TRAILING 0'S)
remarks   BOND ANGLE PARAMETERS INCLUDE THE HARMONIC CONTRIBUTION FROM
remarks   Urey-Bradley F TERM ( (R13-R0)\2 )
remarks   DIHEDRALS INCLUDE MULTIPLICITY FROM CFF
remarks   MODIFIED TO INCLUDE HA AND CT
remarks   HIS AND TRYP INCLUDED 4/13/80 DJS
remarks   hbond parameters as in PARAM7 ( ATB, 3-APR-85 )
remarks 
remarks   modified: increased force constant on cis/trans peptide
remarks   dihedrals
remarks   increased improper force constants
set echo=false end
 
bonds H    NH2    405.0000    0.98
bonds H    NH1    405.0000    0.98
bonds H    OH1    450.0000    0.96
bonds H    S      450.0000    0.96           !added, ATB and GMC, 30-SEP-85
bonds H    OH2    450.0000    0.96
bonds HA   CT     300.0000    1.08
bonds HA   C      350.000     1.08
bonds HC   NC2    405.0000    1.0
bonds HC   NH1    405.000     0.98
bonds HC   NH3    405.0000    1.04
bonds C    C      250.0000    1.38
bonds C    CT     187.0000    1.53
bonds C    N      403.0000    1.305
bonds C    NP     403.0000    1.305
bonds C    NR     403.0000    1.305
bonds C    NH1    403.0000    1.305
bonds C    NH2    403.0000    1.305
bonds C    NC2    403.0000    1.305
bonds C    O      595.0000    1.215
bonds C    OC     450.0000    1.22
bonds C    OH1    450.0000    1.38
bonds CT   CT     110.0000    1.53
bonds CT   N      261.0000    1.49
bonds CT   NH1    261.0000    1.49
bonds CT   NH2    261.0000    1.49
bonds CT   NH3    261.0000    1.49
bonds CT   NC2    261.000     1.49
bonds CT   OH1    400.0000    1.42
bonds CT   S      450.0000    1.8100
bonds S    S      500.0000    2.0200
 
 
angle H    NH1  H       30.000   107.500
angle H    NH1  C       35.3     120.0
angle H    NH1  CT      40.4     120.0
angle H    NH2  H       40.0000  120.0000
angle H    NH2  C       35.3     120.0
angle H    NH2  CT      40.4     120.0
angle H    OH1  CT      35.0000  108.0000
angle H    S    CT      35.0000  108.0000   ! added, ATB and GMC, 30-SEP-85
angle H    OH1  C       35.0000  108.0000
angle H    OH2  H       120.000  104.5000
angle HC   NH3  HC      40.0000  109.5000
angle HC   NH3  CT      40.0000  109.5000
angle HC   NC2  HC      40.0000  120.0000
angle HC   NC2  C       35.3     120.0
angle HC   NC2  CT      40.4     107.5
angle HA   C    C       41.0     120.0
angle HA   C    NH1     41.0     120.0
angle HA   C    NH2     41.0     120.0
angle HA   C    NR      41.0     120.0
angle HA   C    O       41.0     120.0
angle HA   CT   HA      39.3     107.5
angle HA   CT   C       49.3     107.5
angle HA   CT   CT      30.000   108.000
angle HA   CT   N       51.5     107.5
angle HA   CT   NH1     51.5     107.5
angle HA   CT   NH3     51.5     107.5
angle HA   CT   NC2     51.5     107.5
angle HA   CT   OH1     50.0     107.5
angle HA   CT   S       40.000   107.500
angle C    C    C       60.0000  126.5000
angle C    C    CT      45.8     122.3
angle C    C    NH1     50.0000  108.6000
angle C    C    NP      70.0000  109.5000
angle C    C    NR      70.0000  109.5000
angle C    C    OH1     45.2     122.300
angle C    CT   CT      51.8     107.5
angle C    CT   N       43.7     108.6
angle C    CT   NH1     43.7     108.6
angle C    CT   NH2     43.7     108.0
angle C    CT   NH3     43.7     108.0
angle C    N    CT      62.3     119.100
angle C    NH1  C       62.3     119.100
angle C    NH1  CT      62.3     119.1000
angle C    NC2  CT      62.3     107.500
angle C    NP   C       50.0000  106.0000
angle C    NR   C       50.0000  106.0000
angle CT   C    N       49.5     115.7000
angle CT   C    NH1     49.5     117.5000
angle CT   C    NH2     49.5     120.000
angle CT   C    O       64.6     122.3000
angle CT   C    OH1     35.0000  117.5000
angle CT   C    OC      40.0000  118.0000
angle CT   CT   CT      51.1     107.500
angle CT   CT   N       67.7     107.500
angle CT   CT   NH1     67.7     107.1000
angle CT   CT   NH2     67.7     108.9000
angle CT   CT   NH3     67.7     108.9000
angle CT   CT   NC2     67.7     107.500
angle CT   CT   OH1     50.0000  101.5000
angle CT   CT   S       50.0000  117.2000
angle CT   N    CT      60.000   119.000
angle CT   S    CT      50.0000   97.2000
angle CT   S    S       50.0000  104.2000
angle N    C    O       84.1     124.5000
angle NH1  C    NC2     60.0000  117.0000
angle NH1  C    NR      60.0000  117.0000
angle NH1  C    O       60.0000  124.5000
angle NH2  C    O       60.0000  120.6000
angle NC2  C    NC2     60.0000  120.3000
angle O    C    OH1     60.0000  124.5000
angle OC   C    OC      50.0000  129.0000
 
dihedral X    C    CT   X        1.2       3 0.0000
dihedral X    C    N    X       5.0       2 180.0 !reduced to allow cis-PRO
dihedral X    C    NH1  X       300.0      2 180.0 ! increased to avoid cis-pep
dihedral X    C    NH2  X       11.6       2 180.000
dihedral X    C    NC2  X       11.6       2 180.0000
dihedral X    C    OH1  X        0.4       2 180.0000
dihedral X    CT   CT   X        1.8       3 0.0000
dihedral X    CT   N    X        2.1       3 180.0000
dihedral X    CT   NH1  X        2.1       3 180.0000
dihedral X    CT   NH2  X        1.8       2 180.0000
dihedral X    CT   NH3  X        1.8       3 180.0000
dihedral X    CT   NC2  X        2.1       3 0.000
dihedral X    CT   OH1  X        0.6       3 180.0000
dihedral X    CT   S    X        1.5       2 0.0000
dihedral X    S    S    X        4.0000    2 0.0000
 
improper H    X    X    C        500.0 0 0.000   !!
improper H    X    X    NH1      500.0 0 0.0000  !! modified
improper H    X    X    NH2      500.0 0 0.0000  !!
improper H    X    X    O        500.0 0 0.0     !!
improper HC   X    X    NH1      500.0 0 0.000   !!
improper HC   X    X    NC2      500.0 0 0.0000  !!
improper HC   X    X    NH3      500.0 0 0.0     !!
improper HA   X    X    C        500.0 0 0.000   !!
improper CT   C    N    CT       500.0 0 35.26439  !! added for proline 
improper CT   C    NH1  CT       500.0 0 35.26439  !! added
improper CT   C    NH3  CT       500.0 0 35.26439  !! added
improper C    X    X    C        500.0  0 0.0000
improper C    X    X    CT       500.0 0 0.0
improper C    X    X    N        500.0 0 0.000
improper C    X    X    NH1      500.0 0 0.0000
improper C    X    X    NH2      500.0 0 0.000
improper C    X    X    NC2      500.0 0 0.000
improper C    X    X    NR       500.0 0 0.000
improper C    X    X    O        500.0 0 0.0000
improper C    X    X    OH1      500.0 0 0.0000
improper C    X    X    OC       500.0 0 0.0000
improper CT   X    X    N        500.0 0 0.000
improper CT   X    X    NH1      500.0 0 0.000
improper CT   X    X    NH3      500.0 0 0.000
improper CT   X    X    NC2      500.0 0 0.000
improper CT   X    X    NR       500.0 0 0.0000
improper NH1  X    X    NH1      500.0 0 0.0000
improper NH1  X    X    NR       500.0 0 0.0000
!
!! NONBONDED  NBXMOD=5  ATOM RDIEL SWITCH VSWITCH
!!     CUTNB=8.0  CTOFNB=7.5  CTONNB=6.5  EPS=1.0  E14FAC=1.0  WMIN=1.5
!
 !                  eps     sigma       eps(1:4) sigma(1:4)
 NONBonded  C       0.0903   3.2072      0.0903   3.2072
 NONBonded  CT      0.0903   3.2072      0.0903   3.2072
 NONBonded  H       0.0498   1.4254      0.0498   1.4254
 NONBonded  HA      0.0045   2.6157      0.0045   2.6157
 NONBonded  HC      0.0498   1.4254      0.0498   1.4254
 NONBonded  N       0.1592   2.7618      0.1592   2.7618
 NONBonded  NC2     0.1592   2.7618      0.1592   2.7618
 NONBonded  NH1     0.1592   2.7618      0.1592   2.7618
 NONBonded  NH2     0.1592   2.7618      0.1592   2.7618
 NONBonded  NH3     0.1592   2.7618      0.1592   2.7618
 NONBonded  NP      0.1592   2.7618      0.1592   2.7618
 NONBonded  NR      0.1592   2.7618      0.1592   2.7618
 NONBonded  O       0.2342   2.6406      0.2342   2.6406
 NONBonded  OC      1.0244   2.6406      1.0244   2.6406
 NONBonded  OH1     0.2342   2.6406      0.2342   2.6406
 NONBonded  OH2     0.2342   2.6406      0.2342   2.6406
 NONBonded  S       0.0239   3.3854      0.0239   3.3854
 
!! H        0.0440    1.0000    0.8000
!! HC       0.0440    1.0000    0.8000
!! HA       0.1       1.000     1.468
!! C        0.98      5.0000    1.8
!! CT       0.98      5.0000    1.8 
!! NP       0.74      6.0000    1.55
!! NR       0.74      6.0000    1.55
!! N        0.74      6.0000    1.55
!! NH1      0.74      6.0000    1.55
!! NH2      0.74      6.0000    1.55
!! NH3      0.74      6.0000    1.55
!! NC2      0.74      6.0000    1.55
!! O        0.8       6.0000    1.482
!! OH1      0.8       6.0000    1.482
!! OH2      0.8       6.0000    1.482
!! OC       2.1400    6.0000    1.482
!! S        0.3400   16.0000    1.9000
!! FE       0.0100    0.0001    0.6500
!
!! HBOND AEXP=4 REXP=6 HAEX=4 AAEX=2   ACCEPTORS
!!     CUTHB=4.5 CTOFHB=4.0 CTONHB=3.5  CUTHA=90.0  CTOFHA=70.0  CTONHA=50.0
!
AEXP 4
REXP 6
HAEX 4
AAEX 2
!                   Emin      Rmin
!                (Kcal/mol)   (A)
hbond N*+* N%      -3.00      3.0!  VALUES FROM VINOGRADOV AND LINELL FOR
hbond N*+* O*      -3.50      2.9!  TYPICAL LENGTHS AND DEPTHS.
hbond OH*  N%      -4.00      2.85
hbond OH*  O*      -4.25      2.75
hbond S    N%      -3.00      3.0 !! added, ATB
hbond S    O*      -3.50      2.9 !! added, ATB
 
set echo=true end
