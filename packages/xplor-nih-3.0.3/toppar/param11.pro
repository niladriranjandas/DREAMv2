REMARKS PARAM11.PRO ( from PARAM6A )
REMARKS ===========
REMARKS PROTEIN PARAMETERS:
REMARKS PEPTIDE GEOMETRY FROM RAMACHANDRAN ET AL BBA 359:298 (1974)
REMARKS TORSIONS FROM HAGLER ET AL JACS 98:4600 (1976)
REMARKS LENNARD-JONES NONBONDED PARAMETERS WITH SPECIAL TREATMENT OF 1:4
REMARKS CARBON-CARBON INTERACTIONS: JORGENSON ET. AL. 
REMARKS                           JACS 103:3976-3985 WITH 1-4 RC=1.80/0.1
REMARKS
set echo=false end

!             Eb = Kb ( r - r0 ) **2
!
!                 Kb          r0
!            (Kcal/mol A**2)   A
!
BOND C    C      450.0       1.38!  FROM B. R. GELIN THESIS AMIDE AND DIPEPTIDES
BOND C    CH1E   405.0       1.52!  EXCEPT WHERE NOTED.  CH1E,CH2E,CH3E, AND CT
BOND C    CH2E   405.0       1.52!  ALL TREATED THE SAME. UREY BRADLEY TERMS ADDED
BOND C    CH3E   405.0       1.52
BOND C    CR1E   450.0       1.38
!BOND C    CT     405.0       1.53
BOND C    N      471.0       1.33
BOND C    NC2    400.0       1.33!  BOND LENGTH FROM PARMFIX9 FORCE K APROXIMATE
!BOND C    NC2E   400.0       1.33
BOND C    NH1    471.0       1.33
!BOND C    NH1E   471.0       1.33
BOND C    NH2    471.0       1.33
!BOND C    NH2E   471.0       1.33
BOND C    NP     471.0       1.33
BOND C    NR     471.0       1.33
BOND C    O      580.0       1.23
BOND C    OC     580.0       1.23!  FORCE DECREASE AND LENGTH INCREASE FROM C O
BOND C    OH1    450.0       1.38!  FROM PARMFIX9 (NO VALUE IN GELIN THESIS)
!BOND C    OH1E   450.0       1.38
BOND CH1E CH1E   225.0       1.53
BOND CH1E CH2E   225.0       1.52
BOND CH1E CH3E   225.0       1.52
BOND CH1E N      422.0       1.45
BOND CH1E NH1    422.0       1.45
!BOND CH1E NH1E   422.0       1.45
BOND CH1E NH2    422.0       1.45
!BOND CH1E NH2E   422.0       1.45
BOND CH1E NH3    422.0       1.45
!BOND CH1E NH3E   422.0       1.45
BOND CH1E OH1    400.0       1.42!  FROM PARMFIX9 (NO VALUE IN GELIN THESIS)
!BOND CH1E OH1E   400.0       1.42
BOND CH2E CH2E   225.0       1.52
BOND CH2E CH3E   225.0       1.54
BOND CH2E CR1E   250.0       1.45!  FROM WARSHEL AND KARPLUS 1972 JACS 96:5612
BOND CH2E N      422.0       1.45
BOND CH2E NH1    422.0       1.45
!BOND CH2E NH1E   422.0       1.45
BOND CH2E NH2    422.0       1.45
!BOND CH2E NH2E   422.0       1.45
BOND CH2E NH3    422.0       1.45
!BOND CH2E NH3E   422.0       1.45
BOND CH2E OH1    400.0       1.42
!BOND CH2E OH1E   400.0       1.42
BOND CH2E S      450.0       1.81!  FROM PARMFIX9
BOND CH2E SH1E   450.0       1.81
BOND CH3E NH1    422.0       1.49
!BOND CH3E NH1E   422.0       1.49
BOND CH3E S      450.0       1.77!  FROM PARMFIX9
!BOND CM   OM    1115.0       1.128!  FROM CAUGHEY ET AL(1978),CARBON MONOXIDE
BOND CR1E CR1E   450.0       1.38
BOND CR1E NH1    450.0       1.305
!BOND CR1E NH1E   450.0       1.305
BOND CR1E NR     450.0       1.305
!BOND CT   CT     200.0       1.53
!BOND CT   N      422.0       1.45
!BOND CT   NC2    422.0       1.45
!BOND CT   NH1    422.0       1.45
!BOND CT   NH2    422.0       1.45
!BOND CT   NH3    422.0       1.45
!BOND CT   OH1    400.0       1.42
!BOND CT   S      450.0       1.81
!BOND FE   CM     258.0       1.79!   FROM KROEKER ET AL(JCP:72:4846)
!BOND FE   NP     500.0       2.09
!BOND FE   NR      65.0       1.98!   FROM NAGAI ET AL(1980)
!BOND FE   OM     250.0       1.8!    JUST A GUESS.
BOND H    NH1    405.0       0.98!  GELIN AND IR STRETCH 3200 CM 1
BOND H    NH2    405.0       0.98
BOND H    OH1    450.0       0.96!  FROM IR STRETCH 3400 CM 1
!BOND HA   C      350.0       1.08
!BOND HA   CT     300.0       1.08
BOND HC   NC2    405.0       1.00
BOND HC   NH1    405.0       0.98
BOND HC   NH3    405.0       1.04
BOND OC   S      400.0       1.43
!BOND OM   OM     600.0       1.23!   STRETCHING CONSTANT JUST A GUESS.
BOND S    S      500.0       2.02

!               Et = Kt ( theta - theta0 )**2
!
!                       Kt       theta0
!                (Kcal/mol A**2)  deg.
!
ANGLE C    C    C       70.0     106.5!  FROM B. R. GELIN THESIS WITH HARMONIC
ANGLE C    C    CH2E    65.0     126.5!  PART OF F TERMS INCORPORATED. ATOMS
ANGLE C    C    CH3E    65.0     126.5!  WITH EXTENDED H COMPENSATED FOR LACK
ANGLE C    C    CR1E    70.0     122.5!  OF H ANGLES.
!ANGLE C    C    CT      70.0     126.5
!ANGLE C    C    HA      40.0     120.0!  AMIDE PARAMETERS FIT BY LEAST SQUARES
ANGLE C    C    NH1     65.0     109.0!  TO N-METHYL ACETAMIDE VIBRATIONS.
!ANGLE C    C    NH1E    65.0     109.0!  GEOMETRY TO GIVE RAMACHANDRAN GEOM. ON
ANGLE C    C    NP      65.0     112.5!  MINIMIZATION OF N-METHYL ACETAMIDE.
ANGLE C    C    NR      65.0     112.5
ANGLE C    C    OH1     65.0     119.0
ANGLE C    CH1E CH1E    70.0     110.0
ANGLE C    CH1E CH2E    70.0     109.5
ANGLE C    CH1E CH3E    70.0     106.5
ANGLE C    CH1E N       45.0     111.6
ANGLE C    CH1E NH1     45.0     111.6
!ANGLE C    CH1E NH1E    45.0     111.6
ANGLE C    CH1E NH2     45.0     111.6
!ANGLE C    CH1E NH2E    45.0     111.6
ANGLE C    CH1E NH3     45.0     111.6
!ANGLE C    CH1E NH3E    45.0     111.6
ANGLE C    CH2E CH1E    70.0     112.5
ANGLE C    CH2E CH2E    70.0     113.0
ANGLE C    CH2E NH1     70.0     111.6
!ANGLE C    CH2E NH1E    70.0     111.6
ANGLE C    CH2E NH2     70.0     111.6
!ANGLE C    CH2E NH2E    70.0     111.6
ANGLE C    CH2E NH3     70.0     111.6
!ANGLE C    CH2E NH3E    70.0     111.6
ANGLE C    CR1E C       90.0     126.5
ANGLE C    CR1E CH2E    90.0     122.0
ANGLE C    CR1E CR1E    90.0     119.0
ANGLE C    CR1E NH1     90.0     109.5
!ANGLE C    CR1E NH1E    90.0     109.5
ANGLE C    CR1E NR      90.0     106.5
!ANGLE C    CT   CT      70.0     109.5
!ANGLE C    CT   HA      70.0     109.5
!ANGLE C    CT   N       70.0     111.6
!ANGLE C    CT   NH1     70.0     111.6
!ANGLE C    CT   NH2     70.0     111.6
!ANGLE C    CT   NH3     70.0     111.6
ANGLE C    N    CH1E    80.0     120.0
ANGLE C    N    CH2E    80.0     120.0
!ANGLE C    N    CT      80.0     120.0
!ANGLE C    NC2  CT      80.0     120.0
ANGLE C    NC2  HC      35.0     120.0
ANGLE C    NH1  C       60.0     102.5
ANGLE C    NH1  CH1E    77.5     120.0
ANGLE C    NH1  CH2E    77.5     120.0
ANGLE C    NH1  CH3E    77.5     120.0
ANGLE C    NH1  CR1E    60.0     108.0
!ANGLE C    NH1  CT      80.0     120.0
ANGLE C    NH1  H       30.0     120.0
!ANGLE C    NH1E CH1E    77.5     120.0
!ANGLE C    NH1E CH2E    77.5     120.0
!ANGLE C    NH1E CH3E    77.5     120.0
!ANGLE C    NH1E CR1E    60.0     108.0
ANGLE C    NH2  H       30.0     120.0
ANGLE C    NP   C       70.0     102.5
!ANGLE C    NP   FE      50.0     128.0!  FORCE CONSTANT FROM PARMFIX9
ANGLE C    NR   C       70.0     102.5
ANGLE C    NR   CR1E    70.0     109.5
ANGLE C    OH1  H       50.0     109.5
ANGLE CH1E C    N       20.0     117.5
ANGLE CH1E C    NH1     20.0     117.5
!ANGLE CH1E C    NH1E    20.0     117.5
ANGLE CH1E C    O       85.0     121.5
ANGLE CH1E C    OC      85.0     117.5
ANGLE CH1E C    OH1     85.0     120.0
!ANGLE CH1E C    OH1E    85.0     120.0
ANGLE CH1E CH1E CH2E    45.0     112.5
ANGLE CH1E CH1E CH3E    45.0     111.0
ANGLE CH1E CH1E NH1     50.0     110.0
!ANGLE CH1E CH1E NH1E    50.0     110.0
ANGLE CH1E CH1E NH2     50.0     109.5
!ANGLE CH1E CH1E NH2E    50.0     109.5
ANGLE CH1E CH1E NH3     50.0     107.5
!ANGLE CH1E CH1E NH3E    50.0     107.5
ANGLE CH1E CH1E OH1     50.0     104.5
!ANGLE CH1E CH1E OH1E    50.0     104.5
ANGLE CH1E CH2E CH1E    45.0     117.0
ANGLE CH1E CH2E CH2E    45.0     112.5
ANGLE CH1E CH2E CH3E    45.0     113.0
ANGLE CH1E CH2E OH1     45.0     111.0
ANGLE CH3E CH2E OH1     45.0     111.0
!ANGLE CH1E CH2E OH1E    45.0     111.0
!ANGLE CH3E CH2E OH1E    45.0     111.0
ANGLE CH1E CH2E S       50.0     112.5
ANGLE CH1E CH2E SH1E    50.0     112.5
ANGLE CH1E N    CH2E    60.0     110.0
ANGLE CH1E N    CH3E    60.0     110.0
ANGLE CH1E NH1  CH3E    60.0     120.0
ANGLE CH1E NH1  H       35.0     120.0
!ANGLE CH1E NH1E CH3E    65.0     120.0
ANGLE CH1E NH2  CH2E    60.0     120.0
ANGLE CH1E NH2  H       35.0     120.0
!ANGLE CH1E NH2E CH2E    65.0     120.0
ANGLE CH1E NH3  HC      35.0     109.5
ANGLE CH1E NH3  CH2E    35.0     109.5
ANGLE CH1E OH1  H       35.0     109.5
ANGLE CH2E C    CR1E    70.0     121.5
ANGLE CH2E C    N       20.0     117.5
ANGLE CH2E C    NH1     20.0     117.5
!ANGLE CH2E C    NH1E    20.0     117.5
ANGLE CH2E C    NH2     20.0     117.5
!ANGLE CH2E C    NH2E    20.0     117.5
ANGLE CH2E C    NR      60.0     116.0
ANGLE CH2E C    O       85.0     121.6
ANGLE CH2E C    OC      85.0     118.5
ANGLE CH2E C    OH1     85.0     120.0
!ANGLE CH2E C    OH1E    85.0     120.0
ANGLE CH2E CH1E CH3E    50.0     111.5
ANGLE CH2E CH1E N       65.0     104.0
ANGLE CH2E CH1E NH1     65.0     110.0
!ANGLE CH2E CH1E NH1E    65.0     110.0
ANGLE CH2E CH1E NH2     65.0     110.0
!ANGLE CH2E CH1E NH2E    65.0     110.0
ANGLE CH2E CH1E NH3     65.0     110.0
!ANGLE CH2E CH1E NH3E    65.0     110.0
ANGLE CH2E CH2E CH2E    45.0     110.0
ANGLE CH2E CH2E CH3E    45.0     111.0
ANGLE CH2E CH2E N       65.0     105.0
ANGLE CH2E CH2E NH1     65.0     111.0
!ANGLE CH2E CH2E NH1E    65.0     111.0
ANGLE CH2E CH2E NH2     65.0     109.5
!ANGLE CH2E CH2E NH2E    65.0     109.5
ANGLE CH2E CH2E NH3     65.0     110.5
!ANGLE CH2E CH2E NH3E    65.0     110.5
ANGLE CH2E CH2E S       50.0     112.5
ANGLE CH2E N    CH3E    60.0     109.5
ANGLE CH2E NH1  CH3E    60.0     120.0
ANGLE CH2E NH1  H       35.0     120.0
!ANGLE CH2E NH1E CH3E    65.0     120.0
ANGLE CH2E NH2  H       35.0     120.0
ANGLE CH2E NH3  HC      35.0     109.5
ANGLE CH2E OH1  H       35.0     109.5
ANGLE CH2E S    CH3E    50.0      99.5! FROM PARMFIX9, CHECK WITH IR
ANGLE CH2E S    S       50.0     104.2
ANGLE CH3E C    N       20.0     117.5
ANGLE CH3E C    NH1     20.0     117.5
!ANGLE CH3E C    NH1E    20.0     117.5
ANGLE CH3E C    O       85.0     121.5
ANGLE CH3E CH1E CH3E    50.0     111.0
ANGLE CH3E CH1E NH1     65.0     108.5
!ANGLE CH3E CH1E NH1E    65.0     108.5
ANGLE CH3E CH1E NH2     65.0     109.5
!ANGLE CH3E CH1E NH2E    65.0     109.5
ANGLE CH3E CH1E NH3     65.0     109.5
!ANGLE CH3E CH1E NH3E    65.0     109.5
ANGLE CH3E CH1E OH1     60.0     110.5
!ANGLE CH3E CH1E OH1E    60.0     110.5
ANGLE CH3E NH1  H       35.0     120.0
ANGLE CR1E C    CR1E    65.0     120.5
ANGLE CR1E C    NH1     65.0     110.5! USED ONLY IN HIS, NOT IT TRP
!ANGLE CR1E C    NH1E    65.0     110.5
ANGLE CR1E C    NP      65.0     122.5
ANGLE CR1E C    NR      65.0     122.5
ANGLE CR1E C    OH1     65.0     119.0
!ANGLE CR1E C    OH1E    65.0     119.0
ANGLE CR1E CR1E CR1E    90.0     120.5
ANGLE CR1E NH1  CR1E    65.0     110.0
ANGLE CR1E NH1  H       35.0     120.0
!ANGLE CR1E NH1E CR1E    80.0     110.0
ANGLE CR1E NR   CR1E    65.0     110.0
!ANGLE CR1E NR   FE      30.0     124.8!  FORCE CONSTANT FROM PARMFIX9
!ANGLE CT   C    N       20.0     117.5
!ANGLE CT   C    NH1     20.0     117.5
!ANGLE CT   C    NH2     20.0     117.5
!ANGLE CT   C    O       85.0     121.5
!ANGLE CT   C    OC      85.0     118.5
!ANGLE CT   C    OH1     85.0     120.0
!ANGLE CT   CT   CT      45.0     111.00
!ANGLE CT   CT   HA      40.0     109.50
!ANGLE CT   CT   N       65.0     105.00
!ANGLE CT   CT   NC2     65.0     110.00
!ANGLE CT   CT   NH1     65.0     110.00
!ANGLE CT   CT   NH2     65.0     110.00
!ANGLE CT   CT   NH3     65.0     110.00
!ANGLE CT   CT   OH1     50.0     109.50
!ANGLE CT   CT   S       50.0     112.50
!ANGLE CT   N    CT      60.0     110.0
!ANGLE CT   NC2  HC      35.0     120.0
!ANGLE CT   NH1  H       35.0     120.0
!ANGLE CT   NH2  H       35.0     120.0
!ANGLE CT   NH3  HC      35.0     109.5
!ANGLE CT   OH1  H       35.0     109.5
!ANGLE CT   S    CT      50.0      99.5!  FORCE CONSTANTS FROM PARMFIX9
!ANGLE CT   S    S       50.0     104.2
!ANGLE FE   CM   OM       5.0      90.0!       FROM KROEKER ET AL(1980)
!ANGLE FE   OM   OM       0.0     180.0!  DUMMY PARAMETER FOR PATCH AND ANALYSIS.
ANGLE H    NH2  H       40.0     125.0
!ANGLE HA   C    NH1     40.0     120.0
!ANGLE HA   C    NH2     40.0     120.0
!ANGLE HA   C    NR      40.0     120.0
!ANGLE HA   C    O       85.0     121.5
!ANGLE HA   CT   HA      40.0     109.5
!ANGLE HA   CT   N       50.0     109.5
!ANGLE HA   CT   NC2     50.0     109.5
!ANGLE HA   CT   NH1     50.0     109.5
!ANGLE HA   CT   NH3     50.0     109.5
!ANGLE HA   CT   OH1     50.0     109.5
!ANGLE HA   CT   S       40.0     109.5
ANGLE HC   NC2  HC      40.0     120.0
ANGLE HC   NH3  HC      40.0     109.5
ANGLE N    C    O       85.0     121.0
ANGLE NC2  C    NC2     70.0     120.0
ANGLE NC2  C    NH1     70.0     120.0
!ANGLE NC2E C    NC2E    70.0     120.0
!ANGLE NC2E C    NH1E    70.0     120.0
ANGLE NH1  C    NR      70.0     120.0
ANGLE NH1  C    O       65.0     121.0
ANGLE NH1  CR1E NH1     70.0     109.0
ANGLE NH1  CR1E NR      70.0     109.0
!ANGLE NH1E C    O       65.0     121.0
!ANGLE NH1E CR1E NH1E    70.0     109.0
!ANGLE NH1E CR1E NR      70.0     109.0
ANGLE NH2  C    O       65.0     121.0
!ANGLE NH2E C    O       65.0     121.0
!ANGLE NP   FE   CM      5.0       90.0
!ANGLE NP   FE   NP      50.0      90.0
!ANGLE NP   FE   NR      5.0      115.0
!ANGLE NP   FE   OM      5.0       90.0! JUST A GUESS FROM EXISTING FE CM DATA
!ANGLE NR   FE   CM      5.0      180.0
!ANGLE NR   FE   OM      5.0      180.0! JUST A GUESS FROM EXISTING FE CM DATA
ANGLE O    C    OH1     85.0     120.0
!ANGLE O    C    OH1E    85.0     120.0
ANGLE OC   C    OC      85.0     122.5
ANGLE OC   S    OC      85.0     109.5! FORCE CONSTANT JST A GUESS.


!                    Ed = Kd ( 1.0 + cos( n*phi + d )
!
!                            Kd         n       d     ,
!                     (kcal/mol A**2)          deg.
!
!                      where d can be 0.0 or 180.0
!
DIHE CH1E C    N    CH1E    10.0       2     180.0! PRO ISOM. BARRIER 20 KCAL/MOL.
DIHE CH2E C    N    CH1E    10.0       2     180.0
DIHE CR1E C    C    CR1E     5.0       2     180.0! => TRP OOP. VIB 170CM 1
DIHE CR1E C    C    C        2.5       2     180.0! SEE BEHLEN ET AL JCP 75:5685 81
DIHE CR1E C    C    NH1      2.5       2     180.0
DIHE CR1E C    C    NH1E     2.5       2     180.0
DIHE X    C    CH1E X        0.0       3       0.0! FROM GELIN THESIS AMIDES
DIHE X    C    CH2E X        0.0       3       0.0! USING A SINGLE
DIHE X    C    CR1E X       10.0       2     180.0! DIHEDRAL PER BOND RATHER
!DIHE X    C    CT   X        0.0       3       0.0! THAN MULTIPLE TORSIONS.
DIHE X    C    N    X        8.2       2     180.0! ALKANE TORSION REDUCED TO
DIHE X    C    NC2  X        8.2       2     180.0! 1.6 FROM 1.8 TO COINCIDE WITH
DIHE X    C    NH1  X        8.2       2     180.0! THE EXPERIMENTAL BARRIER.
!DIHE X    C    NH1E X        8.2       2     180.0! IRON AND SULFUR TORSIONS
DIHE X    C    NH2  X        8.2       2     180.0
DIHE X    C    OH1  X        1.8       2     180.0
DIHE X    CH1E CH1E X        1.6       3       0.0
DIHE X    CH1E CH2E X        1.6       3       0.0
DIHE X    CH1E N    X        0.3       3       0.0! FROM HAGLER ET AL TABULATION OF
DIHE X    CH1E NH1  X        0.3       3       0.0! EXP. DATA AND 6 31G CALC.
!DIHE X    CH1E NH1E X        0.3       3       0.0
DIHE X    CH1E NH2  X        1.8       3       0.0! PROTONATED SECONDARY AMINE
DIHE X    CH1E NH3  X        0.6       3       0.0! 1/PROTON SO 3 FOR THE BOND
DIHE X    CH1E OH1  X        0.5       3       0.0! CHANGED TO ROUGHLY MEOH
DIHE X    CH2E CH2E X        1.6       3       0.0
DIHE X    CH2E N    X        0.3       3       0.0! SEE CH1E COMMENTS
DIHE X    CH2E NH1  X        0.3       3       0.0
!DIHE X    CH2E NH1E X        0.3       3       0.0
DIHE X    CH2E NH2  X        0.6       3       0.0
DIHE X    CH2E NH3  X        0.6       3       0.0
DIHE X    CH2E OH1  X        0.5       3       0.0
DIHE X    CH2E S    X        1.2       2       0.0
!DIHE X    CT   CT   X        1.6       3       0.0
!DIHE X    CT   N    X        0.3       3       0.0! SEE CH1E COMMENTS
!DIHE X    CT   NC2  X        0.3       3       0.0
!DIHE X    CT   NH1  X        0.3       3       0.0
!DIHE X    CT   NH2  X        0.6       3       0.0
!DIHE X    CT   NH3  X        0.6       3       0.0
!DIHE X    CT   OH1  X        0.5       3       0.0
!DIHE X    CT   S    X        1.2       2       0.0
!DIHE X    FE   NR   X        0.05      4       0.0
!DIHE X    FE   CM   X        0.05      4       0.0
!DIHE X    FE   OM   X        0.00      4       0.0
DIHE X    S    S    X        4.0       2       0.0! FROM EXP. NMR BARRIER

!                     Ei = Ki ( omega - omega0 )**2  (for periodicity=0)
!
!                               Ki    periodicity  omega0
!                         (Kcal/mol A**2)           deg.
!
IMPROPER C    C    CR1E CH2E    90.0    0  0.0! GIVE 220 CM 1 METHYL OOP FOR TOLUENE.
IMPROPER C    CR1E C    CH2E    90.0    0  0.0! USED HERE FOR TRP CG OUT OF PLANE
IMPROPER C    CR1E CR1E CH2E    90.0    0  0.0!               PHE, AND TYR CG OOP
IMPROPER C    CR1E NH1  CH2E    90.0    0  0.0!               HIS CG RING OOP
IMPROPER C    NH1  CR1E CH2E    90.0    0  0.0!
IMPROPER C    CR1E CR1E OH1    150.0    0  0.0! GIVE 249 CM 1 PHENOL OH OOP.
!IMPROPER C    CR1E CR1E OH1E   150.0    0  0.0! USED HERE FOR TYR HYROXYL OOP
IMPROPER C    H    H    NH2     45.0    0  0.0! PRIMARY AMIDES (ASN AND GLN) OOP
IMPROPER C    OC   OC   CH1E   100.0    0  0.0! CARBOXYL OUT OF PLANE.
IMPROPER C    OC   OC   CH2E   100.0    0  0.0!
IMPROPER C    X    X    C       25.0    0  0.0! FROM BENZENE NORMAL MODE ANALYSIS
IMPROPER C    X    X    CH2E    90.0    0  0.0! FROM TOLUENE METHYL OOP. 217 CM 1
IMPROPER C    X    X    CH3E    90.0    0  0.0
IMPROPER C    X    X    CR1E    25.0    0  0.0
IMPROPER C    X    X    H       75.0    0  0.0! FROM BENZENE NORMAL MODE ANALYSIS
!IMPROPER C    X    X    HA      75.0    0  0.0!
IMPROPER C    X    X    NH1    100.0    0  0.0! AMIDES FIT TO N METHYL ACETAMIDE.
!IMPROPER C    X    X    NH1E   100.0    0  0.0
IMPROPER C    X    X    O      100.0    0  0.0
IMPROPER C    X    X    OC     100.0    0  0.0
IMPROPER C    X    X    OH1    150.0    0  0.0! USED FOR TYR HYDROXYL OOP
!IMPROPER C    X    X    OH1E   150.0    0   0.0
IMPROPER CH1E X    X    CH1E    55.0  0  35.26439! CALCULATED TO  BE THE SAME AS THE 3
IMPROPER CH1E X    X    CH2E    55.0  0  35.26439! H CH1E X ANGLES WITH K=40
IMPROPER CH1E X    X    CH3E    55.0  0  35.26439
IMPROPER CR1E X    X    CR1E    25.0  0    0.0! EXTENDED ATOM VERSION OF BENZENE
IMPROPER CR1E X    X    NH1     25.0  0    0.0! SAME AS ABOVE FOR LACK OF BETTER VALUES
!IMPROPER CR1E X    X    NH1E    25.0  0    0.0
!IMPROPER FE   X    X    NP      20.0  0    0.0! FROM PARMFIX9
IMPROPER H    X    X    O       45.0  0    0.0
IMPROPER N    CH1E CH2E C       45.0  0    0.0! PROLINE NITROGENS
IMPROPER N    X    X    CH2E    45.0  0    0.0
!IMPROPER N    X    X    CT      45.0  0    0.0
!IMPROPER NC2  X    X    CT      45.0  0    0.0
IMPROPER NC2  X    X    HC      45.0  0    0.0
IMPROPER NH1  X    X    CH1E    45.0  0    0.0
IMPROPER NH1  X    X    CH2E    45.0  0    0.0
IMPROPER NH1  X    X    CH3E    45.0  0    0.0
!IMPROPER NH1  X    X    CT      45.0  0    0.0
IMPROPER NH1  X    X    H       45.0   0   0.0! AMIDES PROTON OOP
IMPROPER NH1  X    X    NH1     25.0   0   0.0! 
IMPROPER NH1  X    X    NR      25.0   0   0.0
!IMPROPER NH1E X    X    NH1E    25.0   0   0.0
!IMPROPER NH1E X    X    NR      25.0   0   0.0
IMPROPER NH2  X    X    H       45.0   0   0.0
IMPROPER NR   X    X    C       25.0   0   0.0
IMPROPER NR   X    X    CR1E    25.0   0   0.0
!IMPROPER NR   X    X    CT      25.0   0   0.0


!                   Lennard-Jones parameters
!                                       ------1-4-------
!                  epsilon   sigma       epsilon   sigma
!                 (Kcal/mol)   (A)      (Kcal/mol) (A)
!
NONBONDED C         0.1200    3.7418      0.1   3.3854 !carbonyl carbon
NONBONDED CH1E      0.0486    4.2140      0.1   3.3854 !\
NONBONDED CH2E      0.1142    3.9823      0.1   3.3854 ! extended carbons
NONBONDED CH3E      0.1811    3.8576      0.1   3.3854 !/
NONBONDED CR1E      0.1200    3.7418      0.1   3.3854 ! ring carbons
!NONBONDED CT        0.0262    4.4367      0.1   3.3854 ! explicit H carbon
!NONBONDED CM        0.0262    4.4367      0.1   3.3854 !for carbon monoxide
!
!NONBONDED FE        0.0000    1.1582      0.0000 1.1582
!
NONBONDED H         0.0498    1.4254      0.0498 1.4254
!NONBONDED HA        0.0450    2.6157      0.0450 2.6157
NONBONDED HC        0.0498    1.4254      0.0498 1.4254
!
NONBONDED N         0.2384    2.8509      0.2384  2.8509
NONBONDED NC2       0.2384    2.8509      0.2384  2.8509
NONBONDED NH1       0.2384    2.8509      0.2384  2.8509  
NONBONDED NH2       0.2384    2.8509      0.2384  2.8509  
NONBONDED NH3       0.2384    2.8509      0.2384  2.8509
NONBONDED NP        0.2384    2.8509      0.2384  2.8509  
NONBONDED NR        0.2384    2.8509      0.2384  2.8509  
!NONBONDED NC2E      0.3676    3.029       0.3676  3.029   !extended
!NONBONDED NH1E      0.3074    2.940       0.3074  1.6500  !extended
!NONBONDED NH2E      0.3676    3.029       0.3676  3.029   !extended
!NONBONDED NH3E      0.4596    3.1181      0.4596  3.1181  !extended
!
NONBONDED O         0.1591   2.8509       0.1591  2.8509  
NONBONDED OC        0.6469   2.8509       0.6469  2.8509  
NONBONDED OH1       0.1591   2.8509       0.1591  2.8509  
!NONBONDED OM        0.1591   2.8509       0.1591  2.8509  !carbon monoxide
!NONBONDED OH1E      0.2440   2.940        0.2440  2.940   !extended
!NONBONDED OH2E      0.2039   3.029        0.2039  3.029   !extended
!
NONBONDED S         0.0430   3.3676       0.0430  3.3676
NONBONDED SH1E      0.0430   3.3676       0.0430  3.3676
!


! special solute-solute hydrogen bonding potential parameters
AEXP 4
REXP 6 
HAEX 4 
AAEX 2
!                   Emin      Rmin
!                (Kcal/mol)   (A)
 {* note the wildcard specifications *}
HBOND N*   N*      -3.00      3.0 !  VALUES FROM VINOGRADOV AND LINELL FOR
HBOND N*   O*      -3.50      2.9 !  TYPICAL LENGTHS AND DEPTHS.
HBOND O*   N*      -4.00      2.85
HBOND O*   O*      -4.25      2.75


! the following NBFIXes are for hydrogen - acceptor v.d. Waals
! interactions where explicit hydrogen bonding terms are present
!                                         ------1-4------
!                   A           B           A          B
!           [Kcal/(mol A^12)] [Kcal/(mol A^6)]
!
NBFIX H    NP      0.05       0.1      0.05       0.1
NBFIX H    NR      0.05       0.1      0.05       0.1
NBFIX H    O       0.05       0.1      0.05       0.1
NBFIX H    OC      0.05       0.1      0.05       0.1
NBFIX H    OH1     0.05       0.1      0.05       0.1

NBFIX HC   NP      0.05       0.1      0.05       0.1
NBFIX HC   NR      0.05       0.1      0.05       0.1
NBFIX HC   O       0.05       0.1      0.05       0.1
NBFIX HC   OC      0.05       0.1      0.05       0.1
NBFIX HC   OH1     0.05       0.1      0.05       0.1

set echo=true end
