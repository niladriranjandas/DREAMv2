 REMARKS Charmm topology for proteins v22 b5
 REMARKS FILENAME="topallh22x-eefx.pro"  # modified for eefx/immx calculations. 
 
 SET ECHO=FALSE END 

eval ($protein_top_vers="charmm22-eef-1.0")

 {>>>>>>>> Developmental Residue Topology File for Proteins <<<<<
  >>>>>>>>>>>>>>>>> Using All Hydrogens (ALLH) <<<<<<<<<<<<<<<<<<
  >>>>>>>>>>>>>>>>>>>>>>> May 1993 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  >>>>>>> Direct comments to Alexander D. MacKerell Jr. <<<<<<<<<
  >>>>>> 410-706-7442 or bitnet: alex,mmiris.ab.umd.edu <<<<<<<<<
  These files are a beta release; additional parameter development
  and testing may lead to alteration of the contents.}
 
 AUTOGENERATE ANGLES=TRUE  DIHEDRALS=TRUE  END
 
 MASS H      1.0080 ! polar H 
 MASS HC     1.0080 ! N-ter H 
 MASS HA     1.0080 ! nonpolar H 
 MASS HT     1.0080 ! TIPS3P WATER HYDROGEN 
 MASS HP     1.0080 ! aromatic H 
 MASS HB     1.0080 ! backbone H 
 MASS HR1    1.0080 ! his he1, (+) his HG,HD2 
 MASS HR2    1.0080 ! (+) his HE1 
 MASS HR3    1.0080 ! neutral his HG, HD2 
 MASS HS     1.00800 ! thiol hydrogen
 MASS C     12.0110 ! carbonyl C, peptide backbone 
 MASS CA    12.0110 ! aromatic C 
 MASS CT1   12.0110 ! aliphatic sp3 C for CH 
 MASS CT2   12.0110 ! aliphatic sp3 C for CH2 
 MASS CT3   12.0110 ! aliphatic sp3 C for CH3 
 MASS CPH1  12.0110 ! his CG and CD2 carbons 
 MASS CPH2  12.0110 ! his CE1 carbon 
 MASS CPT   12.0110 ! trp C between rings 
 MASS CY    12.0110 ! TRP C in pyrrole ring 
 MASS CP1   12.0110 ! tetrahedral C (proline CA) 
 MASS CP2   12.0110 ! tetrahedral C (proline CB/CG) 
 MASS CP3   12.0110 ! tetrahedral C (proline CD) 
 MASS CC    12.0110 ! carbonyl C, asn,asp,gln,glu,cter,ct2 
 MASS CD    12.0110 ! carbonyl C, pres aspp,glup,ct1 
 MASS CPA   12.0110 ! heme alpha-C 
 MASS CPB   12.0110 ! heme beta-C 
 MASS CPM   12.0110 ! heme meso-C 
 MASS CM    12.0110 ! heme CO carbon 
 MASS CS    12.0110 ! thiolate carbon 
 MASS N     14.0070 ! proline N 
 MASS NR1   14.0070 ! neutral his protonated ring nitrogen 
 MASS NR2   14.0070 ! neutral his unprotonated ring nitrogen 
 MASS NR3   14.0070 ! charged his ring nitrogen 
 MASS NH1   14.0070 ! peptide nitrogen 
 MASS NH2   14.0070 ! amide nitrogen 
 MASS NH3   14.0070 ! ammonium nitrogen 
 MASS NC2   14.0070 ! guanidinium nitroogen 
 MASS NY    14.0070 ! TRP N in pyrrole ring 
 MASS NP    14.0070 ! Proline ring NH2+ (N-terminal) 
 MASS NPH   14.0070 ! heme pyrrole N 
 MASS O     15.9990 ! carbonyl oxygen 
 MASS OB    15.9990 ! carbonyl oxygen in acetic acid 
 MASS OC    15.9990 ! carboxylate oxygen 
 MASS OH1   15.9990 ! hydroxyl oxygen 
 MASS OS    15.9994 ! ester oxygen 
 MASS OT    15.9994 ! TIPS3P WATER OXYGEN 
 MASS OM    15.9990 ! heme CO/O2 oxygen 
 MASS S     32.0600 ! sulphur 
 MASS SM    32.0600 ! sulfur C-S-S-C type 
 MASS SS    32.0600 ! thiolate sulfur 
 MASS CAL   40.0800 ! calcium 2+ 
 MASS ZN    65.3700 ! zinc (II) cation 
 MASS FE    55.8470 ! heme iron 56 
 MASS CLA   34.0000 ! CHLORIDE Anion 
 MASS SOD   22.9898 ! Sodium Ion 
 !-----------------------------------------------------------
 
 RESIDUE ALA   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     | 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |     HB1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |    / 
 GROUP						    !  HA-CA--CB-HB2 
       ATOM  CB    TYPE=CT3  CHARGE=    -.2700  END !     |    \ 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     |     HB3 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   O=C 
       ATOM  HB3   TYPE=HA   CHARGE=     .0900  END !     | 
 GROUP						    ! 
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CB    HB3    
 
      DONOR  HN    N      
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4592  114.44  123.23  111.09   1.5461  
      IC  N     C      *CA   HA       1.4592  114.44 -120.45  106.39   1.0840  
      IC  C     CA     CB    HB1      1.5390  111.09  177.25  109.60   1.1109  
      IC  HB1   CA     *CB   HB2      1.1109  109.60  119.13  111.05   1.1119  
      IC  HB1   CA     *CB   HB3      1.1109  109.60 -119.58  111.61   1.1114  
 
 
 END {ALA }
 !-----------------------------------------------------------
 
 RESIDUE ARG   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     |                      HH11 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N                       | 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1 HG1 HD1 HE     NH1-HH12 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |   |   |   |     / 
 GROUP						    !  HA-CA--CB--CG--CD--NE--CZ 
       ATOM  CB    TYPE=CT2  CHARGE=    -.1800  END !     |   |   |   |         \ 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     |   HB2 HG2 HD2        NH2-HH22 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   O=C                       | 
 GROUP						    !     |                      HH21 
       ATOM  CG    TYPE=CT2  CHARGE=    -.1800  END  
       ATOM  HG1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG2   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CD    TYPE=CT2  CHARGE=    -.3000  END  
       ATOM  HD1   TYPE=HA   CHARGE=     .0500  END  
       ATOM  HD2   TYPE=HA   CHARGE=     .0500  END  
       ATOM  NE    TYPE=NC2  CHARGE=    -.2800  END  
       ATOM  HE    TYPE=HC   CHARGE=     .1200  END  
       ATOM  CZ    TYPE=C    CHARGE=    -.2000  END  
       ATOM  NH1   TYPE=NC2  CHARGE=    -.1210  END  
       ATOM  HH11  TYPE=HC   CHARGE=     .2005  END  
       ATOM  HH12  TYPE=HC   CHARGE=     .2005  END  
       ATOM  NH2   TYPE=NC2  CHARGE=    -.1210  END  
       ATOM  HH21  TYPE=HC   CHARGE=     .2005  END  
       ATOM  HH22  TYPE=HC   CHARGE=     .2005  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  CD    CG     
      BOND  NE    CD     
      BOND  CZ    NE     
      BOND  NH1   CZ     
      BOND  NH2   CZ     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CG    HG1    
      BOND  CG    HG2    
      BOND  CD    HD1    
      BOND  CD    HD2    
      BOND  NE    HE     
      BOND  NH1   HH11   
      BOND  NH1   HH12   
      BOND  NH2   HH21   
      BOND  NH2   HH22   
 
      IMPROPER  CZ    NH1   NH2   NE     
 
      DONOR  HN    N      
      DONOR  HE    NE     
      DONOR  HH11  NH1    
      DONOR  HH12  NH1    
      DONOR  HH21  NH2    
      DONOR  HH22  NH2    
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4544  109.86  123.64  112.26   1.5552  
      IC  N     C      *CA   HA       1.4544  109.86 -117.93  106.61   1.0836  
      IC  N     CA     CB    CG       1.4544  110.70  180.00  115.95   1.5475  
      IC  CG    CA     *CB   HB1      1.5475  115.95  120.05  106.40   1.1163  
      IC  CG    CA     *CB   HB2      1.5475  115.95 -125.81  109.55   1.1124  
      IC  CA    CB     CG    CD       1.5552  115.95  180.00  114.01   1.5384  
      IC  CD    CB     *CG   HG1      1.5384  114.01  125.20  108.55   1.1121  
      IC  CD    CB     *CG   HG2      1.5384  114.01 -120.30  108.96   1.1143  
      IC  CB    CG     CD    NE       1.5475  114.01  180.00  107.09   1.5034  
      IC  NE    CG     *CD   HD1      1.5034  107.09  120.69  109.41   1.1143  
      IC  NE    CG     *CD   HD2      1.5034  107.09 -119.04  111.52   1.1150  
      IC  CG    CD     NE    CZ       1.5384  107.09  180.00  123.05   1.3401  
      IC  CZ    CD     *NE   HE       1.3401  123.05  180.00  113.14   1.0065  
      IC  CD    NE     CZ    NH1      1.5034  123.05  180.00  118.06   1.3311  
      IC  NE    CZ     NH1   HH11     1.3401  118.06 -178.28  120.61    .9903  
      IC  HH11  CZ     *NH1  HH12      .9903  120.61  171.19  116.29   1.0023  
      IC  NH1   NE     *CZ   NH2      1.3311  118.06  178.64  122.14   1.3292  
      IC  NE    CZ     NH2   HH21     1.3401  122.14 -174.14  119.91    .9899  
      IC  HH21  CZ     *NH2  HH22      .9899  119.91  166.16  116.88    .9914  
 
 
 END {ARG }
 !-----------------------------------------------------------
 
 RESIDUE ASN   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     | 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1 OD1    HD21 (cis to OD1) 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |   ||    / 
 GROUP						    !  HA-CA--CB--CG--ND2 
       ATOM  CB    TYPE=CT2  CHARGE=    -.1800  END !     |   |         \ 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     |   HB2        HD22 (trans to OD1) 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   O=C 
 GROUP						    !     | 
       ATOM  CG    TYPE=CC   CHARGE=     .5500  END  
       ATOM  OD1   TYPE=O    CHARGE=    -.5500  END  
 GROUP  
       ATOM  ND2   TYPE=NH2  CHARGE=    -.5000  END  
       ATOM  HD21  TYPE=H    CHARGE=     .2500  END  
       ATOM  HD22  TYPE=H    CHARGE=     .2500  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  OD1   CG     
      BOND  ND2   CG     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  ND2   HD21   
      BOND  ND2   HD22   
 
      IMPROPER  CG    ND2   CB    OD1    
      IMPROPER  CG    CB    ND2   OD1    
      IMPROPER  ND2   CG    HD21  HD22   
      IMPROPER  ND2   CG    HD22  HD21   
 
      DONOR  HN    N      
      DONOR  HD21  ND2    
      DONOR  HD22  ND2    
 
      ACCEPTOR  OD1   CG     
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4510  105.23  121.18  113.04   1.5627  
      IC  N     C      *CA   HA       1.4510  105.23 -115.52  107.63   1.0848  
      IC  N     CA     CB    CG       1.4510  110.91  180.00  114.30   1.5319  
      IC  CG    CA     *CB   HB1      1.5319  114.30  119.17  107.82   1.1120  
      IC  CG    CA     *CB   HB2      1.5319  114.30 -123.74  110.34   1.1091  
      IC  CA    CB     CG    OD1      1.5627  114.30  180.00  122.56   1.2323  
      IC  OD1   CB     *CG   ND2      1.2323  122.56 -179.19  116.15   1.3521  
      IC  CB    CG     ND2   HD21     1.5319  116.15 -179.26  117.35    .9963  
      IC  HD21  CG     *ND2  HD22      .9963  117.35  178.02  120.05    .9951  
 
 
 END {ASN }
 !-----------------------------------------------------------
 
 RESIDUE ASP   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     | 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1   OD1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |    // 
 GROUP						    !  HA-CA--CB--CG 
       ATOM  CB    TYPE=CT2  CHARGE=    -.2800  END !     |   |    \\ 
       ATOM  HB1   TYPE=HA   CHARGE=     .1400  END !     |   HB2   OD2 
       ATOM  HB2   TYPE=HA   CHARGE=     .1400  END !   O=C 
       ATOM  CG    TYPE=CC   CHARGE=    1.0000  END !     | 
       ATOM  OD1   TYPE=OC   CHARGE=    -.5000  END  
       ATOM  OD2   TYPE=OC   CHARGE=    -.5000  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  OD1   CG     
      BOND  OD2   CG     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
 
      IMPROPER  OD1   CB    OD2   CG     
 
      DONOR  HN    N      
 
      ACCEPTOR  OD1   CG     
      ACCEPTOR  OD2   CG     
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4490  105.63  122.33  114.10   1.5619  
      IC  N     C      *CA   HA       1.4490  105.63 -116.40  106.77   1.0841  
      IC  N     CA     CB    CG       1.4490  111.10  180.00  112.60   1.5218  
      IC  CG    CA     *CB   HB1      1.5218  112.60  119.22  109.23   1.1086  
      IC  CG    CA     *CB   HB2      1.5218  112.60 -121.61  110.64   1.1080  
      IC  CA    CB     CG    OD1      1.5619  112.60  180.00  117.99   1.2565  
      IC  OD1   CB     *CG   OD2      1.2565  117.99 -170.23  117.70   1.2541  
 
 
 END {ASP }
 !-----------------------------------------------------------
 
 RESIDUE CYS   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     | 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   | 
 GROUP						    !  HA-CA--CB--SG 
       ATOM  CB    TYPE=CT2  CHARGE=    -.1100  END !     |   |     \ 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     |   HB2    HG1 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   O=C 
       ATOM  SG    TYPE=S    CHARGE=    -.2300  END !     | 
       ATOM  HG    TYPE=HS   CHARGE=     .1600  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  SG    CB     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  SG    HG    

      DIHEDRAL  CA  CB  SG  HG  ! multiple dihedral
      DIHEDRAL  CA  CB  SG  HG  ! multiple dihedral
      DIHEDRAL  CA  CB  SG  HG  ! multiple dihedral

      DONOR  HN    N      
      DONOR  HG   SG     
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4533  105.89  121.79  111.98   1.5584  
      IC  N     C      *CA   HA       1.4533  105.89 -116.34  107.71   1.0837  
      IC  N     CA     CB    SG       1.4533  111.56  180.00  113.87   1.8359  
      IC  SG    CA     *CB   HB1      1.8359  113.87  119.91  107.24   1.1134  
      IC  SG    CA     *CB   HB2      1.8359  113.87 -125.32  109.82   1.1124  
      IC  CA    CB     SG    HG       1.5584  113.87  176.96   97.15   1.3341  
 
 
 END {CYS }
 !-----------------------------------------------------------
 
 RESIDUE GLN   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     | 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1 HG1 OE1   HE21 (cis to OE1) 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |   |   ||    / 
 GROUP						    !  HA-CA--CB--CG--CD--NE2 
       ATOM  CB    TYPE=CT2  CHARGE=    -.1800  END !     |   |   |         \ 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     |   HB2 HG2       HE22 (trans to OE1) 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   O=C 
 GROUP						    !     | 
       ATOM  CG    TYPE=CT2  CHARGE=    -.1800  END  
       ATOM  HG1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG2   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CD    TYPE=CC   CHARGE=     .5500  END  
       ATOM  OE1   TYPE=O    CHARGE=    -.5500  END  
 GROUP  
       ATOM  NE2   TYPE=NH2  CHARGE=    -.5000  END  
       ATOM  HE21  TYPE=H    CHARGE=     .2500  END  
       ATOM  HE22  TYPE=H    CHARGE=     .2500  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  CD    CG     
      BOND  OE1   CD     
      BOND  NE2   CD     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CG    HG1    
      BOND  CG    HG2    
      BOND  NE2   HE21   
      BOND  NE2   HE22   
 
      IMPROPER  CD    NE2   CG    OE1    
      IMPROPER  CD    CG    NE2   OE1    
      IMPROPER  NE2   CD    HE21  HE22   
      IMPROPER  NE2   CD    HE22  HE21   
 
      DONOR  HN    N      
      DONOR  HE21  NE2    
      DONOR  HE22  NE2    
 
      ACCEPTOR  OE1   CD     
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4506  106.57  121.91  111.68   1.5538  
      IC  N     C      *CA   HA       1.4506  106.57 -116.82  107.53   1.0832  
      IC  N     CA     CB    CG       1.4506  111.44  180.00  115.52   1.5534  
      IC  CG    CA     *CB   HB1      1.5534  115.52  120.93  106.80   1.1147  
      IC  CG    CA     *CB   HB2      1.5534  115.52 -124.58  109.34   1.1140  
      IC  CA    CB     CG    CD       1.5538  115.52  180.00  112.50   1.5320  
      IC  CD    CB     *CG   HG1      1.5320  112.50  118.69  110.41   1.1112  
      IC  CD    CB     *CG   HG2      1.5320  112.50 -121.91  110.74   1.1094  
      IC  CB    CG     CD    OE1      1.5534  112.50  180.00  121.52   1.2294  
      IC  OE1   CG     *CD   NE2      1.2294  121.52  179.57  116.84   1.3530  
      IC  CG    CD     NE2   HE21     1.5320  116.84 -179.72  116.86    .9959  
      IC  HE21  CD     *NE2  HE22      .9959  116.86 -178.91  119.83    .9943  
 
 
 END {GLN }
 !-----------------------------------------------------------
 
 RESIDUE GLU   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     | 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1 HG1   OE1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |   |    // 
 GROUP						    !  HA-CA--CB--CG--CD 
       ATOM  CB    TYPE=CT2  CHARGE=    -.1800  END !     |   |   |    \\ 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     |   HB2 HG2   OE2 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   O=C 
 GROUP						    !     | 
       ATOM  CG    TYPE=CT2  CHARGE=    -.2800  END  
       ATOM  HG1   TYPE=HA   CHARGE=     .1400  END  
       ATOM  HG2   TYPE=HA   CHARGE=     .1400  END  
       ATOM  CD    TYPE=CC   CHARGE=    1.0000  END  
       ATOM  OE1   TYPE=OC   CHARGE=    -.5000  END  
       ATOM  OE2   TYPE=OC   CHARGE=    -.5000  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  CD    CG     
      BOND  OE1   CD     
      BOND  OE2   CD     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CG    HG1    
      BOND  CG    HG2    
 
      IMPROPER  OE1   CG    OE2   CD     
 
      DONOR  HN    N      
 
      ACCEPTOR  OE1   CD     
      ACCEPTOR  OE2   CD     
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4512  107.27  121.90  111.71   1.5516  
      IC  N     C      *CA   HA       1.4512  107.27 -118.06  107.26   1.0828  
      IC  N     CA     CB    CG       1.4512  111.04  180.00  115.69   1.5557  
      IC  CG    CA     *CB   HB1      1.5557  115.69  121.22  108.16   1.1145  
      IC  CG    CA     *CB   HB2      1.5557  115.69 -123.65  109.81   1.1131  
      IC  CA    CB     CG    CD       1.5516  115.69  180.00  115.73   1.5307  
      IC  CD    CB     *CG   HG1      1.5307  115.73  117.38  109.50   1.1053  
      IC  CD    CB     *CG   HG2      1.5307  115.73 -121.96  111.00   1.1081  
      IC  CB    CG     CD    OE1      1.5557  115.73  180.00  114.99   1.2590  
      IC  OE1   CG     *CD   OE2      1.2590  114.99 -179.10  120.08   1.2532  
 
 
 END {GLU }
 !-----------------------------------------------------------
 
 RESIDUE GLY   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     | 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !     N-H 
       ATOM  CA    TYPE=CT2  CHARGE=    -.0200  END !     | 
       ATOM  HA1   TYPE=HB   CHARGE=     .0900  END !     | 
       ATOM  HA2   TYPE=HB   CHARGE=     .0900  END ! HA1-CA-HA2 
 GROUP						    !     | 
       ATOM  C     TYPE=C    CHARGE=     .5100  END !     | 
       ATOM  O     TYPE=O    CHARGE=    -.5100  END !     C=O 
      !     | 
 !END GROUP
 
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA1    
      BOND  CA    HA2    
 
      DONOR  HN    N      
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   HA1      1.4553  108.94  117.86  108.03   1.0814  
      IC  N     C      *CA   HA2      1.4553  108.94 -118.12  107.95   1.0817  
 
 
 END {GLY }
 !-----------------------------------------------------------
 
 RESIDUE HSD  ! neutral HIS, proton on ND1 
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     |          HD1    HE1 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N           |     / 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1    ND1--CE1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |     /      | 
 GROUP						    !  HA-CA--CB--CG       | 
       ATOM  ND1   TYPE=NR1  CHARGE=    -.3600  END !     |   |     \      | 
       ATOM  HD1   TYPE=H    CHARGE=     .3200  END !     |   HB2    CD2--NE2 
       ATOM  CG    TYPE=CPH1 CHARGE=    -.0500  END !   O=C           | 
       ATOM  CB    TYPE=CT2  CHARGE=    -.0900  END !     |          HD2 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  NE2   TYPE=NR2  CHARGE=    -.3000  END  
       ATOM  CD2   TYPE=CPH1 CHARGE=     .0500  END  
       ATOM  HD2   TYPE=HR3  CHARGE=     .0500  END  
       ATOM  CE1   TYPE=CPH2 CHARGE=     .1500  END  
       ATOM  HE1   TYPE=HR1  CHARGE=     .0500  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  ND1   CG     
      BOND  CD2   CG     
      BOND  CE1   ND1    
      BOND  NE2   CD2    
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  NE2   CE1    
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  ND1   HD1    
      BOND  CD2   HD2    
      BOND  CE1   HE1    

      DIHEDRAL  CA    CB    CG    CD2  ! multiple dihedral
      DIHEDRAL  CA    CB    CG    CD2  ! multiple dihedral
      DIHEDRAL  CA    CB    CG    CD2  ! multiple dihedral
 
      IMPROPER  ND1   CG    CE1   HD1    
      IMPROPER  CD2   CG    NE2   HD2    
      IMPROPER  CE1   ND1   NE2   HE1    
      IMPROPER  ND1   CE1   CG    HD1    
      IMPROPER  CD2   NE2   CG    HD2    
      IMPROPER  CE1   NE2   ND1   HE1    
 
      DONOR  HN    N      
      DONOR  HD1   ND1    
 
      ACCEPTOR  NE2   NONE   
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4521  107.70  122.46  109.99   1.5519  
      IC  N     C      *CA   HA       1.4521  107.70 -117.49  107.37   1.0830  
      IC  N     CA     CB    CG       1.4521  112.12  180.00  114.05   1.5041  
      IC  CG    CA     *CB   HB1      1.5041  114.05  121.17  109.01   1.1118  
      IC  CG    CA     *CB   HB2      1.5041  114.05 -122.36  109.53   1.1121  
      IC  CA    CB     CG    ND1      1.5519  114.05   90.00  124.10   1.3783  
      IC  ND1   CB     *CG   CD2      1.3783  124.10 -171.29  129.60   1.3597  
      IC  CB    CG     ND1   CE1      1.5041  124.10 -173.21  107.03   1.3549  
      IC  CB    CG     CD2   NE2      1.5041  129.60  171.99  110.03   1.3817  
      IC  NE2   ND1    *CE1  HE1      1.3166  111.63 -179.63  123.89   1.0932  
      IC  CE1   CG     *ND1  HD1      1.3549  107.03 -174.65  126.26   1.0005  
      IC  NE2   CG     *CD2  HD2      1.3817  110.03 -177.85  129.63   1.0834  
 
 
 END {HSD }
 !-----------------------------------------------------------
 
 RESIDUE HSE  ! neutral His, proton on NE2 
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     |                 HE1 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N                 / 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1    ND1--CE1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |     /      | 
 GROUP						    !  HA-CA--CB--CG       | 
       ATOM  NE2   TYPE=NR1  CHARGE=    -.3000  END !     |   |     \      | 
       ATOM  HE2   TYPE=H    CHARGE=     .2000  END !     |   HB2    CD2--NE2 
       ATOM  CD2   TYPE=CPH1 CHARGE=     .0500  END !   O=C           |     \ 
       ATOM  HD2   TYPE=HR3  CHARGE=     .0500  END !     |          HD2    HE2 
 GROUP  
       ATOM  ND1   TYPE=NR2  CHARGE=    -.1000  END  
       ATOM  CG    TYPE=CPH1 CHARGE=     .1000  END  
       ATOM  CE1   TYPE=CPH2 CHARGE=     .0000  END  
       ATOM  HE1   TYPE=HR1  CHARGE=     .0000  END  
 GROUP
       ATOM  CB    TYPE=CT2  CHARGE=    -.0800  END  
       ATOM  HB1   TYPE=HA   CHARGE=     .0400  END  
       ATOM  HB2   TYPE=HA   CHARGE=     .0400  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  ND1   CG     
      BOND  CD2   CG     
      BOND  CE1   ND1    
      BOND  NE2   CD2    
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  NE2   CE1    
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  NE2   HE2    
      BOND  CD2   HD2    
      BOND  CE1   HE1    
 
      DIHEDRAL  CA    CB    CG    CD2  ! multiple dihedral
      DIHEDRAL  CA    CB    CG    CD2  ! multiple dihedral
      DIHEDRAL  CA    CB    CG    CD2  ! multiple dihedral

      IMPROPER  NE2   CD2   CE1   HE2    
      IMPROPER  CD2   CG    NE2   HD2    
      IMPROPER  CE1   ND1   NE2   HE1    
      IMPROPER  NE2   CE1   CD2   HE2    
      IMPROPER  CD2   NE2   CG    HD2    
      IMPROPER  CE1   NE2   ND1   HE1    
 
      DONOR  HN    N      
      DONOR  HE2   NE2    
 
      ACCEPTOR  ND1   NONE   
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4532  106.43  123.52  111.67   1.5578  
      IC  N     C      *CA   HA       1.4532  106.43 -116.49  107.08   1.0833  
      IC  N     CA     CB    CG       1.4532  112.82  180.00  116.94   1.5109  
      IC  CG    CA     *CB   HB1      1.5109  116.94  119.80  107.91   1.1114  
      IC  CG    CA     *CB   HB2      1.5109  116.94 -124.04  109.50   1.1101  
      IC  CA    CB     CG    ND1      1.5578  116.94   90.00  120.17   1.3859  
      IC  ND1   CB     *CG   CD2      1.3859  120.17 -178.26  129.71   1.3596  
      IC  CB    CG     ND1   CE1      1.5109  120.17 -179.20  105.20   1.3170  
      IC  CB    CG     CD2   NE2      1.5109  129.71  178.66  105.80   1.3782  
      IC  NE2   ND1    *CE1  HE1      1.3539  111.76  179.69  124.58   1.0929  
      IC  CE1   CD2    *NE2  HE2      1.3539  107.15 -178.69  125.86    .9996  
      IC  NE2   CG     *CD2  HD2      1.3782  105.80 -179.35  129.89   1.0809  
 
 
 END {HSE }
 !-----------------------------------------------------------
 
 RESIDUE HSP  ! Protonated His 
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     |          HD1    HE1 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N           |     / 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1    ND1--CE1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |     /      | 
 GROUP						    !  HA-CA--CB--CG       | 
       ATOM  ND1   TYPE=NR3  CHARGE=    -.5500  END !     |   |     \      | 
       ATOM  HD1   TYPE=H    CHARGE=     .4500  END !     |   HB2    CD2--NE2 
       ATOM  NE2   TYPE=NR3  CHARGE=    -.5500  END !   O=C           |     \ 
       ATOM  HE2   TYPE=H    CHARGE=     .4500  END !     |          HD2    HE2 
       ATOM  CE1   TYPE=CPH2 CHARGE=     .1000  END  
       ATOM  HE1   TYPE=HR2  CHARGE=     .0000  END  
 GROUP  
       ATOM  CD2   TYPE=CPH1 CHARGE=     .0500  END  
       ATOM  HD2   TYPE=HR1  CHARGE=     .0000  END  
       ATOM  CG    TYPE=CPH1 CHARGE=     .0500  END  
       ATOM  CB    TYPE=CT2  CHARGE=    -.1000  END  
       ATOM  HB1   TYPE=HA   CHARGE=     .0500  END  
       ATOM  HB2   TYPE=HA   CHARGE=     .0500  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  ND1   CG     
      BOND  CD2   CG     
      BOND  CE1   ND1    
      BOND  NE2   CD2    
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  NE2   CE1    
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  ND1   HD1    
      BOND  NE2   HE2    
      BOND  CD2   HD2    
      BOND  CE1   HE1    
 
      DIHEDRAL  CA    CB    CG    CD2  ! multiple dihedral
      DIHEDRAL  CA    CB    CG    CD2  ! multiple dihedral
      DIHEDRAL  CA    CB    CG    CD2  ! multiple dihedral

      IMPROPER  HD1   CG    CE1   ND1    
      IMPROPER  HD1   CE1   CG    ND1    
      IMPROPER  HE2   CD2   CE1   NE2    
      IMPROPER  HE2   CE1   CD2   NE2    
 
      DONOR  HN    N      
      DONOR  HD1   ND1    
      DONOR  HE2   NE2    
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4548  112.03  125.13  109.38   1.5533  
      IC  N     C      *CA   HA       1.4548  112.03 -119.20  106.72   1.0832  
      IC  N     CA     CB    CG       1.4548  112.25  180.00  114.18   1.5168  
      IC  CG    CA     *CB   HB1      1.5168  114.18  122.50  108.99   1.1116  
      IC  CG    CA     *CB   HB2      1.5168  114.18 -121.51  108.97   1.1132  
      IC  CA    CB     CG    ND1      1.5533  114.18   90.00  122.94   1.3718  
      IC  ND1   CB     *CG   CD2      1.3718  122.94 -165.26  128.93   1.3549  
      IC  CB    CG     ND1   CE1      1.5168  122.94 -167.62  108.90   1.3262  
      IC  CB    CG     CD2   NE2      1.5168  128.93  167.13  106.93   1.3727  
      IC  NE2   ND1    *CE1  HE1      1.3256  108.50  178.39  125.76   1.0799  
      IC  CE1   CD2    *NE2  HE2      1.3256  108.82 -172.94  125.52   1.0020  
      IC  CE1   CG     *ND1  HD1      1.3262  108.90  171.49  126.09   1.0018  
      IC  NE2   CG     *CD2  HD2      1.3727  106.93 -174.49  128.41   1.0867  
 
 
 END {HSP }
 !-----------------------------------------------------------
 RESIDUE HIS  ! Protonated His - compatible with standard PDB (PDA 20-7-93)
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     |         HD1    HE1 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N          |    /  
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1   ND1--CE1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |     /     | 
 GROUP                                              !  HA-CA--CB--CG      |
       ATOM  ND1   TYPE=NR3  CHARGE=    -.5500  END !     |   |     \     | 
       ATOM  HD1   TYPE=H    CHARGE=     .4500  END !     |   HB2   CD2--NE2 
       ATOM  NE2   TYPE=NR3  CHARGE=    -.5500  END !   O=C          |	  \ 
       ATOM  HE2   TYPE=H    CHARGE=     .4500  END !     |         HD2    HE2 
       ATOM  CE1   TYPE=CPH2 CHARGE=     .1000  END  
       ATOM  HE1   TYPE=HR2  CHARGE=     .0000  END  
 GROUP  
       ATOM  CD2   TYPE=CPH1 CHARGE=     .0500  END  
       ATOM  HD2   TYPE=HR1  CHARGE=     .0000  END  
       ATOM  CG    TYPE=CPH1 CHARGE=     .0500  END  
       ATOM  CB    TYPE=CT2  CHARGE=    -.1000  END  
       ATOM  HB1   TYPE=HA   CHARGE=     .0500  END  
       ATOM  HB2   TYPE=HA   CHARGE=     .0500  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  ND1   CG     
      BOND  CD2   CG     
      BOND  CE1   ND1    
      BOND  NE2   CD2    
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  NE2   CE1    
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  ND1   HD1    
      BOND  NE2   HE2    
      BOND  CD2   HD2    
      BOND  CE1   HE1    
 
      DIHEDRAL  CA    CB    CG    CD2  ! multiple dihedral
      DIHEDRAL  CA    CB    CG    CD2  ! multiple dihedral
      DIHEDRAL  CA    CB    CG    CD2  ! multiple dihedral

      IMPROPER  HD1   CG    CE1   ND1    
      IMPROPER  HD1   CE1   CG    ND1    
      IMPROPER  HE2   CD2   CE1   NE2    
      IMPROPER  HE2   CE1   CD2   NE2    
 
      DONOR  HN    N      
      DONOR  HD1   ND1    
      DONOR  HE2   NE2    
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4548  112.03  125.13  109.38 1.5533  
      IC  N     C      *CA   HA       1.4548  112.03 -119.20  106.72 1.0832  
      IC  N     CA     CB    CG       1.4548  112.25  180.00  114.18 1.5168  
      IC  CG    CA     *CB   HB1      1.5168  114.18  122.50  108.99 1.1116  
      IC  CG    CA     *CB   HB2      1.5168  114.18 -121.51  108.97 1.1132  
      IC  CA    CB     CG    ND1      1.5533  114.18   90.00  122.94 1.3718  
      IC  ND1   CB     *CG   CD2      1.3718  122.94 -165.26  128.93 1.3549  
      IC  CB    CG     ND1   CE1      1.5168  122.94 -167.62  108.90 1.3262  
      IC  CB    CG     CD2   NE2      1.5168  128.93  167.13  106.93 1.3727  
      IC  NE2   ND1    *CE1  HE1      1.3256  108.50  178.39  125.76 1.0799  
      IC  CE1   CD2    *NE2  HE2      1.3256  108.82 -172.94  125.52 1.0020  
      IC  CE1   CG     *ND1  HD1      1.3262  108.90  171.49  126.09 1.0018  
      IC  NE2   CG     *CD2  HD2      1.3727  106.93 -174.49  128.41 1.0867  
 
 END {HIS }
 !-----------------------------------------------------------

 RESIDUE ILE   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     |    HG21 HG22 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N      | / 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |     CG2--HG23 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |    / 
 GROUP						    !  HA-CA--CB-HB    HD1 
       ATOM  CB    TYPE=CT1  CHARGE=    -.0900  END !     |    \       / 
       ATOM  HB    TYPE=HA   CHARGE=     .0900  END !     |     CG1--CD--HD2 
 GROUP						    !   O=C    / \     \ 
       ATOM  CG2   TYPE=CT3  CHARGE=    -.2700  END !     | HG11 HG12  HD3 
       ATOM  HG21  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG22  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG23  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CG1   TYPE=CT2  CHARGE=    -.1800  END  
       ATOM  HG11  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG12  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CD1    TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  HD11   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HD12   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HD13   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG1   CB     
      BOND  CG2   CB     
      BOND  CD1    CG1    
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB     
      BOND  CG1   HG11   
      BOND  CG1   HG12   
      BOND  CG2   HG21   
      BOND  CG2   HG22   
      BOND  CG2   HG23   
      BOND  CD1    HD11    
      BOND  CD1    HD12    
      BOND  CD1    HD13    
 
      DONOR  HN    N      
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4542  106.35  124.22  112.93   1.5681  
      IC  N     C      *CA   HA       1.4542  106.35 -115.63  106.81   1.0826  
      IC  N     CA     CB    CG1      1.4542  112.79  180.00  113.63   1.5498  
      IC  CG1   CA     *CB   HB       1.5498  113.63  114.55  104.48   1.1195  
      IC  CG1   CA     *CB   CG2      1.5498  113.63 -130.04  113.93   1.5452  
      IC  CA    CB     CG2   HG21     1.5681  113.93 -171.30  110.61   1.1100  
      IC  HG21  CB     *CG2  HG22     1.1100  110.61  119.35  110.90   1.1102  
      IC  HG21  CB     *CG2  HG23     1.1100  110.61 -120.09  110.97   1.1105  
      IC  CA    CB     CG1   CD1       1.5681  113.63  180.00  114.09   1.5381  
      IC  CD1    CB     *CG1  HG11     1.5381  114.09  122.36  109.78   1.1130  
      IC  CD1    CB     *CG1  HG12     1.5381  114.09 -120.59  108.89   1.1141  
      IC  CB    CG1    CD1    HD11      1.5498  114.09 -176.78  110.31   1.1115  
      IC  HD11   CG1    *CD1   HD12      1.1115  110.31  119.75  110.65   1.1113  
      IC  HD11   CG1    *CD1   HD13      1.1115  110.31 -119.70  111.02   1.1103  
 
 
 END {ILE }
 !-----------------------------------------------------------
 
 RESIDUE LEU   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     |        HD11 HD12 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N          | / 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1   CD1--HD13 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |    / 
 GROUP						    !  HA-CA--CB--CG-HG 
       ATOM  CB    TYPE=CT2  CHARGE=    -.1800  END !     |   |    \ 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     |   HB2   CD2--HD23 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   O=C          | \ 
 GROUP						    !     |        HD21 HD22 
       ATOM  CG    TYPE=CT1  CHARGE=    -.0900  END  
       ATOM  HG    TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CD1   TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  HD11  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HD12  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HD13  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CD2   TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  HD21  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HD22  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HD23  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  CD1   CG     
      BOND  CD2   CG     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CG    HG     
      BOND  CD1   HD11   
      BOND  CD1   HD12   
      BOND  CD1   HD13   
      BOND  CD2   HD21   
      BOND  CD2   HD22   
      BOND  CD2   HD23   
 
      DONOR  HN    N      
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4508  106.05  121.52  112.12   1.5543  
      IC  N     C      *CA   HA       1.4508  106.05 -116.50  107.57   1.0824  
      IC  N     CA     CB    CG       1.4508  111.19  180.00  117.46   1.5472  
      IC  CG    CA     *CB   HB1      1.5472  117.46  120.98  107.17   1.1145  
      IC  CG    CA     *CB   HB2      1.5472  117.46 -124.67  108.98   1.1126  
      IC  CA    CB     CG    CD1      1.5543  117.46  180.00  110.48   1.5361  
      IC  CD1   CB     *CG   CD2      1.5361  110.48 -123.75  112.57   1.5360  
      IC  CD1   CB     *CG   HG       1.5361  110.48  116.63  108.68   1.1168  
      IC  CB    CG     CD1   HD11     1.5472  110.48  177.33  110.54   1.1111  
      IC  HD11  CG     *CD1  HD12     1.1111  110.54  119.96  110.62   1.1112  
      IC  HD11  CG     *CD1  HD13     1.1111  110.54 -119.85  110.69   1.1108  
      IC  CB    CG     CD2   HD21     1.5472  112.57  178.96  110.32   1.1116  
      IC  HD21  CG     *CD2  HD22     1.1116  110.32  119.71  111.69   1.1086  
      IC  HD21  CG     *CD2  HD23     1.1116  110.32 -119.61  110.49   1.1115  
 
 
 END {LEU }
 !-----------------------------------------------------------
 
 RESIDUE LYS   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     | 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1 HG1 HD1 HE1    HZ1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |   |   |   |     / 
 GROUP						    !  HA-CA--CB--CG--CD--CE--NZ--HZ2 
       ATOM  CB    TYPE=CT2  CHARGE=    -.1800  END !     |   |   |   |   |     \ 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     |   HB2 HG2 HD2 HE2    HZ3 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   O=C 
 GROUP						    !     | 
       ATOM  CG    TYPE=CT2  CHARGE=    -.1800  END  
       ATOM  HG1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG2   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CD    TYPE=CT2  CHARGE=    -.1800  END  
       ATOM  HD1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HD2   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CE    TYPE=CT2  CHARGE=     .0000  END  
       ATOM  HE1   TYPE=HA   CHARGE=     .0000  END  
       ATOM  HE2   TYPE=HA   CHARGE=     .0000  END  
       ATOM  NZ    TYPE=NH3  CHARGE=    -.9000  END  
       ATOM  HZ1   TYPE=HC   CHARGE=     .3000  END  
       ATOM  HZ2   TYPE=HC   CHARGE=     .3000  END  
       ATOM  HZ3   TYPE=HC   CHARGE=     .3000  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  CD    CG     
      BOND  CE    CD     
      BOND  NZ    CE     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CG    HG1    
      BOND  CG    HG2    
      BOND  CD    HD1    
      BOND  CD    HD2    
      BOND  CE    HE1    
      BOND  CE    HE2    
      BOND  NZ    HZ1    
      BOND  NZ    HZ2    
      BOND  NZ    HZ3    
 
      DONOR  HN    N      
      DONOR  HZ1   NZ     
      DONOR  HZ2   NZ     
      DONOR  HZ3   NZ     
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4504  107.29  122.23  111.36   1.5568  
      IC  N     C      *CA   HA       1.4504  107.29 -116.88  107.36   1.0833  
      IC  N     CA     CB    CG       1.4504  111.47  180.00  115.76   1.5435  
      IC  CG    CA     *CB   HB1      1.5435  115.76  120.90  107.11   1.1146  
      IC  CG    CA     *CB   HB2      1.5435  115.76 -124.48  108.99   1.1131  
      IC  CA    CB     CG    CD       1.5568  115.76  180.00  113.28   1.5397  
      IC  CD    CB     *CG   HG1      1.5397  113.28  120.74  109.10   1.1138  
      IC  CD    CB     *CG   HG2      1.5397  113.28 -122.34  108.99   1.1143  
      IC  CB    CG     CD    CE       1.5435  113.28  180.00  112.33   1.5350  
      IC  CE    CG     *CD   HD1      1.5350  112.33  122.25  108.41   1.1141  
      IC  CE    CG     *CD   HD2      1.5350  112.33 -121.59  108.13   1.1146  
      IC  CG    CD     CE    NZ       1.5397  112.33  180.00  110.46   1.4604  
      IC  NZ    CD     *CE   HE1      1.4604  110.46  119.91  110.51   1.1128  
      IC  NZ    CD     *CE   HE2      1.4604  110.46 -120.02  110.57   1.1123  
      IC  CD    CE     NZ    HZ1      1.5350  110.46  179.92  110.02   1.0404  
      IC  HZ1   CE     *NZ   HZ2      1.0404  110.02  120.27  109.50   1.0402  
      IC  HZ1   CE     *NZ   HZ3      1.0404  110.02 -120.13  109.40   1.0401  
 
 
 END {LYS }
 !-----------------------------------------------------------
 
 RESIDUE MET   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     | 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1 HG1     HE1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |   |       | 
 GROUP						    !  HA-CA--CB--CG--SD--CE--HE3 
       ATOM  CB    TYPE=CT2  CHARGE=    -.1800  END !     |   |   |       | 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     |   HB2 HG2     HE2 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   O=C 
 GROUP						    !     | 
       ATOM  CG    TYPE=CT2  CHARGE=    -.1400  END  
       ATOM  HG1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG2   TYPE=HA   CHARGE=     .0900  END  
       ATOM  SD    TYPE=S    CHARGE=    -.0900  END  
       ATOM  CE    TYPE=CT3  CHARGE=    -.2200  END  
       ATOM  HE1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HE2   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HE3   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  SD    CG     
      BOND  CE    SD     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CG    HG1    
      BOND  CG    HG2    
      BOND  CE    HE1    
      BOND  CE    HE2    
      BOND  CE    HE3    

      DIHEDRAL CB  CG  SD  CE  ! multiple dihedral 
      DIHEDRAL CB  CG  SD  CE  ! multiple dihedral 

      DONOR  HN    N      
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4510  106.31  121.62  111.88   1.5546  
      IC  N     C      *CA   HA       1.4510  106.31 -116.98  107.57   1.0832  
      IC  N     CA     CB    CG       1.4510  111.25  180.00  115.92   1.5460  
      IC  CG    CA     *CB   HB1      1.5460  115.92  120.56  106.90   1.1153  
      IC  CG    CA     *CB   HB2      1.5460  115.92 -124.80  109.38   1.1129  
      IC  CA    CB     CG    SD       1.5546  115.92  180.00  110.28   1.8219  
      IC  SD    CB     *CG   HG1      1.8219  110.28  120.50  110.34   1.1106  
      IC  SD    CB     *CG   HG2      1.8219  110.28 -121.16  109.64   1.1119  
      IC  CB    CG     SD    CE       1.5460  110.28  180.00   98.94   1.8206  
      IC  CG    SD     CE    HE1      1.8219   98.94 -179.42  110.91   1.1111  
      IC  HE1   SD     *CE   HE2      1.1111  110.91  119.95  111.03   1.1115  
      IC  HE1   SD     *CE   HE3      1.1111  110.91 -119.95  111.09   1.1112  
 
 
 END {MET }
 !-----------------------------------------------------------
 
 RESIDUE PHE   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     |        HD1  HE1 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N         |    | 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1  CD1--CE1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |    /      \ 
 GROUP						    !  HA-CA--CB--CG      CZ--HZ 
       ATOM  CB    TYPE=CT2  CHARGE=    -.1800  END !     |   |    \      / 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     |   HB2  CD2--CE2 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   O=C         |    | 
 GROUP						    !     |        HD2  HE2 
       ATOM  CG    TYPE=CA   CHARGE=     .0000  END  
 GROUP  
       ATOM  CD1   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HD1   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CD2   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HD2   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CE1   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HE1   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CE2   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HE2   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CZ    TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HZ    TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  CD1   CG     
      BOND  CD2   CG     
      BOND  CE1   CD1    
      BOND  CE2   CD2    
      BOND  CZ    CE1    
      BOND  CZ    CE2    
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CD1   HD1    
      BOND  CD2   HD2    
      BOND  CE1   HE1    
      BOND  CE2   HE2    
      BOND  CZ    HZ     
 
      DONOR  HN    N      
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4504  106.38  122.49  112.45   1.5594  
      IC  N     C      *CA   HA       1.4504  106.38 -115.63  107.05   1.0832  
      IC  N     CA     CB    CG       1.4504  111.63  180.00  112.76   1.5109  
      IC  CG    CA     *CB   HB1      1.5109  112.76  118.27  109.10   1.1119  
      IC  CG    CA     *CB   HB2      1.5109  112.76 -123.83  111.11   1.1113  
      IC  CA    CB     CG    CD1      1.5594  112.76   90.00  120.32   1.4059  
      IC  CD1   CB     *CG   CD2      1.4059  120.32 -177.96  120.76   1.4062  
      IC  CB    CG     CD1   CE1      1.5109  120.32 -177.37  120.63   1.4006  
      IC  CE1   CG     *CD1  HD1      1.4006  120.63  179.70  119.65   1.0814  
      IC  CB    CG     CD2   CE2      1.5109  120.76  177.20  120.62   1.4002  
      IC  CE2   CG     *CD2  HD2      1.4002  120.62 -178.69  119.99   1.0811  
      IC  CG    CD1    CE1   CZ       1.4059  120.63    -.12  119.93   1.4004  
      IC  CZ    CD1    *CE1  HE1      1.4004  119.93 -179.69  120.01   1.0808  
      IC  CZ    CD2    *CE2  HE2      1.4000  119.96 -179.93  119.87   1.0811  
      IC  CE1   CE2    *CZ   HZ       1.4004  119.98  179.51  119.97   1.0807  
 
 
 END {PHE }
 !-----------------------------------------------------------
 
 RESIDUE PRO   
 
 GROUP						    !       HD1 HD2 
       ATOM  N     TYPE=N    CHARGE=    -.2900  END !     |   \ / 
       ATOM  CA    TYPE=CP1  CHARGE=     .0200  END !     N---CD   HG1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |     \  / 
       ATOM  CD    TYPE=CP3  CHARGE=     .0000  END !     |      CG 
       ATOM  HD1   TYPE=HA   CHARGE=     .0900  END !     |     /  \ 
       ATOM  HD2   TYPE=HA   CHARGE=     .0900  END !  HA-CA--CB   HG2 
 GROUP						    !     |   / \ 
       ATOM  CB    TYPE=CP2  CHARGE=    -.1800  END !     | HB1 HB2 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !   O=C 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !     | 
 GROUP  
       ATOM  CG    TYPE=CP2  CHARGE=    -.1800  END  
       ATOM  HG1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG2   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  O     C      
      BOND  C     CA     
      BOND  N     CA     
      BOND  CA    CB     
      BOND  CB    CG     
      BOND  CG    CD     
      BOND  N     CD     
      BOND  HA    CA     
      BOND  HG1   CG     
      BOND  HG2   CG     
      BOND  HD1   CD     
      BOND  HD2   CD     
      BOND  HB1   CB     
      BOND  HB2   CB     

      DIHEDRAL CB  CA  C   O   ! multiple dihedral
      DIHEDRAL CB  CA  C   O   ! multiple dihedral
      DIHEDRAL HA  CA  C   O   ! multiple dihedral
      DIHEDRAL HA  CA  C   O   ! multiple dihedral

      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4585  110.86  113.74  111.74   1.5399  
      IC  N     C      *CA   HA       1.4585  110.86 -122.40  109.09   1.0837  
      IC  N     CA     CB    CG       1.4585  102.56   31.61  104.39   1.5322  
      IC  CA    CB     CG    CD       1.5399  104.39  -34.59  103.21   1.5317  
      IC  N     CA     CB    HB1      1.4585  102.56  -84.94  109.02   1.1131  
      IC  N     CA     CB    HB2      1.4585  102.56  153.93  112.74   1.1088  
      IC  CA    CB     CG    HG1      1.5399  104.39 -156.72  112.95   1.1077  
      IC  CA    CB     CG    HG2      1.5399  104.39   81.26  109.22   1.1143  
      IC  CB    CG     CD    HD1      1.5322  103.21  -93.55  110.03   1.1137  
      IC  CB    CG     CD    HD2      1.5322  103.21  144.52  110.00   1.1144  

 END {PRO }
 !-----------------------------------------------------------
 
 RESIDUE SER   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     | 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   | 
 GROUP						    !  HA-CA--CB--OG 
       ATOM  CB    TYPE=CT2  CHARGE=     .0500  END !     |   |     \ 
       ATOM  HB1   TYPE=HA   CHARGE=     .1000  END !     |   HB2    HG1 
       ATOM  HB2   TYPE=HA   CHARGE=     .1000  END !   O=C 
       ATOM  OG    TYPE=OH1  CHARGE=    -.5500  END !     | 
       ATOM  HG    TYPE=H    CHARGE=     .3000  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  OG    CB     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  OG    HG    

      DIHEDRAL CA  CB  OG  HG ! multiple dihedral
      DIHEDRAL CA  CB  OG  HG ! multiple dihedral
      DIHEDRAL CA  CB  OG  HG ! multiple dihedral
 
      DONOR  HN    N      
      DONOR  HG1   OG     
 
      ACCEPTOR  OG    NONE   
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4579  105.81  124.75  111.40   1.5585  
      IC  N     C      *CA   HA       1.4579  105.81 -115.56  107.30   1.0821  
      IC  N     CA     CB    OG       1.4579  114.28  180.00  112.45   1.4341  
      IC  OG    CA     *CB   HB1      1.4341  112.45  119.32  108.10   1.1140  
      IC  OG    CA     *CB   HB2      1.4341  112.45 -123.86  110.38   1.1136  
      IC  CA    CB     OG    HG1      1.5585  112.45  165.96  107.08    .9655  
 
 
 END {SER }
 !-----------------------------------------------------------
 
 RESIDUE THR   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     | 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |     OG1--HG1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |    / 
 GROUP						    !  HA-CA--CB-HB 
       ATOM  CB    TYPE=CT1  CHARGE=     .1600  END !     |    \ 
       ATOM  HB    TYPE=HA   CHARGE=     .0900  END !     |     CG2--HG21 
       ATOM  OG1   TYPE=OH1  CHARGE=    -.5500  END !   O=C    / \ 
       ATOM  HG1   TYPE=H    CHARGE=     .3000  END !     | HG21 HG22 
 GROUP  
       ATOM  CG2   TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  HG21  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG22  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG23  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  OG1   CB     
      BOND  CG2   CB     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB     
      BOND  OG1   HG1    
      BOND  CG2   HG21   
      BOND  CG2   HG22   
      BOND  CG2   HG23   

      DIHEDRAL CA  CB  OG1 HG1 ! multiple dihedral
      DIHEDRAL CA  CB  OG1 HG1 ! multiple dihedral
      DIHEDRAL CA  CB  OG1 HG1 ! multiple dihedral
      DIHEDRAL CG2 CB  OG1 HG1 ! multiple dihedral
      DIHEDRAL CG2 CB  OG1 HG1 ! multiple dihedral
      DIHEDRAL CG2 CB  OG1 HG1 ! multiple dihedral

      DONOR  HN    N      
      DONOR  HG1   OG1    
 
      ACCEPTOR  OG1   NONE   
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4607  106.09  126.46  112.74   1.5693  
      IC  N     C      *CA   HA       1.4607  106.09 -114.92  106.53   1.0817  
      IC  N     CA     CB    OG1      1.4607  114.81  180.00  112.16   1.4252  
      IC  OG1   CA     *CB   HB       1.4252  112.16  116.39  106.11   1.1174  
      IC  OG1   CA     *CB   CG2      1.4252  112.16 -124.13  115.91   1.5324  
      IC  CA    CB     OG1   HG1      1.5693  112.16 -179.28  105.45    .9633  
      IC  CA    CB     CG2   HG21     1.5693  115.91 -173.65  110.85   1.1104  
      IC  HG21  CB     *CG2  HG22     1.1104  110.85  119.51  110.41   1.1109  
      IC  HG21  CB     *CG2  HG23     1.1104  110.85 -120.39  111.11   1.1113  
 
 
 END {THR }
 !-----------------------------------------------------------
 
 RESIDUE TRP   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     |                  HE3 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N                   | 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1            CE3 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |             /   \ 
 GROUP						    !  HA-CA--CB---CG-----CD2   CZ3-HZ3 
       ATOM  CB    TYPE=CT2  CHARGE=    -.1800  END !     |   |     |      |     | 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     |   HB2  CD1    CE2   CH2-HH2 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   O=C       /   \   / \   / 
 GROUP						    !     |     HD1    NE1   CZ2 
       ATOM  CG    TYPE=CY   CHARGE=    -.0300  END !                   |     | 
       ATOM  CD2   TYPE=CPT  CHARGE=    -.0200  END !                  HE1   HZ2 
       ATOM  CD1   TYPE=CA   CHARGE=     .0350  END  
       ATOM  HD1   TYPE=HP   CHARGE=     .1150  END  
       ATOM  NE1   TYPE=NY   CHARGE=    -.6100  END  
       ATOM  HE1   TYPE=H    CHARGE=     .3800  END  
       ATOM  CE2   TYPE=CPT  CHARGE=     .1300  END  
 GROUP  
       ATOM  CE3   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HE3   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CZ2   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HZ2   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CZ3   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HZ3   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CH2   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HH2   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  CD1   CG     
      BOND  CD2   CG     
      BOND  NE1   CD1    
      BOND  CE2   CD2    
      BOND  CZ2   CE2    
      BOND  CZ3   CE3    
      BOND  CH2   CZ2    
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CZ3   CH2    
      BOND  CD2   CE3    
      BOND  NE1   CE2    
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CD1   HD1    
      BOND  NE1   HE1    
      BOND  CE3   HE3    
      BOND  CZ2   HZ2    
      BOND  CZ3   HZ3    
      BOND  CH2   HH2    
 
      DONOR  HN    N      
      DONOR  HE1   NE1    
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4507  107.69  122.68  111.23   1.5560  
      IC  N     C      *CA   HA       1.4507  107.69 -117.02  106.92   1.0835  
      IC  N     CA     CB    CG       1.4507  111.68  180.00  115.14   1.5233  
      IC  CG    CA     *CB   HB1      1.5233  115.14  119.17  107.84   1.1127  
      IC  CG    CA     *CB   HB2      1.5233  115.14 -124.73  109.87   1.1118  
      IC  CA    CB     CG    CD2      1.5560  115.14   90.00  123.95   1.4407  
      IC  CD2   CB     *CG   CD1      1.4407  123.95 -172.81  129.18   1.3679  
      IC  CD1   CG     CD2   CE2      1.3679  106.57    -.08  106.65   1.4126  
      IC  CG    CD2    CE2   NE1      1.4407  106.65     .14  107.87   1.3746  
      IC  CE2   CG     *CD2  CE3      1.4126  106.65  179.21  132.54   1.4011  
      IC  CE2   CD2    CE3   CZ3      1.4126  120.80    -.20  118.16   1.4017  
      IC  CD2   CE3    CZ3   CH2      1.4011  118.16     .10  120.97   1.4019  
      IC  CE3   CZ3    CH2   CZ2      1.4017  120.97     .01  120.87   1.4030  
      IC  CZ3   CD2    *CE3  HE3      1.4017  118.16 -179.62  121.84   1.0815  
      IC  CH2   CE3    *CZ3  HZ3      1.4019  120.97 -179.82  119.45   1.0811  
      IC  CZ2   CZ3    *CH2  HH2      1.4030  120.87 -179.92  119.57   1.0811  
      IC  CE2   CH2    *CZ2  HZ2      1.3939  118.42  179.87  120.08   1.0790  
      IC  CD1   CE2    *NE1  HE1      1.3752  108.81  177.78  124.68    .9767  
      IC  CG    NE1    *CD1  HD1      1.3679  110.10  178.10  125.43   1.0820  
 
 
 END {TRP }
 !-----------------------------------------------------------
 
 RESIDUE TYR   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     |        HD1  HE1 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N         |    | 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |   HB1  CD1--CE1 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |   |    /      \ 
 GROUP						    !  HA-CA--CB--CG      CZ--OH 
       ATOM  CB    TYPE=CT2  CHARGE=    -.1800  END !     |   |    \      /     \ 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     |   HB2  CD2--CE2     HH 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   O=C         |    | 
 GROUP						    !     |        HD2  HE2 
       ATOM  CG    TYPE=CA   CHARGE=     .0000  END  
 GROUP  
       ATOM  CD1   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HD1   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CD2   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HD2   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CE1   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HE1   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CE2   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HE2   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CZ    TYPE=CA   CHARGE=     .2500  END  
       ATOM  OH    TYPE=OH1  CHARGE=    -.5500  END  
       ATOM  HH    TYPE=H    CHARGE=     .3000  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG    CB     
      BOND  CD1   CG     
      BOND  CD2   CG     
      BOND  CE1   CD1    
      BOND  CE2   CD2    
      BOND  CZ    CE1    
      BOND  CZ    CE2    
      BOND  OH    CZ     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CD1   HD1    
      BOND  CD2   HD2    
      BOND  CE1   HE1    
      BOND  CE2   HE2    
      BOND  OH    HH     
 
      DONOR  HN    N      
      DONOR  HH    OH     
 
      ACCEPTOR  OH    NONE   
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4501  106.52  122.27  112.34   1.5606  
      IC  N     C      *CA   HA       1.4501  106.52 -116.04  107.15   1.0833  
      IC  N     CA     CB    CG       1.4501  111.43  180.00  112.94   1.5113  
      IC  CG    CA     *CB   HB1      1.5113  112.94  118.89  109.12   1.1119  
      IC  CG    CA     *CB   HB2      1.5113  112.94 -123.36  110.70   1.1115  
      IC  CA    CB     CG    CD1      1.5606  112.94   90.00  120.49   1.4064  
      IC  CD1   CB     *CG   CD2      1.4064  120.49 -176.46  120.46   1.4068  
      IC  CB    CG     CD1   CE1      1.5113  120.49 -175.49  120.40   1.4026  
      IC  CE1   CG     *CD1  HD1      1.4026  120.40  178.94  119.80   1.0814  
      IC  CB    CG     CD2   CE2      1.5113  120.46  175.32  120.56   1.4022  
      IC  CE2   CG     *CD2  HD2      1.4022  120.56 -177.57  119.98   1.0813  
      IC  CG    CD1    CE1   CZ       1.4064  120.40    -.19  120.09   1.3978  
      IC  CZ    CD1    *CE1  HE1      1.3978  120.09  179.64  120.58   1.0799  
      IC  CZ    CD2    *CE2  HE2      1.3979  119.92 -178.69  119.76   1.0798  
      IC  CE1   CE2    *CZ   OH       1.3978  120.05 -178.98  120.25   1.4063  
      IC  CE1   CZ     OH    HH       1.3978  119.68  175.45  107.47    .9594  
 
 
 END {TYR }
 !-----------------------------------------------------------
 
 RESIDUE VAL   
 
 GROUP  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END !     |    HG11 HG12 
       ATOM  HN    TYPE=H    CHARGE=     .3100  END !  HN-N      | / 
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END !     |     CG1--HG13 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !     |    / 
 GROUP						    !  HA-CA--CB-HB 
       ATOM  CB    TYPE=CT1  CHARGE=    -.0900  END !     |    \ 
       ATOM  HB    TYPE=HA   CHARGE=     .0900  END !     |     CG2--HG21 
 GROUP						    !   O=C    / \ 
       ATOM  CG1   TYPE=CT3  CHARGE=    -.2700  END !     | HG21 HG22 
       ATOM  HG11  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG12  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG13  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CG2   TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  HG21  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG22  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HG23  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
 !END GROUP
 
      BOND  CB    CA     
      BOND  CG1   CB     
      BOND  CG2   CB     
      BOND  N     HN     
      BOND  N     CA     
      BOND  O     C      
      BOND  C     CA     
      BOND  CA    HA     
      BOND  CB    HB     
      BOND  CG1   HG11   
      BOND  CG1   HG12   
      BOND  CG1   HG13   
      BOND  CG2   HG21   
      BOND  CG2   HG22   
      BOND  CG2   HG23   
 
      DONOR  HN    N      
 
      ACCEPTOR  O     C      
 
      IC  N     C      *CA   CB       1.4570  105.54  122.95  111.23   1.5660  
      IC  N     C      *CA   HA       1.4570  105.54 -117.24  107.46   1.0828  
      IC  N     CA     CB    CG1      1.4570  113.05  180.00  113.97   1.5441  
      IC  CG1   CA     *CB   CG2      1.5441  113.97  123.99  112.17   1.5414  
      IC  CG1   CA     *CB   HB       1.5441  113.97 -119.17  107.57   1.1178  
      IC  CA    CB     CG1   HG11     1.5660  113.97  177.83  110.30   1.1114  
      IC  HG11  CB     *CG1  HG12     1.1114  110.30  119.25  111.67   1.1097  
      IC  HG11  CB     *CG1  HG13     1.1114  110.30 -119.49  110.70   1.1110  
      IC  CA    CB     CG2   HG21     1.5660  112.17 -177.78  110.71   1.1108  
      IC  HG21  CB     *CG2  HG22     1.1108  110.71  120.08  110.56   1.1115  
      IC  HG21  CB     *CG2  HG23     1.1108  110.71 -119.55  111.23   1.1098  
 
 
 END {VAL }
 !-----------------------------------------------------------
 
 RESIDUE TIP3 ! tip3p water model, generate using noangle nodihedral 
 
 GROUP  
       ATOM  OH2   TYPE=OT   CHARGE=    -.8340  END  
       ATOM  H1    TYPE=HT   CHARGE=     .4170  END  
       ATOM  H2    TYPE=HT   CHARGE=     .4170  END  
 !END GROUP
 
      BOND  OH2   H1     
      BOND  OH2   H2     
!      BOND  H1    H2    ! the last bond is needed for shake (not in Xplor!!!)
 
      ANGLE  H1    OH2   H2    ! required 
 
      ACCEPTOR  OH2   NONE   
 
 
 END {TIP3}
 !-----------------------------------------------------------
 
 RESIDUE CAL  ! Calcium ion 
 
 GROUP  
       ATOM  CAL   TYPE=CAL  CHARGE=    2.0000  END  
 !END GROUP
 
 
 END {CAL }
 !-----------------------------------------------------------
 
 RESIDUE ZN2  ! Zinc ion 
 
 GROUP  
       ATOM  ZN    TYPE=ZN   CHARGE=    2.0000  END  
 !END GROUP
 
 
 END {ZN2 }
 !-----------------------------------------------------------
 
 RESIDUE HEME ! 6-liganded planar heme 
 
 GROUP  
       ATOM  FE    TYPE=FE   CHARGE=     .2400  END !       O2A   O1A             O2D  O1D 
       ATOM  NA    TYPE=NPH  CHARGE=    -.1800  END !         \\ //                 \\ // 
       ATOM  NB    TYPE=NPH  CHARGE=    -.1800  END !          CGA                   CGD 
       ATOM  NC    TYPE=NPH  CHARGE=    -.1800  END !           |                     | 
       ATOM  ND    TYPE=NPH  CHARGE=    -.1800  END !    HBA1--CBA--HBA2  HA   HBD1--CBD--HBD2 
       ATOM  C1A   TYPE=CPA  CHARGE=     .1200  END !           |          |          | 
       ATOM  C2A   TYPE=CPB  CHARGE=    -.0600  END !    HAA1--CAA-HAA2  _CHA_ HAD1--CAD--HAD2 
       ATOM  C3A   TYPE=CPB  CHARGE=    -.0600  END !           |       /     \       | 
       ATOM  C4A   TYPE=CPA  CHARGE=     .1200  END !          C2A---C1A       C4D---C3D 
       ATOM  C1B   TYPE=CPA  CHARGE=     .1200  END !           |     |         |     | 
       ATOM  C2B   TYPE=CPB  CHARGE=    -.0600  END !HMA1\      |     |         |     |      /HMD1 
       ATOM  C3B   TYPE=CPB  CHARGE=    -.0600  END !HMA2-CMA--C3A    NA       ND    C2D--CMD-HMD2 
       ATOM  C4B   TYPE=CPA  CHARGE=     .1200  END !HMA3/       \   / \       / \   /       \HMD3 
       ATOM  C1C   TYPE=CPA  CHARGE=     .1200  END !             C4A   \     /   C1D 
       ATOM  C2C   TYPE=CPB  CHARGE=    -.0600  END !            /       \   /       \ 
       ATOM  C3C   TYPE=CPB  CHARGE=    -.0600  END !      HB--CHB        FE         CHD--HD 
       ATOM  C4C   TYPE=CPA  CHARGE=     .1200  END !            \       /   \       / 
       ATOM  C1D   TYPE=CPA  CHARGE=     .1200  END !             C1B   /     \   C4C        HAC 
       ATOM  C2D   TYPE=CPB  CHARGE=    -.0600  END !HMB1\       /   \ /       \ /   \       / 
       ATOM  C3D   TYPE=CPB  CHARGE=    -.0600  END !HMB2-CMB--C2B    NB       NC    C3C--CAC 
       ATOM  C4D   TYPE=CPA  CHARGE=     .1200  END !HMB3/      |     |         |     |     \\  /HBC1 
 GROUP						    !           |     |         |     |      CBC 
       ATOM  CHA   TYPE=CPM  CHARGE=    -.1000  END !          C3B---C4B       C1C---C2C        \HBC2 
       ATOM  HA    TYPE=HA   CHARGE=     .1000  END !           |       \_CHC_/       | 
 GROUP						    !          CAB         |         CMC--HMC3 
       ATOM  CHB   TYPE=CPM  CHARGE=    -.1000  END !         //  \        HC       /  | 
       ATOM  HB    TYPE=HA   CHARGE=     .1000  END !        CBB  HAB           HMC1  HMC2 
 GROUP						    !       /   \ 
       ATOM  CHC   TYPE=CPM  CHARGE=    -.1000  END !    HBB1  HBB2 
       ATOM  HC    TYPE=HA   CHARGE=     .1000  END ! 
 GROUP  
       ATOM  CHD   TYPE=CPM  CHARGE=    -.1000  END  
       ATOM  HD    TYPE=HA   CHARGE=     .1000  END  
 GROUP  
       ATOM  CMA   TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  HMA1  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HMA2  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HMA3  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CAA   TYPE=CT2  CHARGE=    -.1800  END  
       ATOM  HAA1  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HAA2  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CBA   TYPE=CT2  CHARGE=    -.2800  END  
       ATOM  HBA1  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HBA2  TYPE=HA   CHARGE=     .0900  END  
       ATOM  CGA   TYPE=CC   CHARGE=     .6200  END  
       ATOM  O1A   TYPE=OC   CHARGE=    -.7600  END  
       ATOM  O2A   TYPE=OC   CHARGE=    -.7600  END  
 GROUP  
       ATOM  CMB   TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  HMB1  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HMB2  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HMB3  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CAB   TYPE=C    CHARGE=    -.2000  END  
       ATOM  HAB   TYPE=HA   CHARGE=     .2000  END  
 GROUP  
       ATOM  CBB   TYPE=C    CHARGE=    -.2000  END  
       ATOM  HBB1  TYPE=HA   CHARGE=     .1000  END  
       ATOM  HBB2  TYPE=HA   CHARGE=     .1000  END  
 GROUP  
       ATOM  CMC   TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  HMC1  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HMC2  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HMC3  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CAC   TYPE=C    CHARGE=    -.2000  END  
       ATOM  HAC   TYPE=HA   CHARGE=     .2000  END  
 GROUP  
       ATOM  CBC   TYPE=C    CHARGE=    -.2000  END  
       ATOM  HBC1  TYPE=HA   CHARGE=     .1000  END  
       ATOM  HBC2  TYPE=HA   CHARGE=     .1000  END  
 GROUP  
       ATOM  CMD   TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  HMD1  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HMD2  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HMD3  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CAD   TYPE=CT2  CHARGE=    -.1800  END  
       ATOM  HAD1  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HAD2  TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CBD   TYPE=CT2  CHARGE=    -.2800  END  
       ATOM  HBD1  TYPE=HA   CHARGE=     .0900  END  
       ATOM  HBD2  TYPE=HA   CHARGE=     .0900  END  
       ATOM  CGD   TYPE=CC   CHARGE=     .6200  END  
       ATOM  O1D   TYPE=OC   CHARGE=    -.7600  END  
       ATOM  O2D   TYPE=OC   CHARGE=    -.7600  END  
 !END GROUP
 
      BOND  FE    NA     
      BOND  FE    NB     
      BOND  FE    NC     
      BOND  FE    ND     
      BOND  NA    C1A    
      BOND  C1A   C2A    
      BOND  C2A   C3A    
      BOND  C3A   C4A    
      BOND  NA    C4A    
      BOND  C2A   CAA    
      BOND  CAA   CBA    
      BOND  CBA   CGA    
      BOND  CGA   O1A    
      BOND  CGA   O2A    
      BOND  C3A   CMA    
      BOND  CHB   C4A    
      BOND  CHB   C1B    
      BOND  NB    C1B    
      BOND  C1B   C2B    
      BOND  C2B   C3B    
      BOND  C3B   C4B    
      BOND  NB    C4B    
      BOND  C2B   CMB    
      BOND  C3B   CAB    
      BOND  CAB   CBB    
      BOND  CHC   C4B    
      BOND  CHC   C1C    
      BOND  NC    C1C    
      BOND  C1C   C2C    
      BOND  C2C   C3C    
      BOND  C3C   C4C    
      BOND  NC    C4C    
      BOND  C2C   CMC    
      BOND  C3C   CAC    
      BOND  CAC   CBC    
      BOND  CHD   C4C    
      BOND  CHD   C1D    
      BOND  ND    C1D    
      BOND  C1D   C2D    
      BOND  C2D   C3D    
      BOND  C3D   C4D    
      BOND  ND    C4D    
      BOND  C2D   CMD    
      BOND  C3D   CAD    
      BOND  CAD   CBD    
      BOND  CBD   CGD    
      BOND  CGD   O1D    
      BOND  CGD   O2D    
      BOND  CHA   C4D    
      BOND  CHA   C1A    
      BOND  CHA   HA     
      BOND  CHB   HB     
      BOND  CHC   HC     
      BOND  CHD   HD     
      BOND  CAA   HAA1   
      BOND  CAA   HAA2   
      BOND  CBA   HBA1   
      BOND  CBA   HBA2   
      BOND  CMA   HMA1   
      BOND  CMA   HMA2   
      BOND  CMA   HMA3   
      BOND  CMB   HMB1   
      BOND  CMB   HMB2   
      BOND  CMB   HMB3   
      BOND  CAB   HAB    
      BOND  CBB   HBB1   
      BOND  CBB   HBB2   
      BOND  CMC   HMC1   
      BOND  CMC   HMC2   
      BOND  CMC   HMC3   
      BOND  CAC   HAC    
      BOND  CBC   HBC1   
      BOND  CBC   HBC2   
      BOND  CMD   HMD1   
      BOND  CMD   HMD2   
      BOND  CMD   HMD3   
      BOND  CAD   HAD1   
      BOND  CAD   HAD2   
      BOND  CBD   HBD1   
      BOND  CBD   HBD2   
 
      IMPROPER  C2A   C1A   C3A   CAA    
      IMPROPER  C3A   C2A   C4A   CMA    
      IMPROPER  C2B   C1B   C3B   CMB    
      IMPROPER  C3B   C2B   C4B   CAB    
      IMPROPER  C2C   C1C   C3C   CMC    
      IMPROPER  C3C   C2C   C4C   CAC    
      IMPROPER  C2D   C1D   C3D   CMD    
      IMPROPER  C3D   C2D   C4D   CAD    
      IMPROPER  O1A   CBA   O2A   CGA    
      IMPROPER  O1D   CBD   O2D   CGD    
      IMPROPER  C4A   NA    C1A   C2A    
      IMPROPER  C1A   NA    C4A   C3A    
      IMPROPER  C4B   NB    C1B   C2B    
      IMPROPER  C1B   NB    C4B   C3B    
      IMPROPER  C4C   NC    C1C   C2C    
      IMPROPER  C1C   NC    C4C   C3C    
      IMPROPER  C4D   ND    C1D   C2D    
      IMPROPER  C1D   ND    C4D   C3D    
      IMPROPER  NA    C1A   C2A   C3A    
      IMPROPER  NA    C4A   C3A   C2A    
      IMPROPER  NB    C1B   C2B   C3B    
      IMPROPER  NB    C4B   C3B   C2B    
      IMPROPER  NC    C1C   C2C   C3C    
      IMPROPER  NC    C4C   C3C   C2C    
      IMPROPER  ND    C1D   C2D   C3D    
      IMPROPER  ND    C4D   C3D   C2D    
      IMPROPER  NA    C1A   CHA   C4D    
      IMPROPER  NA    C4A   CHB   C1B    
      IMPROPER  NB    C1B   CHB   C4A    
      IMPROPER  NB    C4B   CHC   C1C    
      IMPROPER  NC    C1C   CHC   C4B    
      IMPROPER  NC    C4C   CHD   C1D    
      IMPROPER  ND    C1D   CHD   C4C    
      IMPROPER  ND    C4D   CHA   C1A    
      IMPROPER  CHA   C1A   C4D   HA     
      IMPROPER  CHB   C1B   C4A   HB     
      IMPROPER  CHC   C1C   C4B   HC     
      IMPROPER  CHD   C1D   C4C   HD     
      IMPROPER  C1A   C2A   CHA   NA     
      IMPROPER  C4A   C3A   CHB   NA     
      IMPROPER  C1B   C2B   CHB   NB     
      IMPROPER  C4B   C3B   CHC   NB     
      IMPROPER  C1C   C2C   CHC   NC     
      IMPROPER  C4C   C3C   CHD   NC     
      IMPROPER  C1D   C2D   CHD   ND     
      IMPROPER  C4D   C3D   CHA   ND     
      IMPROPER  NA    C1A   C4A   FE     
      IMPROPER  NB    C1B   C4B   FE     
      IMPROPER  NC    C1C   C4C   FE     
      IMPROPER  ND    C1D   C4D   FE     
      IMPROPER  CAB   CBB   C3B   HAB    
      IMPROPER  HAB   CAB   CBB   HBB2   
      IMPROPER  CAB   CBB   HBB2  HBB1   
      IMPROPER  CAC   CBC   C3C   HAC    
      IMPROPER  HAC   CAC   CBC   HBC2   
      IMPROPER  CAC   CBC   HBC2  HBC1   
 
      ACCEPTOR  NA    NONE   
      ACCEPTOR  O1A   NONE   
      ACCEPTOR  O2A   NONE   
      ACCEPTOR  NB    NONE   
      ACCEPTOR  NC    NONE   
      ACCEPTOR  ND    NONE   
      ACCEPTOR  O1D   NONE   
      ACCEPTOR  O2D   NONE   
 
      IC  FE    NA     C4A   C3A       .0000     .00  180.00     .00    .0000  
      IC  C3A   C4A    NA    C1A       .0000     .00     .00     .00    .0000  
      IC  FE    NA     C1A   C2A       .0000     .00  180.00     .00    .0000  
      IC  C4A   NA     FE    NB        .0000     .00     .00     .00    .0000  
      IC  NA    FE     NB    C1B       .0000     .00     .00     .00    .0000  
      IC  FE    NB     C1B   C2B       .0000     .00  180.00     .00    .0000  
      IC  C2B   C1B    NB    C4B       .0000     .00     .00     .00    .0000  
      IC  FE    NB     C4B   C3B       .0000     .00  180.00     .00    .0000  
      IC  C4B   NB     FE    NC        .0000     .00     .00     .00    .0000  
      IC  NB    FE     NC    C1C       .0000     .00     .00     .00    .0000  
      IC  FE    NC     C1C   C2C       .0000     .00  180.00     .00    .0000  
      IC  C2C   C1C    NC    C4C       .0000     .00     .00     .00    .0000  
      IC  FE    NC     C4C   C3C       .0000     .00  180.00     .00    .0000  
      IC  C4C   NC     FE    ND        .0000     .00     .00     .00    .0000  
      IC  NC    FE     ND    C1D       .0000     .00     .00     .00    .0000  
      IC  FE    ND     C1D   C2D       .0000     .00  180.00     .00    .0000  
      IC  C2D   C1D    ND    C4D       .0000     .00     .00     .00    .0000  
      IC  FE    ND     C4D   C3D       .0000     .00  180.00     .00    .0000  
      IC  C4A   NA     FE    CHB       .0000     .00     .00     .00    .0000  
      IC  NA    FE     CHB   HB        .0000     .00     .00     .00    .0000  
      IC  C4B   NB     FE    CHC       .0000     .00     .00     .00    .0000  
      IC  NB    FE     CHC   HC        .0000     .00     .00     .00    .0000  
      IC  C4C   NC     FE    CHD       .0000     .00     .00     .00    .0000  
      IC  NC    FE     CHD   HD        .0000     .00     .00     .00    .0000  
      IC  C4D   ND     FE    CHA       .0000     .00     .00     .00    .0000  
      IC  ND    FE     CHA   HA        .0000     .00     .00     .00    .0000  
      IC  C3B   C1B    *C2B  CMB       .0000     .00  180.00     .00    .0000  
      IC  C4B   C2B    *C3B  CAB       .0000     .00  180.00     .00    .0000  
      IC  C2B   C3B    CAB   CBB       .0000     .00  -45.00     .00    .0000 ! PREVENTS VINYL COLLISION 
      IC  CHC   C1C    C2C   CMC       .0000     .00     .00     .00    .0000  
      IC  C4C   C2C    *C3C  CAC       .0000     .00  180.00     .00    .0000  
      IC  C2C   C3C    CAC   CBC       .0000     .00  -45.00     .00    .0000 ! PREVENTS VINYL COLLISION 
      IC  C4D   C2D    *C3D  CAD       .0000     .00  180.00     .00    .0000  
      IC  C3D   C1D    *C2D  CMD       .0000     .00  180.00     .00    .0000  
      IC  C2D   C3D    CAD   CBD       .0000     .00 -120.00     .00    .0000  
      IC  C3D   CAD    CBD   CGD       .0000     .00  180.00     .00    .0000  
      IC  CAD   CBD    CGD   O1D       .0000     .00     .00     .00    .0000  
      IC  CAD   CBD    CGD   O2D       .0000     .00  180.00     .00    .0000  
      IC  C4A   C2A    *C3A  CMA       .0000     .00  180.00     .00    .0000  
      IC  C3A   C1A    *C2A  CAA       .0000     .00  180.00     .00    .0000  
      IC  C1A   C2A    CAA   CBA       .0000     .00  120.00     .00    .0000  
      IC  C2A   CAA    CBA   CGA       .0000     .00  180.00     .00    .0000  
      IC  CAA   CBA    CGA   O1A       .0000     .00     .00     .00    .0000  
      IC  CAA   CBA    CGA   O2A       .0000     .00  180.00     .00    .0000  
      IC  C3C   C1C    *C2C  CMC       .0000     .00  180.00     .00    .0000  
      IC  CBA   O1A    *CGA  O2A       .0000     .00  180.00     .00    .0000  
      IC  CBD   O1D    *CGD  O2D       .0000     .00  180.00     .00    .0000  
 
 
 END {HEME}
 !-----------------------------------------------------------
 
 PRESIDUE NTER ! standard N-terminus 
 
 GROUP  ! use in generate statement 
      MODI ATOM  +N     TYPE=NH3  CHARGE=    -.9000  END ! 
      ADD  ATOM  +HT1   TYPE=HC   CHARGE=     .2000  END !         HT1 
      ADD  ATOM  +HT2   TYPE=HC   CHARGE=     .2000  END !        / 
      ADD  ATOM  +HT3   TYPE=HC   CHARGE=     .2000  END ! --CA--N--HT2 
      MODI ATOM  +CA    TYPE=CT1  CHARGE=     .2000  END !   |    \ 
      MODI ATOM  +HA    TYPE=HB   CHARGE=     .1000  END !   HA    HT3 
      DELE ATOM  +HN       END  
 !END GROUP
 
 ADD  BOND  +HT1   +N      
 ADD  BOND  +HT2   +N      
 ADD  BOND  +HT3   +N      
 
 ADD ANGLe +HT1  +N    +HT2
 ADD ANGLe +HT1  +N    +HT3
 ADD ANGLe +HT2  +N    +HT3
 ADD ANGLe +HT2  +N    +CA
 ADD ANGLe +HT1  +N    +CA
 ADD ANGLe +HT3  +N    +CA

 ADD DIHEdral +HT1  +N    +CA   +C
 ADD DIHEdral +HT2  +N    +CA   +C
 ADD DIHEdral +HT3  +N    +CA   +C
 ADD DIHEdral +HT1  +N    +CA   +CB
 ADD DIHEdral +HT2  +N    +CA   +CB
 ADD DIHEdral +HT3  +N    +CA   +CB
 ADD DIHEdral +HT1  +N    +CA   +HA
 ADD DIHEdral +HT2  +N    +CA   +HA
 ADD DIHEdral +HT3  +N    +CA   +HA

 ADD  DONOR  +HT1   +N      
 ADD  DONOR  +HT2   +N      
 ADD  DONOR  +HT3   +N      
 
 ADD  IC  +HT1   +N      +CA    +C         .0000     .00  180.00     .00    .0000  
 ADD  IC  +HT2   +CA     *+N    +HT1       .0000     .00  120.00     .00    .0000  
 ADD  IC  +HT3   +CA     *+N    +HT2       .0000     .00  120.00     .00    .0000  
 
 END {NTER}
 !-----------------------------------------------------------
 
 PRESIDUE GLYP ! Glycine N-terminus 
 
 GROUP   ! use in generate statement 
      MODI ATOM  +N     TYPE=NH3  CHARGE=    -.3000  END ! 
      ADD  ATOM  +HT1   TYPE=HC   CHARGE=     .3300  END !   HA1   HT1 
      ADD  ATOM  +HT2   TYPE=HC   CHARGE=     .3300  END !   |    / 
      ADD  ATOM  +HT3   TYPE=HC   CHARGE=     .3300  END ! --CA--N--HT2 
      MODI ATOM  +CA    TYPE=CT2  CHARGE=     .1300  END !   |    \ 
      MODI ATOM  +HA1   TYPE=HB   CHARGE=     .0900  END !   HA2   HT3 
      MODI ATOM  +HA2   TYPE=HB   CHARGE=     .0900  END ! 
      DELE ATOM  +HN       END  
 !END GROUP
 
 ADD  BOND  +HT1   +N      
 ADD  BOND  +HT2   +N      
 ADD  BOND  +HT3   +N      

 ADD  ANGLe +HT1  +N    +HT2
 ADD  ANGLe +HT1  +N    +HT3
 ADD  ANGLe +HT2  +N    +HT3
 ADD  ANGLe +HT1  +N    +CA
 ADD  ANGLe +HT2  +N    +CA
 ADD  ANGLe +HT3  +N    +CA

 ADD  DIHEdral +HT1  +N    +CA   +C
 ADD  DIHEdral +HT2  +N    +CA   +C
 ADD  DIHEdral +HT3  +N    +CA   +C
 ADD  DIHEdral +HT1  +N    +CA   +HA1
 ADD  DIHEdral +HT2  +N    +CA   +HA1
 ADD  DIHEdral +HT3  +N    +CA   +HA1
 ADD  DIHEdral +HT1  +N    +CA   +HA2
 ADD  DIHEdral +HT2  +N    +CA   +HA2
 ADD  DIHEdral +HT3  +N    +CA   +HA2

 ADD  DONOR  +HT1   +N      
 ADD  DONOR  +HT2   +N      
 ADD  DONOR  +HT3   +N      
 
 ADD  IC  +HT1   +N     +CA    +C         .0000     .00  180.00     .00    .0000  
 ADD  IC  +HT2   +CA    *+N    +HT1       .0000     .00  120.00     .00    .0000  
 ADD  IC  +HT3   +CA    *+N    +HT2       .0000     .00  120.00     .00    .0000  
 
 
 END {GLYP}
 !-----------------------------------------------------------
 
 PRESIDUE PROP ! Proline N-Terminal 
 
 GROUP  ! use in generate statement 
      MODI ATOM  +N     TYPE=NP   CHARGE=    -.0700  END !   HA 
      ADD  ATOM  +HN1   TYPE=HC   CHARGE=     .2400  END !   | 
      ADD  ATOM  +HN2   TYPE=HC   CHARGE=     .2400  END !  -CA   HN1 
      MODI ATOM  +CD    TYPE=CP3  CHARGE=     .1600  END !  /  \ / 
      MODI ATOM  +CA    TYPE=CP1  CHARGE=     .1600  END !      N 
      MODI ATOM  +HA    TYPE=HB   CHARGE=     .0900  END !     / \ 
      MODI ATOM  +HD1   TYPE=HA   CHARGE=     .0900  END !  -CD   HN2 
      MODI ATOM  +HD2   TYPE=HA   CHARGE=     .0900  END !   | \ 
 !END GROUP
 
 ADD  BOND  +HN1   +N      
 ADD  BOND  +HN2   +N    
 
 ADD ANGLe +HN1  +N    +HN2
 ADD ANGLe +HN1  +N    +CA
 ADD ANGLe +HN2  +N    +CA
 ADD ANGLe +HN1  +N    +CD
 ADD ANGLe +HN2  +N    +CD

 DELE DIHEdral +N   +CA    +C    +O

 ADD DIHEdral +HN1  +N    +CA   +C
 ADD DIHEdral +HN2  +N    +CA   +C
 ADD DIHEdral +HN1  +N    +CA   +HA
 ADD DIHEdral +HN2  +N    +CA   +HA
 ADD DIHEdral +HN1  +N    +CA   +CB
 ADD DIHEdral +HN2  +N    +CA   +CB
 ADD DIHEdral +HN1  +N    +CD   +CG
 ADD DIHEdral +HN2  +N    +CD   +CG
 ADD DIHEdral +HN1  +N    +CD   +HD1
 ADD DIHEdral +HN2  +N    +CD   +HD1
 ADD DIHEdral +HN1  +N    +CD   +HD2
 ADD DIHEdral +HN2  +N    +CD   +HD2

 ADD  DONOR  +HN1   +N      
 ADD  DONOR  +HN2   +N      
 
 ADD  IC  +HN1   +CA     *+N    +CD        .0000     .00  120.00     .00    .0000  
 ADD  IC  +HN2   +CA     *+N    +HN1       .0000     .00  120.00     .00    .0000  
 
 
 END {PROP}
 !-----------------------------------------------------------
 ! 08-16-96 added CY-N-CA, CY-N-HN, and CAY-CY-N as angles - were missing
 PRESIDUE ACE  ! acetylated N-terminus 
 
 GROUP	! use in generate statement 
      ADD  ATOM  +CAY   TYPE=CT3  CHARGE=    -.2700  END ! 
      ADD  ATOM  +HY1   TYPE=HA   CHARGE=     .0900  END ! HY1 HY2 HY3 
      ADD  ATOM  +HY2   TYPE=HA   CHARGE=     .0900  END !    \ | / 
      ADD  ATOM  +HY3   TYPE=HA   CHARGE=     .0900  END !     CAY 
 GROUP	   					         !      | 
      ADD  ATOM  +CY    TYPE=C    CHARGE=     .5100  END !      CY=OY 
      ADD  ATOM  +OY    TYPE=O    CHARGE=    -.5100  END !      | 
      ! 
 !END GROUP
 
 ADD  BOND  +CY    +CAY    
 ADD  BOND  +OY    +CY     
 ADD  BOND  +CY    +N      
 ADD  BOND  +CAY   +HY1    
 ADD  BOND  +CAY   +HY2    
 ADD  BOND  +CAY   +HY3    

 ADD  ANGLe +HY1  +CAY    +CY
 ADD  ANGLe +HY2  +CAY    +CY
 ADD  ANGLe +HY2  +CAY    +HY1
 ADD  ANGLe +HY3  +CAY    +CY
 ADD  ANGLe +HY3  +CAY    +HY1
 ADD  ANGLe +HY3  +CAY    +HY2
 ADD  ANGLE +CAY  +CY     +N
 ADD  ANGLe +CAY  +CY     +OY
 ADD  ANGLe +OY   +CY     +N
 ADD  ANGLe +CY   +N      +CA
 ADD  ANGLe +CY   +N      +HN

 ADD  DIHEdral +CAY   +CY    +N    +CA
 ADD  DIHEdral +CAY   +CY    +N    +HN
 ADD  DIHEdral +HY1   +CAY   +CY   +N
 ADD  DIHEdral +HY2   +CAY   +CY   +N
 ADD  DIHEdral +HY3   +CAY   +CY   +N
 ADD  DIHEdral +HY1   +CAY   +CY   +OY
 ADD  DIHEdral +HY2   +CAY   +CY   +OY
 ADD  DIHEdral +HY3   +CAY   +CY   +OY
 ADD  DIHEdral +CY    +N     +CA   +C
 ADD  DIHEdral +CY    +N     +CA   +HA
 ADD  DIHEdral +CY    +N     +CA   +CB
 ADD  DIHEdral +OY    +CY    +N    +CA
 ADD  DIHEdral +OY    +CY    +N    +HN
 
 ADD  IMPROPER  +CY    +CAY   +N     +OY     
 ADD  IMPROPER  +N     +CY    +CA    +HN     
 
 ADD  ACCEPTOR  +OY    +CY     
 
 ADD  IC  +C    +CA    +N    +CY       .0000     .00  180.00     .00    .0000  
 ADD  IC  +CY   +CA    *+N   +HN       .0000     .00  180.00     .00    .0000  
 ADD  IC  +CAY  +CY    +N    +CA       .0000     .00  180.00     .00    .0000  
 ADD  IC  +N    +CAY   *+CY  +OY       .0000     .00  180.00     .00    .0000  
 ADD  IC  +OY   +CY    +CAY  +HY1      .0000     .00  180.00     .00    .0000  
 ADD  IC  +OY   +CY    +CAY  +HY2      .0000     .00   60.00     .00    .0000  
 ADD  IC  +OY   +CY    +CAY  +HY3      .0000     .00  -60.00     .00    .0000  
 
 
 END {ACE }
 !-----------------------------------------------------------
 
 PRESIDUE ACP  ! acetylated N-terminus for proline 
 ! 08-16-96 added CY-N-CA, CY-N-D, CAY-CY-N and OY-CY-N as angles - were missing
 
 GROUP  ! use in generate statement 
      ADD  ATOM  +CAY   TYPE=CT3  CHARGE=    -.2700  END ! 
      ADD  ATOM  +HY1   TYPE=HA   CHARGE=     .0900  END ! HY1 HY2 HY3 
      ADD  ATOM  +HY2   TYPE=HA   CHARGE=     .0900  END !    \ | / 
      ADD  ATOM  +HY3   TYPE=HA   CHARGE=     .0900  END !     CAY 
 GROUP						         !      | 
      ADD  ATOM  +CY    TYPE=C    CHARGE=     .5100  END !      CY=OY 
      ADD  ATOM  +OY    TYPE=O    CHARGE=    -.5100  END !      | 
      ! 
 !END GROUP
 
 ADD  BOND  +CY    +CAY    
 ADD  BOND  +OY    +CY     
 ADD  BOND  +CY    +N      
 ADD  BOND  +CAY   +HY1    
 ADD  BOND  +CAY   +HY2    
 ADD  BOND  +CAY   +HY3    

 ADD  ANGLe +HY1  +CAY    +CY
 ADD  ANGLe +HY2  +CAY    +CY
 ADD  ANGLe +HY2  +CAY    +HY1
 ADD  ANGLe +HY3  +CAY    +CY
 ADD  ANGLe +HY3  +CAY    +HY1
 ADD  ANGLe +HY3  +CAY    +HY2
 ADD  ANGLE +CAY  +CY     +N
 ADD  ANGLe +CAY  +CY     +OY
 ADD  ANGLe +OY   +CY     +N
 ADD  ANGLe +CY   +N      +CA
 ADD  ANGLe +CY   +N      +CD

 ADD  DIHEdral +CAY   +CY    +N    +CA
 ADD  DIHEdral +CAY   +CY    +N    +CA  ! multiple dihedral
 ADD  DIHEdral +CAY   +CY    +N    +CD
 ADD  DIHEdral +CAY   +CY    +N    +CD  ! multiple dihedral
 ADD  DIHEdral +HY1   +CAY   +CY   +N
 ADD  DIHEdral +HY2   +CAY   +CY   +N
 ADD  DIHEdral +HY3   +CAY   +CY   +N
 ADD  DIHEdral +HY1   +CAY   +CY   +OY
 ADD  DIHEdral +HY2   +CAY   +CY   +OY
 ADD  DIHEdral +HY3   +CAY   +CY   +OY
 ADD  DIHEdral +CY    +N     +CA   +C
 ADD  DIHEdral +CY    +N     +CA   +CB
 ADD  DIHEdral +CY    +N     +CA   +HA
 ADD  DIHEdral +CY    +N     +CD   +CG
 ADD  DIHEdral +CY    +N     +CD   +HD1
 ADD  DIHEdral +CY    +N     +CD   +HD2
 ADD  DIHEdral +OY    +CY    +N    +CA
 ADD  DIHEdral +OY    +CY    +N    +CA  ! multiple dihedral
 ADD  DIHEdral +OY    +CY    +N    +CD
 ADD  DIHEdral +OY    +CY    +N    +CD  ! multiple dihedral
 
 ADD  IMPROPER  +CY    +CAY   +N     +OY     
 ADD  IMPROPER  +N     +CY    +CA    +CD     
 
 ADD  ACCEPTOR  +OY    +CY     
 
 ADD  IC  +CY    +N      +CA    +C         .0000     .00  -60.00     .00    .0000  
 ADD  IC  +CY    +CA     *N     +CD        .0000     .00  180.00     .00    .0000  
 ADD  IC  +CAY   +CY     +N     +CA        .0000     .00  180.00     .00    .0000  
 ADD  IC  +N     +CAY    *CY    +OY        .0000     .00  180.00     .00    .0000  
 ADD  IC  +OY    +CY     +CAY   +HY1       .0000     .00  180.00     .00    .0000  
 ADD  IC  +OY    +CY     +CAY   +HY2       .0000     .00   60.00     .00    .0000  
 ADD  IC  +OY    +CY     +CAY   +HY3       .0000     .00  -60.00     .00    .0000  
 
 
 END {ACP }
 !-----------------------------------------------------------
 
 PRESIDUE CTER ! standard C-terminus 
 
 GROUP  ! use in generate statement 
      MODI ATOM  -C     TYPE=CC   CHARGE=    1.0000  END !   OT2 
      ADD  ATOM  -OT1   TYPE=OC   CHARGE=    -.5000  END !  // 
      ADD  ATOM  -OT2   TYPE=OC   CHARGE=    -.5000  END ! -C 
      DELE ATOM  -O        END 				 !  \\ 
 !END GROUP						 !   OT1 
 
 ADD  BOND  -C     -OT1    
 ADD  BOND  -C     -OT2  

 ADD ANGLe -CA   -C   -OT1
 ADD ANGLe -CA   -C   -OT2
 ADD ANGLe -OT1  -C   -OT2

 ADD DIHEdral -N    -CA    -C   -OT2
 ADD DIHEdral -N    -CA    -C   -OT1
 ADD DIHEdral -CB   -CA    -C   -OT2
 ADD DIHEdral -CB   -CA    -C   -OT1
 ADD DIHEdral -HA   -CA    -C   -OT2
 ADD DIHEdral -HA   -CA    -C   -OT1

 ADD  IMPROPER  -OT1   -CA    -OT2   -C      
 
 ADD  ACCEPTOR  -OT1   -C      
 ADD  ACCEPTOR  -OT2   -C      
 
 ADD  IC  -N     -CA     -C     -OT2       .0000     .00  180.00     .00    .0000  
 ADD  IC  -OT2   -CA     *-C    -OT1       .0000     .00  180.00     .00    .0000  
 
 
 END {CTER }
 !-----------------------------------------------------------
 
 PRESIDUE CT1  ! methylated C-terminus from methyl acetate 
 
 GROUP  ! use in generate statement 
      MODI ATOM  -N     TYPE=NH1  CHARGE=    -.4700  END ! 
      MODI ATOM  -HN    TYPE=H    CHARGE=     .3100  END !          OT1 
      MODI ATOM  -CA    TYPE=CT1  CHARGE=     .1700  END !     |   // 
      MODI ATOM  -HA    TYPE=HB   CHARGE=     .0900  END ! -N--CA--C       HT1 
      MODI ATOM  -C     TYPE=CD   CHARGE=     .6300  END !  |  |    \      / 
      ADD  ATOM  -OT1   TYPE=OB   CHARGE=    -.5200  END ! HN  HA    OT2--CT--HT2 
      ADD  ATOM  -OT2   TYPE=OS   CHARGE=    -.3400  END !                 \ 
      ADD  ATOM  -CT    TYPE=CT3  CHARGE=    -.1400  END !                 HT3 
      ADD  ATOM  -HT1   TYPE=HA   CHARGE=     .0900  END ! 
      ADD  ATOM  -HT2   TYPE=HA   CHARGE=     .0900  END ! 
      ADD  ATOM  -HT3   TYPE=HA   CHARGE=     .0900  END ! 
      DELE ATOM  -O        END  
 !END GROUP
 
 ADD  BOND  -C     -OT1    
 ADD  BOND  -C     -OT2    
 ADD  BOND  -OT2   -CT     
 ADD  BOND  -CT    -HT1    
 ADD  BOND  -CT    -HT2    
 ADD  BOND  -CT    -HT3  

 ADD  ANGLe -CA   -C   -OT1
 ADD  ANGLe -CA   -C   -OT2
 ADD  ANGLe -OT1  -C   -OT2
 ADD  ANGLe -C    -OT2 -CT
 ADD  ANGLe -HT1  -CT  -OT2
 ADD  ANGLe -HT2  -CT  -OT2
 ADD  ANGLe -HT3  -CT  -OT2
 ADD  ANGLe -HT1  -CT  -HT2
 ADD  ANGLe -HT1  -CT  -HT3
 ADD  ANGLe -HT2  -CT  -HT3

 ADD  DIHEDRAL  -N     -CA    -C     -OT1
 ADD  DIHEDRAL  -N     -CA    -C     -OT2
 ADD  DIHEDRAL  -CB    -CA    -C     -OT1
 ADD  DIHEDRAL  -CB    -CA    -C     -OT2
 ADD  DIHEDRAL  -HA    -CA    -C     -OT1
 ADD  DIHEDRAL  -HA    -CA    -C     -OT2
 ADD  DIHEDRAL  -CA    -C     -OT2   -CT     
 ADD  DIHEDRAL  -OT1   -C     -OT2   -CT  ! multiple dihedral
 ADD  DIHEDRAL  -OT1   -C     -OT2   -CT  ! multiple dihedral
 ADD  DIHEDRAL  -HT1   -CT    -OT2   -C     
 ADD  DIHEDRAL  -HT2   -CT    -OT2   -C     
 ADD  DIHEDRAL  -HT3   -CT    -OT2   -C     
  
 ADD  IMPROPER  -OT1   -CA    -OT2   -C      
 
 ADD  ACCEPTOR  -OT1   -C      
 ADD  ACCEPTOR  -OT2   -C      
 
 ADD  IC  -N     -CA     -C     -OT2       .0000     .00  180.00     .00    .0000  
 ADD  IC  -OT2   -CA     *-C    -OT1       .0000     .00  180.00     .00    .0000  
 ADD  IC  -CA    -C      -OT2   -CT        .0000     .00  180.00     .00    .0000  
 ADD  IC  -C     -OT2    -CT    -HT1       .0000     .00     .00     .00    .0000  
 ADD  IC  -C     -OT2    -CT    -HT2       .0000     .00  120.00     .00    .0000  
 ADD  IC  -C     -OT2    -CT    -HT3       .0000     .00  240.00     .00    .0000  
 
 
 END {CT1 }
 !-----------------------------------------------------------
 
 PRESIDUE CT2  ! amidated C-terminus 
 
 GROUP  ! use in generate or patch statement 
      MODI ATOM  -C     TYPE=CC   CHARGE=     .5500  END ! 
      MODI ATOM  -O     TYPE=O    CHARGE=    -.5500  END !     | 
 GROUP						         !   O=C 
      ADD  ATOM  -NT    TYPE=NH2  CHARGE=    -.6200  END !     | 
      ADD  ATOM  -HT1   TYPE=H    CHARGE=     .3200  END !     NT 
      ADD  ATOM  -HT2   TYPE=H    CHARGE=     .3000  END !    / \ 
 !END GROUP						 !  HT1 HT2 
 
 ADD  BOND  -C     -NT    
 ADD  BOND  -NT    -HT1    
 ADD  BOND  -NT    -HT2   ! 

 ADD ANGLe -CA   -C   -NT
 ADD ANGLe -O    -C   -NT
 ADD ANGLe -C    -NT  -HT1   
 ADD ANGLe -C    -NT  -HT2   
 ADD ANGLe -HT1  -NT  -HT2   
     
 ADD DIHEdral -N   -CA    -C   -NT
 ADD DIHEdral -N   -CA    -C   -NT   ! multiple dihedral
 ADD DIHEdral -HA  -CA    -C   -NT
 ADD DIHEdral -HA  -CA    -C   -NT   ! multiple dihedral
 ADD DIHEdral -CB  -CA    -C   -NT
 ADD DIHEdral -CB  -CA    -C   -NT   ! multiple dihedral
 ADD DIHEdral -O   -C     -NT  -HT1  !  (HT1 is cis to O)
 ADD DIHEdral -O   -C     -NT  -HT2
 ADD DIHEdral -CA  -C     -NT  -HT1
 ADD DIHEdral -CA  -C     -NT  -HT2
 
 ADD  IMPROPER  -C     -NT    -CA    -O      
 ADD  IMPROPER  -C     -CA    -NT    -O      
 ADD  IMPROPER  -NT    -C     -HT1   -HT2    
 ADD  IMPROPER  -NT    -C     -HT2   -HT1    
 
 ADD  DONOR  -HT1   -NT     
 ADD  DONOR  -HT2   -NT     
 
 ADD  IC  -N     -CA     -C     -O         .0000     .00  180.00     .00    .0000  
 ADD  IC  -NT    -CA     *-C    -O         .0000     .00  180.00     .00    .0000  
 ADD  IC  -CA    -C      -NT    -HT1       .0000     .00  180.00     .00    .0000  
 ADD  IC  -HT1   -C      *-NT   -HT2       .0000     .00  180.00     .00    .0000  
 
 
 END {CT2 }
 !-----------------------------------------------------------
 
 PRESIDUE CT3  ! N-Methylamide C-terminus 
 
 GROUP   ! use in generate or patch statement 
      MODI ATOM  -C     TYPE=C    CHARGE=     .5100  END ! 
      MODI ATOM  -O     TYPE=O    CHARGE=    -.5100  END !      | 
 GROUP 						         !      C=O 
      ADD  ATOM  -NT    TYPE=NH1  CHARGE=    -.4700  END !      | 
      ADD  ATOM  -HNT   TYPE=H    CHARGE=     .3100  END !      NT-HNT 
      ADD  ATOM  -CAT   TYPE=CT3  CHARGE=    -.1100  END !      | 
      ADD  ATOM  -HT1   TYPE=HA   CHARGE=     .0900  END ! HT1-CAT-HT3 
      ADD  ATOM  -HT2   TYPE=HA   CHARGE=     .0900  END !      | 
      ADD  ATOM  -HT3   TYPE=HA   CHARGE=     .0900  END !     HT2 
      ! 
 !END GROUP
 
 ADD  BOND  -C     -NT     
 ADD  BOND  -NT    -HNT    
 ADD  BOND  -NT    -CAT    
 ADD  BOND  -CAT   -HT1    
 ADD  BOND  -CAT   -HT2    
 ADD  BOND  -CAT   -HT3    

 ADD  ANGLe -CA   -C   -NT
 ADD  ANGLe -O    -C   -NT
 ADD  ANGLe -HNT  -NT  -C
 ADD  ANGLe -HNT  -NT  -CAT
 ADD  ANGLe -C    -NT  -CAT
 ADD  ANGLe -HT1  -CAT -NT
 ADD  ANGLe -HT2  -CAT -NT
 ADD  ANGLe -HT3  -CAT -NT
 ADD  ANGLe -HT1  -CAT -HT2
 ADD  ANGLe -HT1  -CAT -HT3
 ADD  ANGLe -HT2  -CAT -HT3
 
 ADD  DIHEDRAL  -N     -CA    -C     -NT
 ADD  DIHEDRAL  -N     -CA    -C     -NT   ! multiple dihedral
 ADD  DIHEDRAL  -HA    -CA    -C     -NT
 ADD  DIHEDRAL  -HA    -CA    -C     -NT   ! multiple dihedral
 ADD  DIHEDRAL  -CB    -CA    -C     -NT
 ADD  DIHEDRAL  -CB    -CA    -C     -NT   ! multiple dihedral
 ADD  DIHEDRAL  -CA    -C     -NT    -CAT     
 ADD  DIHEDRAL  -O     -C     -NT    -CAT  
 ADD  DIHEDRAL  -CA    -C     -NT    -HNT
 ADD  DIHEDRAL  -O     -C     -NT    -HNT
 ADD  DIHEDRAL  -HT1   -CAT   -NT    -C     
 ADD  DIHEDRAL  -HT2   -CAT   -NT    -C     
 ADD  DIHEDRAL  -HT3   -CAT   -NT    -C
 ADD  DIHEDRAL  -HT1   -CAT   -NT    -HNT     
 ADD  DIHEDRAL  -HT2   -CAT   -NT    -HNT   
 ADD  DIHEDRAL  -HT3   -CAT   -NT    -HNT
 
 ADD  IMPROPER  -NT    -C     -CAT   -HNT    
 ADD  IMPROPER  -C     -CA    -NT    -O      
 
 ADD  DONOR  -HNT   -NT     
 
 ADD  IC  -N     -CA     -C     -O         .0000     .00  180.00     .00    .0000  
 ADD  IC  -NT    -CA     *-C    -O         .0000     .00  180.00     .00    .0000  
 ADD  IC  -C     -CAT    *-NT   -HNT       .0000     .00  180.00     .00    .0000  
 ADD  IC  -CA    -C      -NT    -CAT       .0000     .00  180.00     .00    .0000  
 ADD  IC  -C     -NT     -CAT   -HT1       .0000     .00   60.00     .00    .0000  
 ADD  IC  -C     -NT     -CAT   -HT2       .0000     .00  180.00     .00    .0000  
 ADD  IC  -C     -NT     -CAT   -HT3       .0000     .00  -60.00     .00    .0000  
 
 
 END {CT3 }
 !-----------------------------------------------------------
 
 PRESIDUE ASPP ! patch for protonated aspartic acid, proton on od2 
 
 GROUP	! via acetic acid, use in a patch statement 
      MODI ATOM  CB    TYPE=CT2  CHARGE=    -.2100  END ! 
      MODI ATOM  HB1   TYPE=HA   CHARGE=     .0900  END ! HB1    OD1 
      MODI ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !  |    // 
      MODI ATOM  CG    TYPE=CD   CHARGE=     .7500  END ! -CB--CG 
      MODI ATOM  OD1   TYPE=OB   CHARGE=    -.5500  END !  |     \ 
      MODI ATOM  OD2   TYPE=OH1  CHARGE=    -.6100  END ! HB2     OD2-HD2 
      ADD  ATOM  HD2   TYPE=H    CHARGE=     .4400  END ! 
 !END GROUP
 
 ADD  BOND  OD2   HD2    
 
 ADD  ANGLE  HD2   OD2   CG     
 
 ADD  DIHEDRAL  HD2   OD2   CG    OD1    
 ADD  DIHEDRAL  HD2   OD2   CG    CB     
 
 ADD  DONOR  HD2   OD2    
 
 ADD  IC  HD2   OD2    CG    OD1       .0000     .00     .00     .00    .0000  
 
 
 END {ASPP}
 !-----------------------------------------------------------
 
 PRESIDUE GLUP ! patch for protonated glutamic acid, proton on oe2 
 
 GROUP  ! via acetic acid, use in a patch statement 
      MODI ATOM  CG    TYPE=CT2  CHARGE=    -.2100  END ! 
      MODI ATOM  HG1   TYPE=HA   CHARGE=     .0900  END ! HG1    OE1 
      MODI ATOM  HG2   TYPE=HA   CHARGE=     .0900  END !  |    // 
      MODI ATOM  CD    TYPE=CD   CHARGE=     .7500  END ! -CG--CD 
      MODI ATOM  OE1   TYPE=OB   CHARGE=    -.5500  END !  |     \ 
      MODI ATOM  OE2   TYPE=OH1  CHARGE=    -.6100  END ! HG2     OE2-HE2 
      ADD  ATOM  HE2   TYPE=H    CHARGE=     .4400  END ! 
 !END GROUP
 
 ADD  BOND  OE2   HE2    
 
 ADD  ANGLE  HE2   OE2   CD     
 
 ADD  DIHEDRAL  HE2   OE2   CD    OE1    
 ADD  DIHEDRAL  HE2   OE2   CD    CG     
 
 ADD  DONOR  HE2   OE2    
 
 ADD  IC  HE2   OE2    CD    OE1       .0000     .00     .00     .00    .0000  
 

 END {GLUP}
 !-----------------------------------------------------------
 
 PRESIDUE LINK ! linkage for IMAGES or for joining segments 
	       ! 1 refers to previous (N terminal) 
	       ! 2 refers to next (C terminal) 
	       ! use in a patch st        
 
 
 ADD  BOND  1C    2N     
 
 ADD  ANGLE  1C    2N    2CA    
 ADD  ANGLE  1C    2N    2HN    
 ADD  ANGLE  1CA   1C    2N     
 ADD  ANGLE  1O    1C    2N     
 
 ADD  DIHEDRAL  1C    2N    2CA   2C     
 ADD  DIHEDRAL  1C    2N    2CA   2HA    
 ADD  DIHEDRAL  1C    2N    2CA   2CB    
 ADD  DIHEDRAL  1CA   1C    2N    2HN    
 ADD  DIHEDRAL  1CA   1C    2N    2CA    ! for multiple dihedral
 ADD  DIHEDRAL  1CA   1C    2N    2CA    ! for multiple dihedral
 ADD  DIHEDRAL  1HA   1CA   1C    2N     
 ADD  DIHEDRAL  1N    1CA   1C    2N     
 ADD  DIHEDRAL  1CB   1CA   1C    2N     
 ADD  DIHEDRAL  1O    1C    2N    2HN    
 ADD  DIHEDRAL  1O    1C    2N    2CA    
 
 ADD  IMPROPER  2N    1C    2CA   2HN    
 ADD  IMPROPER  1C    1CA   2N    1O     
 
 ADD  IC  1N    1CA    1C    2N        .0000     .00  180.00     .00    .0000  
 ADD  IC  2N    1CA    *1C   1O        .0000     .00  180.00     .00    .0000  
 ADD  IC  1CA   1C     2N    2CA       .0000     .00  180.00     .00    .0000  
 ADD  IC  1C    2N     2CA   2C        .0000     .00  180.00     .00    .0000  
 ADD  IC  1C    2CA    *2N   2HN       .0000     .00  180.00     .00    .0000  
 
 
 END {LINK}
 !-----------------------------------------------------------
 
 PRESIDUE DISU ! patch for disulfides. Patch must be 1-CYS and 2-CYS.

 GROUP ! use in a patch statement
      delete    atom 1HG        end
      MODIFY ATOM  1CB   TYPE=CT2  CHARGE=    -.1000  END !
      MODIFY ATOM  1SG   TYPE=SM   CHARGE=    -.0800  END !           2SG--2CB--
 GROUP                                                  !          /
      MODIFY ATOM  2CB   TYPE=CT2  CHARGE=    -.1000  END !
      MODIFY ATOM  2SG   TYPE=SM   CHARGE=    -.0800  END ! -1CB--1SG
      delete    atom 2HG        end
 !END GROUP
 
 ADD  BOND  1SG   2SG    
 
 ADD  ANGLE  1CB   1SG   2SG    
 ADD  ANGLE  1SG   2SG   2CB    
 
 ADD  DIHEDRAL  1HB1  1CB   1SG   2SG    
 ADD  DIHEDRAL  1HB2  1CB   1SG   2SG    
 ADD  DIHEDRAL  2HB1  2CB   2SG   1SG    
 ADD  DIHEDRAL  2HB2  2CB   2SG   1SG    
 ADD  DIHEDRAL  1CA   1CB   1SG   2SG    
 ADD  DIHEDRAL  1SG   2SG   2CB   2CA    
 ADD  DIHEDRAL  1CB   1SG   2SG   2CB    ! multiple dihedral
 ADD  DIHEDRAL  1CB   1SG   2SG   2CB    ! multiple dihedral
 ADD  DIHEDRAL  1CB   1SG   2SG   2CB    ! multiple dihedral
 
 ADD  IC  1CA   1CB    1SG   2SG       .0000     .00  180.00     .00    .0000  
 ADD  IC  1CB   1SG    2SG   2CB       .0000     .00   90.00     .00    .0000  
 ADD  IC  1SG   2SG    2CB   2CA       .0000     .00  180.00     .00    .0000  
 
 
 END {DISU}
 !-----------------------------------------------------------
 
 PRESIDUE HS2  ! Patch for neutral His, move proton from ND1 to NE2 
 
 GROUP	! use in a patch statement 
      MODI ATOM  CE1   TYPE=CPH2 CHARGE=     .2500  END !                 HE1 
      MODI ATOM  HE1   TYPE=HR1  CHARGE=     .1300  END !                 / 
      MODI ATOM  ND1   TYPE=NR2  CHARGE=    -.7000  END !   HB1    ND1--CE1 
      MODI ATOM  CG    TYPE=CPH1 CHARGE=     .2200  END !   |     /      | 
      MODI ATOM  CB    TYPE=CT2  CHARGE=    -.0800  END !  -CB--CG       | 
      MODI ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !   |     \      | 
      MODI ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !   HB2    CD2--NE2 
 GROUP						        !           |     \ 
      MODI ATOM  NE2   TYPE=NR1  CHARGE=    -.3600  END !          HD2    HE2 
      ADD  ATOM  HE2   TYPE=H    CHARGE=     .3200  END  
      MODI ATOM  CD2   TYPE=CPH1 CHARGE=    -.0500  END  
      MODI ATOM  HD2   TYPE=HR3  CHARGE=     .0900  END  
      DELE ATOM  HD1      END  
 !END GROUP
 
 ADD  BOND  NE2   HE2    
 
 ADD  ANGLE  CD2   NE2   HE2    
 ADD  ANGLE  CE1   NE2   HE2    
 
 ADD  DIHEDRAL  HE2   NE2   CE1   HE1    
 ADD  DIHEDRAL  HE2   NE2   CE1   ND1    
 ADD  DIHEDRAL  HE2   NE2   CD2   HD2    
 ADD  DIHEDRAL  HE2   NE2   CD2   CG     
 
 ADD  IMPROPER  NE2   CD2   CE1   HE2    
 ADD  IMPROPER  NE2   CE1   CD2   HE2    
 
 ADD  DONOR  HE2   NE2    
 
 DELE ACCEPTOR  NE2   NONE   
 ADD  ACCEPTOR  ND1   NONE   
 
 ADD  IC  CE1   CD2   *NE2   HE2       .0000     .00  180.00     .00    .0000  
 
 END {HS2 }
 !-----------------------------------------------------------
 
 RESIDUE O2   ! O2 ligand for heme 
 
 GROUP  
       ATOM  O1    TYPE=OM   CHARGE=     .0210  END  
       ATOM  O2    TYPE=OM   CHARGE=    -.0210  END  
 !END GROUP
 
      BOND  O1    O2     
 
 
 END {O2  }
 !-----------------------------------------------------------
 
 RESIDUE CO   ! CO ligand for heme 
 
 GROUP  
       ATOM  C     TYPE=CM   CHARGE=     .0210  END  
       ATOM  O     TYPE=OM   CHARGE=    -.0210  END  
 !END GROUP
 
      BOND  C     O      
 
 
 END {CO  }
 !-----------------------------------------------------------
 
 PRESIDUE FHEM ! FIX UP THE HEME BY DELETING UNWANTED AUTOGENERATED ANGLES 
	       ! unliganded heme patch 
 
 
 DELE ANGLE  1NA   1FE   1NC    
 DELE ANGLE  1NB   1FE   1ND    
 
 
 END {FHEM}
 !-----------------------------------------------------------
 
 PRESIDUE PHEM ! Patch for HEME to His link. 
	       ! Patch residues must be 1-HIS, and 2-HEME. 
 
 
 ADD  BOND  1NE2  2FE    
 
 DELE ANGLE  2NA   2FE   2NC    
 DELE ANGLE  2NB   2FE   2ND    

 ADD  ANGLE  1CD2  1NE2  2FE    
 ADD  ANGLE  1CE1  1NE2  2FE    
 ADD  ANGLE  1NE2  2FE   2NA    
 ADD  ANGLE  1NE2  2FE   2NB    
 ADD  ANGLE  1NE2  2FE   2NC    
 ADD  ANGLE  1NE2  2FE   2ND    


 ADD  IC  1CD2  1NE2   2FE   2NA       .0000     .00     .00     .00    .0000  
 ADD  IC  1CD2  1NE2   2FE   2NB       .0000     .00     .00     .00    .0000  
 ADD  IC  1CD2  1NE2   2FE   2NB       .0000     .00     .00     .00    .0000  
 ADD  IC  1CD2  1NE2   2FE   2NB       .0000     .00     .00     .00    .0000  
 ADD  IC  1CE1  1NE2   2FE   2NA       .0000     .00     .00     .00    .0000  
 
 
 END {PHEM}
 !-----------------------------------------------------------
 
 PRESIDUE PLO2 ! Patch residue for Heme ligand.  Residues must be 
	       ! 1-O2 , 2-HEME , and 3-HIS 
	       ! O1 of the oxygen is bonded to the iron. 
 
 
 ADD  BOND  1O1   2FE    
 
 ADD  ANGLE  1O2   1O1   2FE    
 ADD  ANGLE  1O1   2FE   2NA    
 ADD  ANGLE  1O1   2FE   2NB    
 ADD  ANGLE  1O1   2FE   2NC    
 ADD  ANGLE  1O1   2FE   2ND    
 
 ADD  DIHEDRAL  1O2   1O1   2FE   2NA    
 
 ADD  IC  1O2   1O1    2FE   2NA       .0000     .00     .00     .00    .0000  
 ADD  IC  1O2   1O1    2FE   3NE2      .0000     .00     .00     .00    .0000  
 
 
 END {PLO2}
 !-----------------------------------------------------------
 
 PRESIDUE PLIG ! Patch residue for Heme ligand. Residues must be, 
	       ! 1-CO , 2-HEME , and 3-HIS 
 
 
 ADD  BOND  1C    2FE    
 
 ADD  ANGLE  1C    2FE   3NE2   
 ADD  ANGLE  1O    1C    2FE    
 ADD  ANGLE  1C    2FE   2NA    
 ADD  ANGLE  1C    2FE   2NB    
 ADD  ANGLE  1C    2FE   2NC    
 ADD  ANGLE  1C    2FE   2ND    
 
 ADD  IC  1O    1C     2FE   2NA       .0000     .00     .00     .00    .0000  
 ADD  IC  1O    1C     2FE   3NE2      .0000     .00     .00     .00    .0000  
 
 
 END {PLIG}
 !-----------------------------------------------------------
 
 PRESIDUE PLWA ! Patch residue for Heme ligand. Residues must be, 
	       ! 1-TIP3 , 2-HEME , and 3-HIS 
 
 
 ADD  BOND  1OH2  2FE    
 
 ADD  ANGLE  1OH2  2FE   3NE2   
 ADD  ANGLE  1H1   1OH2  FE     
 ADD  ANGLE  1H2   1OH2  FE     
 ADD  ANGLE  1O    1OH2  2FE    
 ADD  ANGLE  1OH2  2FE   2NA    
 ADD  ANGLE  1OH2  2FE   2NB    
 ADD  ANGLE  1OH2  2FE   2NC    
 ADD  ANGLE  1OH2  2FE   2ND    
 
 ADD  IC  1H1   1OH2   2FE   2NA       .0000     .00     .00     .00    .0000  
 ADD  IC  1H1   1OH2   2FE   3NE2      .0000     .00     .00     .00    .0000  
 
 
 END {PLWA}
 !------------------------------------------------------------------
 
 PRESidue PEPT { PEPTide bond link, for all 
                amino acids ...*(-)     (+)*...
                                 \ PEPT /

                except the  *(-) - (+)PRO link        }

      ADD BOND -C     +N

      ADD ANGLE -CA -C +N
      ADD ANGLE -O  -C +N
      ADD ANGLE -C  +N +CA
      ADD ANGLE -C  +N +HN

      ADD DIHEdral  -C    +N    +CA   +C
!      ADD DIHEdral  -C    +N    +CA   +HA
!      ADD DIHEdral  -C    +N    +CA   +CB
!      ADD DIHEdral  -HA   -CA   -C    +N
      ADD DIHEdral  -N    -CA   -C    +N
!      ADD DIHEdral  -CB   -CA   -C    +N
      ADD DIHEdral  -CA   -C    +N    +HN
      ADD DIHEdral  -CA   -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -CA   -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -O    -C    +N    +HN
      ADD DIHEdral  -O    -C    +N    +CA

      ADD IMPR +N    -C    +CA   +HN
      ADD IMPR -C    -CA   +N    -O

 ADD  IC  -N    -CA    -C    +N      .0000     .00  180.00     .00    .0000
 ADD  IC  -CA    -C    +N    +CA     .0000     .00  180.00     .00    .0000
 ADD  IC  -C    +N     +CA   +C      .0000     .00  180.00     .00    .0000
 ADD  IC  +N    -CA    *-C   -O      .0000     .00  180.00     .00    .0000
 ADD  IC  -C    +CA    *+N   +HN     .0000     .00  180.00     .00    .0000

 END {PEPT}

 !------------------------------------------------------------------
 
 PRESidue PPG1 { for ...*(-) - (+) GLY LINK
               same as PEPT except replacement HA,CB with HA1 HA2
               at the (+) positions, required for proper dihedral setup }


      ADD BOND -C     +N

      ADD ANGLE -CA -C +N
      ADD ANGLE -O  -C +N
      ADD ANGLE -C  +N +CA
      ADD ANGLE -C  +N +HN

!      ADD DIHEdral  -C    +N    +CA   +C
!      ADD DIHEdral  -C    +N    +CA   +HA1
!      ADD DIHEdral  -C    +N    +CA   +HA2
!      ADD DIHEdral  -HA   -CA   -C    +N
!      ADD DIHEdral  -N    -CA   -C    +N
!      ADD DIHEdral  -CB   -CA   -C    +N
!      ADD DIHEdral  -CA   -C    +N    +HN
!      ADD DIHEdral  -CA   -C    +N    +CA    ! multiple dihedral
!      ADD DIHEdral  -CA   -C    +N    +CA    ! multiple dihedral
!      ADD DIHEdral  -O    -C    +N    +HN
!      ADD DIHEdral  -O    -C    +N    +CA

      ADD IMPR +N    -C    +CA   +HN
      ADD IMPR -C    -CA   +N    -O

 ADD  IC  -N    -CA    -C    +N      .0000     .00  180.00     .00    .0000
 ADD  IC  -CA    -C    +N    +CA     .0000     .00  180.00     .00    .0000
 ADD  IC  -C    +N     +CA   +C      .0000     .00  180.00     .00    .0000
 ADD  IC  +N    -CA    *-C   -O      .0000     .00  180.00     .00    .0000
 ADD  IC  -C    +CA    *+N   +HN     .0000     .00  180.00     .00    .0000

 END {PPG1}
 !------------------------------------------------------------------
 
 PRESidue PPG2 { for ... GLY(-) - (+)* LINK
               same as PEPT except replacement HA,CB with HA1 HA2
               at the (-) positions, required for proper dihedral setup }


      ADD BOND -C     +N

      ADD ANGLE -CA -C +N
      ADD ANGLE -O  -C +N
      ADD ANGLE -C  +N +CA
      ADD ANGLE -C  +N +HN

      ADD DIHEdral  -C    +N    +CA   +C
      ADD DIHEdral  -C    +N    +CA   +HA
      ADD DIHEdral  -C    +N    +CA   +CB
      ADD DIHEdral  -HA1  -CA   -C    +N
      ADD DIHEdral  -N    -CA   -C    +N
      ADD DIHEdral  -HA2  -CA   -C    +N
      ADD DIHEdral  -CA   -C    +N    +HN
      ADD DIHEdral  -CA   -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -CA   -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -O    -C    +N    +HN
      ADD DIHEdral  -O    -C    +N    +CA

      ADD IMPR +N    -C    +CA   +HN
      ADD IMPR -C    -CA   +N    -O

 ADD  IC  -N    -CA    -C    +N      .0000     .00  180.00     .00    .0000
 ADD  IC  -CA    -C    +N    +CA     .0000     .00  180.00     .00    .0000
 ADD  IC  -C    +N     +CA   +C      .0000     .00  180.00     .00    .0000
 ADD  IC  +N    -CA    *-C   -O      .0000     .00  180.00     .00    .0000
 ADD  IC  -C    +CA    *+N   +HN     .0000     .00  180.00     .00    .0000

 END {PPG2}
 !------------------------------------------------------------------
 
 PRESidue PEPP  { for  ...*(-) - (+)PRO  link
               same as PEPT except replacement H by CD
               and improper +N +CA +CD -C              }

      ADD BOND -C +N 

      ADD ANGLE -CA -C +N
      ADD ANGLE -O  -C +N
      ADD ANGLE -C  +N +CA
      ADD ANGLE -C  +N +CD

      ADD DIHEdral  -C    +N    +CA   +C
      ADD DIHEdral  -C    +N    +CA   +HA
      ADD DIHEdral  -C    +N    +CA   +CB
      ADD DIHEdral  -HA   -CA   -C    +N     
      ADD DIHEdral  -N    -CA   -C    +N     
      ADD DIHEdral  -CB   -CA   -C    +N     
      ADD DIHEdral  -CA   -C    +N    +CD    ! multiple dihedral
      ADD DIHEdral  -CA   -C    +N    +CD    ! multiple dihedral
      ADD DIHEdral  -CA   -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -CA   -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -O    -C    +N    +CD    ! multiple dihedral
      ADD DIHEdral  -O    -C    +N    +CD    ! multiple dihedral
      ADD DIHEdral  -O    -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -O    -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -C    +N    +CD   +CG    ! for proline
      ADD DIHEdral  -C    +N    +CD   +HD1   ! for proline
      ADD DIHEdral  -C    +N    +CD   +HD2   ! for proline

      ADD IMPRoper  -C -CA +N -O  {planar -C}
! CDS: force constant is zero!
!      ADD IMPRoper  +N -C +CA +CD  !!! +N  +CA +CD -C  {planar +N} MODIFIED

 ADD  IC  -N    -CA    -C    +N      .0000     .00  180.00     .00    .0000
 ADD  IC  -CA    -C    +N    +CA     .0000     .00  180.00     .00    .0000
 ADD  IC  -C    +N     +CA   +C      .0000     .00  -80.00     .00    .0000
 ADD  IC  +N    -CA    *-C   -O      .0000     .00  180.00     .00    .0000
 ADD  IC  -C    +CA    *+N   +CD     .0000     .00  180.00     .00    .0000

 END {PEPP}
 !------------------------------------------------------------------

 PRESidue PPGG { for ... GLY(-) - (+) GLY LINK
               same as PEPT except replacement HA,CB with HA1 HA2
               at the (+) and (-) positions, 
	       required for proper dihedral setup } 
						{ PDA 5-5-94 }


      ADD BOND -C     +N

      ADD ANGLE -CA -C +N
      ADD ANGLE -O  -C +N
      ADD ANGLE -C  +N +CA
      ADD ANGLE -C  +N +HN

      ADD DIHEdral  -C    +N    +CA   +C
      ADD DIHEdral  -C    +N    +CA   +HA1
      ADD DIHEdral  -C    +N    +CA   +HA2
      ADD DIHEdral  -HA1  -CA   -C    +N
      ADD DIHEdral  -N    -CA   -C    +N
      ADD DIHEdral  -HA2  -CA   -C    +N
      ADD DIHEdral  -CA   -C    +N    +HN
      ADD DIHEdral  -CA   -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -CA   -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -O    -C    +N    +HN
      ADD DIHEdral  -O    -C    +N    +CA

      ADD IMPR +N    -C    +CA   +HN
      ADD IMPR -C    -CA   +N    -O

 ADD  IC  -N    -CA    -C    +N      .0000     .00  180.00     .00    .0000
 ADD  IC  -CA    -C    +N    +CA     .0000     .00  180.00     .00    .0000
 ADD  IC  -C    +N     +CA   +C      .0000     .00  180.00     .00    .0000
 ADD  IC  +N    -CA    *-C   -O      .0000     .00  180.00     .00    .0000
 ADD  IC  -C    +CA    *+N   +HN     .0000     .00  180.00     .00    .0000

 END {PPGG}
 !------------------------------------------------------------------

 PRESidue PPGP  { for  ... GLY(-) - (+)PRO  link
               same as PEPT except replacement +H by +CD
                            and -HA by -HA1 and -CB by -HA2
               and improper +N +CA +CD -C              }
							 { PDA 5-5-94 }

      ADD BOND -C +N 

      ADD ANGLE -CA -C +N
      ADD ANGLE -O  -C +N
      ADD ANGLE -C  +N +CA
      ADD ANGLE -C  +N +CD

      ADD DIHEdral  -C    +N    +CA   +C
      ADD DIHEdral  -C    +N    +CA   +HA
      ADD DIHEdral  -C    +N    +CA   +CB
      ADD DIHEdral  -HA1  -CA   -C    +N     
      ADD DIHEdral  -N    -CA   -C    +N     
      ADD DIHEdral  -HA2  -CA   -C    +N     
      ADD DIHEdral  -CA   -C    +N    +CD    ! multiple dihedral
      ADD DIHEdral  -CA   -C    +N    +CD    ! multiple dihedral
      ADD DIHEdral  -CA   -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -CA   -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -O    -C    +N    +CD    ! multiple dihedral
      ADD DIHEdral  -O    -C    +N    +CD    ! multiple dihedral
      ADD DIHEdral  -O    -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -O    -C    +N    +CA    ! multiple dihedral
      ADD DIHEdral  -C    +N    +CD   +CG    ! for proline
      ADD DIHEdral  -C    +N    +CD   +HD1   ! for proline
      ADD DIHEdral  -C    +N    +CD   +HD2   ! for proline

      ADD IMPRoper  -C -CA +N -O  {planar -C}
      ADD IMPRoper  +N -C +CA +CD  !!! +N  +CA +CD -C  {planar +N} MODIFIED

 ADD  IC  -N    -CA    -C    +N      .0000     .00  180.00     .00    .0000
 ADD  IC  -CA    -C    +N    +CA     .0000     .00  180.00     .00    .0000
 ADD  IC  -C    +N     +CA   +C      .0000     .00  -80.00     .00    .0000
 ADD  IC  +N    -CA    *-C   -O      .0000     .00  180.00     .00    .0000
 ADD  IC  -C    +CA    *+N   +CD     .0000     .00  180.00     .00    .0000

 END {PPGP}
 !------------------------------------------------------------------

 PRESidue LIG1  { linkage for cyclic peptide
                  1 refers to the C terminus which is a glycine
                  2 refers to the N terminus
                  use in a patch statement, perform initial
                  generation using first NONE last NONE }

 ADD  BOND  1C 2N

 ADD  ANGLE  1C    2N    2CA
 ADD  ANGLE  1CA   1C    2N
 ADD  ANGLE  1O    1C    2N
 ADD  ANGLE  1C    2N    2HN

 ADD  DIHEDRAL  1C    2N    2CA   2C    
 ADD  DIHEDRAL  1C    2N    2CA   2HA
 ADD  DIHEDRAL  1C    2N    2CA   2CB   
 ADD  DIHEDRAL  1HA1  1CA   1C    2N
 ADD  DIHEDRAL  1N    1CA   1C    2N    
 ADD  DIHEDRAL  1HA2  1CA   1C    2N
 ADD  DIHEDRAL  1CA   1C    2N    2HN   
 ADD  DIHEDRAL  1CA   1C    2N    2CA   ! multiple dihedral
 ADD  DIHEDRAL  1CA   1C    2N    2CA   ! multiple dihedral
 ADD  DIHEDRAL  1O    1C    2N    2HN
 ADD  DIHEDRAL  1O    1C    2N    2CA   
 ADD  DIHEDRAL  1CA   1C    2N    2CA   

 ADD  IMPROPER  2N    1C    2CA   2HN
 ADD  IMPROPER  1C    1CA   2N    1O

 ADD  IC  1N    1CA    1C    2N        .0000     .00  180.00     .00   .0000
 ADD  IC  2N    1CA    *1C   1O        .0000     .00  180.00     .00   .0000
 ADD  IC  1CA   1C     2N    2CA       .0000     .00  180.00     .00   .0000
 ADD  IC  1C    2N     2CA   2C        .0000     .00  180.00     .00   .0000
 ADD  IC  1C    2CA   *2N    2HN       .0000     .00  180.00     .00   .0000

 END {LIG1}
 !-----------------------------------------------------------

 PRESidue LIG2  { linkage for cyclic peptide
                  1 refers to the C terminus
                  2 refers to the N terminus which is a glycine
                  use in a patch statement, perform initial
                  generation using first NONE last NONE }

 ADD  BOND  1C 2N

 ADD  ANGLE  1C    2N    2CA
 ADD  ANGLE  1CA   1C    2N
 ADD  ANGLE  1O    1C    2N
 ADD  ANGLE  1C    2N    2HN

 ADD  DIHEDRAL  1C    2N    2CA   2C    
 ADD  DIHEDRAL  1C    2N    2CA   2HA1
 ADD  DIHEDRAL  1C    2N    2CA   2HA2  
 ADD  DIHEDRAL  1HA   1CA   1C    2N
 ADD  DIHEDRAL  1N    1CA   1C    2N    
 ADD  DIHEDRAL  1CB   1CA   1C    2N
 ADD  DIHEDRAL  1CA   1C    2N    2HN   
 ADD  DIHEDRAL  1CA   1C    2N    2CA   ! multiple dihedral
 ADD  DIHEDRAL  1CA   1C    2N    2CA   ! multiple dihedral
 ADD  DIHEDRAL  1O    1C    2N    2HN
 ADD  DIHEDRAL  1O    1C    2N    2CA   
 ADD  DIHEDRAL  1CA   1C    2N    2CA   

 ADD  IMPROPER  2N    1C    2CA   2HN
 ADD  IMPROPER  1C    1CA   2N    1O

 ADD  IC  1N    1CA    1C    2N        .0000     .00  180.00     .00   .0000
 ADD  IC  2N    1CA    *1C   1O        .0000     .00  180.00     .00   .0000
 ADD  IC  1CA   1C     2N    2CA       .0000     .00  180.00     .00   .0000
 ADD  IC  1C    2N     2CA   2C        .0000     .00  180.00     .00   .0000
 ADD  IC  1C    2CA    *2N   2HN       .0000     .00  180.00     .00   .0000

 END {LIG2}
 !-----------------------------------------------------------

 PRESidue LIG3  { linkage for cyclic peptide
                  1 refers to the C terminus which is a glycine
                  2 refers to the N terminus which is a glycine
                  use in a patch statement, perform initial
                  generation using first NONE last NONE }

 ADD  BOND  1C 2N

 ADD  ANGLE  1C    2N    2CA
 ADD  ANGLE  1CA   1C    2N
 ADD  ANGLE  1O    1C    2N
 ADD  ANGLE  1C    2N    2HN

 ADD  DIHEDRAL  1C    2N    2CA   2C    
 ADD  DIHEDRAL  1C    2N    2CA   2HA1
 ADD  DIHEDRAL  1C    2N    2CA   2HA2  
 ADD  DIHEDRAL  1HA1  1CA   1C    2N
 ADD  DIHEDRAL  1N    1CA   1C    2N    
 ADD  DIHEDRAL  1HA2  1CA   1C    2N
 ADD  DIHEDRAL  1CA   1C    2N    2HN   
 ADD  DIHEDRAL  1CA   1C    2N    2CA   ! multiple dihedral
 ADD  DIHEDRAL  1CA   1C    2N    2CA   ! multiple dihedral
 ADD  DIHEDRAL  1O    1C    2N    2HN
 ADD  DIHEDRAL  1O    1C    2N    2CA   
 ADD  DIHEDRAL  1CA   1C    2N    2CA   

 ADD  IMPROPER  2N    1C    2CA   2HN
 ADD  IMPROPER  1C    1CA   2N    1O

 ADD  IC  1N    1CA    1C    2N        .0000     .00  180.00     .00   .0000
 ADD  IC  2N    1CA    *1C   1O        .0000     .00  180.00     .00   .0000
 ADD  IC  1CA   1C     2N    2CA       .0000     .00  180.00     .00   .0000
 ADD  IC  1C    2N     2CA   2C        .0000     .00  180.00     .00   .0000
 ADD  IC  1C    2CA    *2N   2HN       .0000     .00  180.00     .00   .0000

 END {LIG3}
 !-----------------------------------------------------------
 
 RESIDUE SOD  ! Sodium Ion 
 
 GROUP  
       ATOM  SOD   TYPE=SOD  CHARGE=    1.0000  END  
 !END GROUP
 
 
 END {SOD }
 !-----------------------------------------------------------
 
 RESIDUE CL   ! Chloride Anion 
 
 GROUP  
       ATOM  CL1   TYPE=CLA  CHARGE=   -1.0000  END  
 !END GROUP
 
 
 END {CL  }

! linkage statement comment must begin with !LINK
! order is significant
!LINK PPGP HEAD - GLY TAIL + PRO END
!LINK PPGG HEAD - GLY TAIL + GLY END
!LINK PPG1 HEAD - *   TAIL + GLY END
!LINK PPG2 HEAD - GLY TAIL + *   END
!LINK PEPP HEAD - *   TAIL + PRO END
!LINK PEPT  HEAD - *     TAIL + *       END

 SET ECHO=TRUE END 
