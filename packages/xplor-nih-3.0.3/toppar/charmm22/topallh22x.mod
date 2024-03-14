 REMARKS Charmm topology file for miscellaneous molecules v22 b4
 REMARKS FILENAME="topallh22x.mod"
 
 SET ECHO=FALSE END 

 {>> All-hydrogen topology for small model compounds used in the <<
  >> development of the CHARMM22 protein all-hydrogen parameters <<
  >>>>>>>>>>>>>>>>>>>>>>> June 1992<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
  >>>>>>>> Direct comments to Alexander D. MacKerell Jr. <<<<<<<<<<
  >>>>>>> 410-706-7442 or bitnet: alex,tammy.harvard.edu <<<<<<<<<<
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
 MASS HS     1.0080 ! thiol hydrogen 
 MASS C     12.0110 ! polar C 
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
 MASS CC    12.0110 ! carbonyl C for sidechains asn,asp,gln,glu 
 MASS CD    12.0110 ! carbonyl C for none amides, asp,glu,cter 
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
 !-----------------------------------------------------------
 
 RESIDUE ACET ! acetate, K. Kuczera 
 
 GROUP  
       ATOM  C1    TYPE=CT3  CHARGE=    -.3700  END ! 
       ATOM  C2    TYPE=CC   CHARGE=     .6200  END !     H1    O1 
       ATOM  H1    TYPE=HA   CHARGE=     .0900  END !     |    // 
       ATOM  H2    TYPE=HA   CHARGE=     .0900  END ! H2--C1--C2 
       ATOM  H3    TYPE=HA   CHARGE=     .0900  END !     |    \\ 
       ATOM  O1    TYPE=OC   CHARGE=    -.7600  END !     H3    O2 
       ATOM  O2    TYPE=OC   CHARGE=    -.7600  END ! 
 !END GROUP
 
      BOND  C1    H1     
      BOND  C1    H2     
      BOND  C1    H3     
      BOND  C1    C2     
      BOND  C2    O1     
      BOND  C2    O2     
 
      IMPROPER  O1    C1    O2    C2     
 
      IC  H1    C1     C2    O1        .0000     .00     .00     .00    .0000  
      IC  H2    C1     C2    O1        .0000     .00     .00     .00    .0000  
      IC  H3    C1     C2    O1        .0000     .00     .00     .00    .0000  
      IC  H1    C1     C2    O2        .0000     .00     .00     .00    .0000  
      IC  H2    C1     C2    O2        .0000     .00     .00     .00    .0000  
      IC  H3    C1     C2    O2        .0000     .00     .00     .00    .0000  
      IC  H1    C1     H2    H3        .0000     .00     .00     .00    .0000  
      IC  H2    C1     H3    H1        .0000     .00     .00     .00    .0000  
      IC  H3    C1     H1    H2        .0000     .00     .00     .00    .0000  
      IC  O1    O2     *C2   C1        .0000     .00     .00     .00    .0000  
 
 
 END {ACET}
 !-----------------------------------------------------------
 
 RESIDUE GUAN ! guandinium, K. Kuczera 
 
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .6400  END !      H11  H12 
       ATOM  N1    TYPE=NC2  CHARGE=    -.8000  END !        \  / 
       ATOM  H11   TYPE=HC   CHARGE=     .4600  END !         N1 
       ATOM  H12   TYPE=HC   CHARGE=     .4600  END !         | 
       ATOM  N2    TYPE=NC2  CHARGE=    -.8000  END !         C 
       ATOM  H21   TYPE=HC   CHARGE=     .4600  END !        / \ 
       ATOM  H22   TYPE=HC   CHARGE=     .4600  END !  H21-N2   N3-H31 
       ATOM  N3    TYPE=NC2  CHARGE=    -.8000  END !       |   | 
       ATOM  H31   TYPE=HC   CHARGE=     .4600  END !     H22   H32 
       ATOM  H32   TYPE=HC   CHARGE=     .4600  END ! 
 !END GROUP
 
      BOND  C     N1     
      BOND  C     N2     
      BOND  C     N3     
      BOND  N1    H11    
      BOND  N1    H12    
      BOND  N2    H21    
      BOND  N2    H22    
      BOND  N3    H31    
      BOND  N3    H32    
 
      IMPROPER  C     N2    N1    N3     
 
      IC  H11   N1     C     N2        .0000     .00     .00     .00    .0000  
      IC  H11   N1     C     N3        .0000     .00  180.00     .00    .0000  
      IC  H12   N1     C     N3        .0000     .00     .00     .00    .0000  
      IC  H21   N2     C     N1        .0000     .00     .00     .00    .0000  
      IC  H22   N2     C     N3        .0000     .00     .00     .00    .0000  
      IC  H31   N3     C     N1        .0000     .00     .00     .00    .0000  
      IC  H32   N3     C     N2        .0000     .00     .00     .00    .0000  
 
 
 END {GUAN}
 !-----------------------------------------------------------
 
 RESIDUE MGUA ! methyl-guanidinium 
 
 GROUP  
       ATOM  C     TYPE=C    CHARGE=     .6400  END !      H11  H12 
       ATOM  N1    TYPE=NC2  CHARGE=    -.8000  END !        \  / 
       ATOM  H11   TYPE=HC   CHARGE=     .4600  END !         N1 
       ATOM  H12   TYPE=HC   CHARGE=     .4600  END !         | 
       ATOM  N2    TYPE=NC2  CHARGE=    -.8000  END !         C       HC1 
       ATOM  H21   TYPE=HC   CHARGE=     .4600  END !        / \     / 
       ATOM  H22   TYPE=HC   CHARGE=     .4600  END !  H21-N2   N3--C2-HC2 
       ATOM  N3    TYPE=NC2  CHARGE=    -.7000  END !       |   |    \ 
       ATOM  H31   TYPE=HC   CHARGE=     .4400  END !     H22   H31   HC3 
       ATOM  C2    TYPE=CT3  CHARGE=     .1100  END ! 
       ATOM  HC1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HC2   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HC3   TYPE=HA   CHARGE=     .0900  END  
 !END GROUP
 
      BOND  C     N1     
      BOND  C     N2     
      BOND  C     N3     
      BOND  N1    H11    
      BOND  N1    H12    
      BOND  N2    H21    
      BOND  N2    H22    
      BOND  N3    H31    
      BOND  N3    C2     
      BOND  C2    HC1    
      BOND  C2    HC2    
      BOND  C2    HC3    
 
      IMPROPER  C     N2    N1    N3    !IMPH  N3 C2 C H31 
 
      IC  H11   N1     C     N2        .0000     .00     .00     .00    .0000  
      IC  H11   N1     C     N3        .0000     .00  180.00     .00    .0000  
      IC  H12   N1     C     N3        .0000     .00     .00     .00    .0000  
      IC  H21   N2     C     N1        .0000     .00     .00     .00    .0000  
      IC  H22   N2     C     N3        .0000     .00     .00     .00    .0000  
      IC  H31   N3     C     N1        .0000     .00     .00     .00    .0000  
      IC  C2    N3     C     N1        .0000     .00  180.00     .00    .0000  
      IC  C     N3     C2    HC1       .0000     .00     .00     .00    .0000  
      IC  C     N3     C2    HC2       .0000     .00  120.00     .00    .0000  
      IC  C     N3     C2    HC3       .0000     .00  240.00     .00    .0000  
 
 
 END {MGUA}
 !-----------------------------------------------------------
 
 RESIDUE MAMM ! methylammonium, K. Kuczera 
 
 GROUP  
       ATOM  CE    TYPE=CT3  CHARGE=     .1600  END ! 
       ATOM  NZ    TYPE=NH3  CHARGE=    -.3000  END ! HE1     HZ1 
       ATOM  HE1   TYPE=HA   CHARGE=     .0500  END !    \     | 
       ATOM  HE2   TYPE=HA   CHARGE=     .0500  END !HE2-CE---NZ--HZ2 
       ATOM  HE3   TYPE=HA   CHARGE=     .0500  END !    /     | 
       ATOM  HZ1   TYPE=HC   CHARGE=     .3300  END ! HE3     HZ3 
       ATOM  HZ2   TYPE=HC   CHARGE=     .3300  END ! 
       ATOM  HZ3   TYPE=HC   CHARGE=     .3300  END ! 
 !END GROUP
 
      BOND  CE    HE1    
      BOND  CE    HE2    
      BOND  CE    HE3    
      BOND  CE    NZ     
      BOND  NZ    HZ1    
      BOND  NZ    HZ2    
      BOND  NZ    HZ3    
 
      IC  HE1   CE     NZ    HZ1       .0000     .00     .00     .00    .0000  
      IC  HE1   CE     NZ    HZ2       .0000     .00  120.00     .00    .0000  
      IC  HE1   CE     NZ    HZ3       .0000     .00  240.00     .00    .0000  
      IC  HE2   CE     NZ    HZ1       .0000     .00  120.00     .00    .0000  
      IC  HE2   CE     NZ    HZ2       .0000     .00     .00     .00    .0000  
      IC  HE2   CE     NZ    HZ3       .0000     .00  240.00     .00    .0000  
      IC  HE3   CE     NZ    HZ1       .0000     .00  240.00     .00    .0000  
      IC  HE3   CE     NZ    HZ2       .0000     .00  120.00     .00    .0000  
      IC  HE3   CE     NZ    HZ3       .0000     .00     .00     .00    .0000  
      IC  CE    NZ     HZ1   HZ3       .0000     .00     .00     .00    .0000  
 
 
 END {MAMM}
 !-----------------------------------------------------------
 
 RESIDUE MAMI ! methylamine, adm jr. 
 
 GROUP  
       ATOM  N1    TYPE=NH2  CHARGE=    -.9600  END ! 
       ATOM  C1    TYPE=CT3  CHARGE=     .1300  END !  HC1       HN1 
       ATOM  HN1   TYPE=HC   CHARGE=     .3400  END !     \      / 
       ATOM  HC1   TYPE=HA   CHARGE=     .0500  END ! HC2-C1---NZ 
       ATOM  HC2   TYPE=HA   CHARGE=     .0500  END !     /      \ 
       ATOM  HC3   TYPE=HA   CHARGE=     .0500  END !  HC3       HN2 
       ATOM  HN2   TYPE=HC   CHARGE=     .3400  END ! 
 !END GROUP
 
      BOND  N1    C1     
      BOND  N1    HN1    
      BOND  N1    HN2    
      BOND  C1    HC1    
      BOND  C1    HC2    
      BOND  C1    HC3    
 
      IC  BLNK  HN1    N1    C1        .0000     .00  180.00     .00    .0000  
      IC  HN1   N1     C1    HC1       .0000     .00  180.00     .00    .0000  
      IC  HN1   N1     C1    HC2       .0000     .00   60.00     .00    .0000  
      IC  HN1   N1     C1    HC3       .0000     .00  300.00     .00    .0000  
      IC  HC1   C1     N1    HN2       .0000     .00     .00     .00    .0000  
 
 
 END {MAMI}
 !-----------------------------------------------------------
 
 RESIDUE ACEH ! acetic acid, ADM jr. 
 
 GROUP  
       ATOM  C2    TYPE=CT3  CHARGE=    -.3000  END ! 
       ATOM  C1    TYPE=CD   CHARGE=     .7500  END !   H21      O2 
       ATOM  H21   TYPE=HA   CHARGE=     .0900  END !      \    // 
       ATOM  H22   TYPE=HA   CHARGE=     .0900  END !  H22-C2--C1 
       ATOM  H23   TYPE=HA   CHARGE=     .0900  END !      /     \ 
       ATOM  O2    TYPE=OB   CHARGE=    -.5500  END !   H23       O1-HO1 
       ATOM  O1    TYPE=OH1  CHARGE=    -.6100  END ! 
       ATOM  HO1   TYPE=H    CHARGE=     .4400  END ! 
 !END GROUP
 
      BOND  C1    O1     
      BOND  O1    HO1    
      BOND  C1    O2     
      BOND  C1    C2     
      BOND  C2    H21    
      BOND  C2    H22    
      BOND  C2    H23    
 
      IMPROPER  O2    C2    O1    C1     
 
      DONOR  BLNK  HO1   ! O1 
 
      ACCEPTOR  O1    NONE   
      ACCEPTOR  O2    NONE   
 
      IC  O2    C1     C2    H21       .0000     .00     .00     .00    .0000  
      IC  HO1   O1     C1    O2        .0000     .00     .00     .00    .0000  
      IC  HO1   O1     C1    C2        .0000     .00  180.00     .00    .0000  
      IC  O1    C1     C2    H21       .0000     .00  180.00     .00    .0000  
      IC  O1    C1     C2    H22       .0000     .00   60.00     .00    .0000  
      IC  O1    C1     C2    H23       .0000     .00  -60.00     .00    .0000  
 
 
 END {ACEH}
 !-----------------------------------------------------------
 
 RESIDUE MAS  ! methylacetate 
 
 GROUP ! 
       ATOM  C1    TYPE=CT3  CHARGE=    -.1700  END !           H22 
       ATOM  C     TYPE=CD   CHARGE=     .6300  END !           | 
       ATOM  OM    TYPE=OS   CHARGE=    -.3400  END !       H21-C2-H23 
       ATOM  C2    TYPE=CT3  CHARGE=    -.1400  END !             \ 
       ATOM  O     TYPE=OB   CHARGE=    -.5200  END !             OM 
       ATOM  H11   TYPE=HA   CHARGE=     .0900  END !            / 
       ATOM  H12   TYPE=HA   CHARGE=     .0900  END !         O=C 
       ATOM  H13   TYPE=HA   CHARGE=     .0900  END !           | 
       ATOM  H21   TYPE=HA   CHARGE=     .0900  END !       H11-C1-H13 
       ATOM  H22   TYPE=HA   CHARGE=     .0900  END !           | 
       ATOM  H23   TYPE=HA   CHARGE=     .0900  END !           H12 
 !END GROUP
 
      BOND  C1    C      
      BOND  C     OM     
      BOND  C     O      
      BOND  OM    C2     
      BOND  C1    H11    
      BOND  C1    H12    
      BOND  C1    H13    
      BOND  C2    H21    
      BOND  C2    H22    
      BOND  C2    H23    
 
      IMPROPER  O     C1    OM    C

      ! internal coordinates from experiment for heavy atoms 

      IC  C1    C      OM    C2       1.5200  109.00  180.00  114.80   1.4370  
      IC  O     C      OM    C2       1.2000  125.90     .00  114.80   1.4370  
      IC  H11   C1     C     OM       1.1000  108.90  180.00  109.00   1.3340  
      IC  H12   C1     C     OM       1.1000  109.75   60.40  109.00   1.3340  
      IC  H13   C1     C     OM       1.1000  109.75  -60.40  109.00   1.3340  
      IC  H21   C2     OM    C        1.0788  109.94  180.00  114.80   1.3340  
      IC  H22   C2     OM    C        1.0802  110.50   60.50  114.80   1.3340  
      IC  H23   C2     OM    C        1.0802  110.50  -60.50  114.80   1.3340  
 
 
 END {MAS }
 !-----------------------------------------------------------
 
 RESIDUE MEOH ! methanol, adm jr. 
 
 GROUP ! order of atoms to match that used in ab initio 
       ATOM  CB    TYPE=CT3  CHARGE=    -.0400  END !  H11 
       ATOM  OG    TYPE=OH1  CHARGE=    -.6600  END !     \ 
       ATOM  HG1   TYPE=H    CHARGE=     .4300  END ! H12--C1--O1 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !     /      \ 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !  H13       HO1 
       ATOM  HB3   TYPE=HA   CHARGE=     .0900  END ! 
 !END GROUP
 
      BOND  CB    OG     
      BOND  OG    HG1    
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CB    HB3    
 
      DONOR  HG1   OG     
 
      ACCEPTOR  OG    NONE   
 
      IC  HG1   OG     CB    HB1       .0000     .00  180.00     .00    .0000  
      IC  HG1   OG     CB    HB2       .0000     .00  180.00     .00    .0000  
      IC  HG1   OG     CB    HB3       .0000     .00  180.00     .00    .0000  
      IC  OG    CB     HB1   HB2       .0000     .00  180.00     .00    .0000  
      IC  HG1   OG     CB    HB1       .0000     .00  180.00     .00    .0000  
 
 
 END {MEOH}
 !-----------------------------------------------------------
 
 RESIDUE METO ! methoxide, adm jr. 
 
 GROUP ! order of atoms to match that used in ab initio 
       ATOM  CB    TYPE=CT3  CHARGE=    -.4100  END !  HB1 
       ATOM  OG    TYPE=OC   CHARGE=    -.9200  END !     \ 
					            ! HB2--CB--OG 
       ATOM  HB1   TYPE=HA   CHARGE=     .1100  END !     / 
       ATOM  HB2   TYPE=HA   CHARGE=     .1100  END !  HB3 
       ATOM  HB3   TYPE=HA   CHARGE=     .1100  END ! 
 !END GROUP
 
      BOND  CB    OG     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CB    HB3    

      ACCEPTOR  OG    NONE   
 
      IC  HB2   OG     CB    HB1       .0000     .00  180.00     .00    .0000  
      IC  OG    HB2    *CB   HB1       .0000     .00  180.00     .00    .0000  
      IC  OG    CB     HB1   HB2       .0000     .00  180.00     .00    .0000  
 
 
 END {METO}
 !-----------------------------------------------------------
 
 RESIDUE ETOH ! Ethanol, adm jr. 
 
 GROUP  
       ATOM  C1    TYPE=CT2  CHARGE=     .0500  END !  H21  H11 H12 
       ATOM  O1    TYPE=OH1  CHARGE=    -.6600  END !     \   \  / 
       ATOM  HO1   TYPE=H    CHARGE=     .4300  END ! H22--C2--C1 
       ATOM  H11   TYPE=HA   CHARGE=     .0900  END !     /      \ 
       ATOM  H12   TYPE=HA   CHARGE=     .0900  END !  H23        O1--HO1 
 GROUP  
       ATOM  C2    TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  H21   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H22   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H23   TYPE=HA   CHARGE=     .0900  END  
 !END GROUP
 
      BOND  C1    C2     
      BOND  C1    O1     
      BOND  C1    H11    
      BOND  C1    H12    
      BOND  O1    HO1    
      BOND  C2    H21    
      BOND  C2    H22    
      BOND  C2    H23    
 
      DIHEDRAL C2  C1  O1  HO1 ! multiple dihedral
      DIHEDRAL C2  C1  O1  HO1 ! multiple dihedral
      DIHEDRAL C2  C1  O1  HO1 ! multiple dihedral

      DONOR  HO1   O1     
 
      ACCEPTOR  O1    NONE   
 
      IC  O1    C1     C2    H21       .0000     .00  180.00     .00    .0000  
      IC  O1    C1     C2    H22       .0000     .00   60.00     .00    .0000  
      IC  O1    C1     C2    H23       .0000     .00  300.00     .00    .0000  
      IC  H21   C2     C1    H11       .0000     .00  120.00     .00    .0000  
      IC  H21   C2     C1    H12       .0000     .00  240.00     .00    .0000  
      IC  C2    C1     O1    HO1       .0000     .00  180.00     .00    .0000 ! only for analysis 
      !IC HO1  O1   C1   H11   0.0000  0.0000 180.0000  0.0000  0.0000 
      !IC HO1  O1   C1   H12   0.0000  0.0000 180.0000  0.0000  0.0000
      !IC HO1  O1   C1   H13   0.0000  0.0000 180.0000  0.0000  0.0000
      !IC O1   C1   H11  H12   0.0000  0.0000 180.0000  0.0000  0.0000
 
 END {ETOH}
 !-----------------------------------------------------------
 
 RESIDUE ETO  ! Ethoxide, adm jr. 
 
 GROUP  
       ATOM  OG    TYPE=OC   CHARGE=    -.9200  END !  HA1  HB1 HB2 
       ATOM  CB    TYPE=CT2  CHARGE=    -.3000  END !     \   \  / 
       ATOM  CA    TYPE=CT3  CHARGE=    -.2700  END ! HA2--CA--CB 
       ATOM  HB1   TYPE=HA   CHARGE=     .1100  END !     /      \ 
       ATOM  HB2   TYPE=HA   CHARGE=     .1100  END !  HA3        OG 
      !GROUP 
       ATOM  HA1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HA2   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HA3   TYPE=HA   CHARGE=     .0900  END  
 !END GROUP
 
      BOND  CA    CB     
      BOND  CB    OG     
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CA    HA1    
      BOND  CA    HA2    
      BOND  CA    HA3    
 
      ACCEPTOR  OG    NONE   
 
      IC  CB    CA     HA1   HA2       .0000     .00  180.00     .00    .0000  
      IC  CA    OG     CB    HB1       .0000     .00  180.00     .00    .0000  
      IC  CA    OG     CB    HB2       .0000     .00  180.00     .00    .0000  
      IC  CA    OG     CB    HB3       .0000     .00  180.00     .00    .0000  
      IC  OG    CB     CA    HA1       .0000     .00  180.00     .00    .0000  
      IC  OG    CB     CA    HA2       .0000     .00   60.00     .00    .0000  
      IC  OG    CB     CA    HA3       .0000     .00  300.00     .00    .0000  
      IC  HA1   CA     CB    HB1       .0000     .00  120.00     .00    .0000  
      IC  HA1   CA     CB    HB2       .0000     .00  240.00     .00    .0000  
 
 
 END {ETO }
 !-----------------------------------------------------------
 
 RESIDUE PRO2 ! 2-proponal, adm jr. 
 
 GROUP  
       ATOM  C2    TYPE=CT1  CHARGE=     .1400  END !  H12  H13  H33 H32 
       ATOM  O2    TYPE=OH1  CHARGE=    -.6600  END !     \ /      \ / 
       ATOM  HO2   TYPE=H    CHARGE=     .4300  END ! H11--C1      C3--H31 
       ATOM  H21   TYPE=HA   CHARGE=     .0900  END !        \    / 
 GROUP						    !          C2 
       ATOM  C1    TYPE=CT3  CHARGE=    -.2700  END !        /    \ 
       ATOM  H11   TYPE=HA   CHARGE=     .0900  END !      O2     H21 
       ATOM  H12   TYPE=HA   CHARGE=     .0900  END !       | 
       ATOM  H13   TYPE=HA   CHARGE=     .0900  END !      HO2 
 GROUP ! 
       ATOM  C3    TYPE=CT3  CHARGE=    -.2700  END ! 
       ATOM  H31   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H32   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H33   TYPE=HA   CHARGE=     .0900  END  
 !END GROUP
 
      BOND  C1    C2     
      BOND  C2    C3     
      BOND  C2    O2     
      BOND  C2    H21    
      BOND  O2    HO2    
      BOND  C1    H11    
      BOND  C1    H12    
      BOND  C1    H13    
      BOND  C3    H31    
      BOND  C3    H32    
      BOND  C3    H33    

      DIHEDRAL C1  C2  O2  HO2 ! multiple dihedral
      DIHEDRAL C1  C2  O2  HO2 ! multiple dihedral
      DIHEDRAL C1  C2  O2  HO2 ! multiple dihedral
      DIHEDRAL C3  C2  O2  HO2 ! multiple dihedral
      DIHEDRAL C3  C2  O2  HO2 ! multiple dihedral
      DIHEDRAL C3  C2  O2  HO2 ! multiple dihedral

      DONOR  HO2   O2     
 
      ACCEPTOR  O2    NONE   
 
      IC  C1    C2     C3    H31       .0000     .00     .00     .00    .0000  
      IC  C1    C2     C3    H32       .0000     .00  120.00     .00    .0000  
      IC  C1    C2     C3    H33       .0000     .00  240.00     .00    .0000  
      IC  C3    C2     C1    H11       .0000     .00  180.00     .00    .0000  
      IC  C3    C2     C1    H12       .0000     .00   60.00     .00    .0000  
      IC  C3    C2     C1    H13       .0000     .00  300.00     .00    .0000  
      IC  H31   C3     C2    O2        .0000     .00  120.00     .00    .0000  
      IC  H31   C3     C2    H21       .0000     .00  240.00     .00    .0000  
      IC  C3    C2     O2    HO2       .0000     .00  180.00     .00    .0000
      ! only for analysis 
      !IC HO1  O1   C1   H11   0.0000  0.0000 180.0000  0.0000  0.0000 
      !IC HO1  O1   C1   H12   0.0000  0.0000 180.0000  0.0000  0.0000
      !IC HO1  O1   C1   H13   0.0000  0.0000 180.0000  0.0000  0.0000
      !IC O1   C1   H11  H12   0.0000  0.0000 180.0000  0.0000  0.0000
 
 END {PRO2}
 !-----------------------------------------------------------
 
 RESIDUE FORM ! formamide, adm jr. 
 
 GROUP  
       ATOM  HA    TYPE=HA   CHARGE=     .0800  END ! 
       ATOM  C     TYPE=CC   CHARGE=     .4200  END !  O     Hc 
       ATOM  N     TYPE=NH2  CHARGE=    -.6900  END !  \\   / 
       ATOM  HC    TYPE=H    CHARGE=     .3500  END !   C--N 
       ATOM  HT    TYPE=H    CHARGE=     .3500  END !  /    \ 
       ATOM  O     TYPE=O    CHARGE=    -.5100  END ! HA     Ht 
      ! 
 !END GROUP
 
      BOND  C     O      
      BOND  C     HA     
      BOND  C     N      
      BOND  N     HC     
      BOND  N     HT    
      !BOND  O   DUM 
 
      IMPROPER  C     HA    N     O      
      IMPROPER  C     N     HA    O
      ! GAS PHASE GEOMETRY 
      !IC  DUM O   C   N      1.00    90.0     0. 124.70 1.352 
      IC  O     C      N     HC       1.2190  124.70     .00  118.50   1.0016 ! variable 1 
      IC  HC    N      C     HA       1.0016  118.50  180.00  112.70   1.0980 ! variable 2 
      IC  O     C      N     HT       1.2190  124.70  180.00  120.00   1.0015 ! variable 3 
 
 
 END {FORM}
 !-----------------------------------------------------------
 
 RESIDUE ACEM ! acetamide, adm jr. 
	      ! the amide charges listed below are used in asparagine and glutamine 
	      ! if Hc and Ht are made equivalent use N                         
 
 GROUP  
       ATOM  CC    TYPE=CT3  CHARGE=    -.2700  END ! 
       ATOM  C     TYPE=CC   CHARGE=     .5500  END !  HC1           Ht 
       ATOM  N     TYPE=NH2  CHARGE=    -.6200  END !     \         / 
       ATOM  HC    TYPE=H    CHARGE=     .3200  END ! HC2--CC--C---N 
       ATOM  HT    TYPE=H    CHARGE=     .3000  END !     /    ||   \ 
       ATOM  O     TYPE=O    CHARGE=    -.5500  END !  HC3     O     Hc 
       ATOM  HC1   TYPE=HA   CHARGE=     .0900  END ! 
       ATOM  HC2   TYPE=HA   CHARGE=     .0900  END ! 
       ATOM  HC3   TYPE=HA   CHARGE=     .0900  END ! 
      !atom dum  dum    0.0  ! dummy for ic build 
 !END GROUP
 
      BOND  C     O      
      BOND  C     N      
      BOND  N     HC     
      BOND  N     HT     
      BOND  C     CC     
      BOND  CC    HC1    
      BOND  CC    HC2    
      BOND  CC    HC3   
      !BOND  O   DUM 
 
      IMPROPER  C     CC    N     O      
      IMPROPER  C     N     CC    O      
      IMPROPER  N     C     HC    HT     
      IMPROPER  N     C     HT    HC

      ! GAS PHASE GEOMETRY 
      IC  HC1   CC     C     N        1.0832  113.80  180.00  124.70   1.3523  
      IC  O     C      N     HC       1.2012  124.70     .00  118.63    .9929 ! variable 1 
      IC  O     C      N     HT       1.2012  124.70  180.00  120.92    .9960 ! variable 3 
      IC  HC    N      C     CC       1.0016  118.63  180.00  115.65   1.5150 ! variable 2 
      IC  O     C      CC    HC1      1.2012  123.27     .00  113.80   1.0832  
      IC  O     C      CC    HC2      1.2012  123.27  120.00  108.52   1.0836  
      IC  O     C      CC    HC3      1.2012  123.27  300.00  108.52   1.0836  
 
 
 END {ACEM}
 !-----------------------------------------------------------
 
 RESIDUE NMA  ! N-methylacetamide, Louis Kuchnir 
 
 GROUP  
       ATOM  CL    TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  HL1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HL2   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HL3   TYPE=HA   CHARGE=     .0900  END  
       ATOM  C     TYPE=C    CHARGE=     .5100  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
       ATOM  N     TYPE=NH1  CHARGE=    -.4700  END  
       ATOM  H     TYPE=H    CHARGE=     .3100  END  
       ATOM  CR    TYPE=CT3  CHARGE=    -.1100  END  
       ATOM  HR1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HR2   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HR3   TYPE=HA   CHARGE=     .0900  END  
 !END GROUP
 
      BOND  HL1   CL     
      BOND  HL2   CL     
      BOND  HL3   CL    !          N-Methylacetamide: 
      BOND  CL    C     !      HL1\       O            /HR1  
      BOND  C     N     !    HL2-- CL  -- C -- N -- CR --HR2  
      BOND  N     CR    !      HL3/            H       \HR3 
      BOND  C     O      
      BOND  N     H     
      BOND  HR1   CR     
      BOND  HR2   CR     
      BOND  HR3   CR 
   
      DIHEDRAL  CL    C     N     CR    ! multiple dihedral
      DIHEDRAL  CL    C     N     CR    ! multiple dihedral
 
      IMPROPER  N     C     CR    H      
      IMPROPER  C     CL    N     O      
 
      IC  O     C      N     H        1.2233  122.84  180.00  119.23    .9933 ! variable 1 
      IC  H     N      C     CL        .9933  119.23     .00  116.25   1.5118 ! variable 8 
      IC  O     C      N     CR       1.2233  122.84     .00  122.57   1.4488 ! variable 8 
      IC  N     C      CL    HL1      1.3418  116.25  180.00  109.30   1.1090 ! variable 2 
      IC  N     C      CL    HL2      1.3418  116.25   60.00  109.30   1.1090 ! variable 3 
      IC  N     C      CL    HL3      1.3418  116.25  300.00  109.30   1.1090 ! variable 4 
      IC  C     N      CR    HR1      1.3418  122.57  180.00  110.70   1.1130 ! variable 5 
      IC  C     N      CR    HR2      1.3418  122.57   60.00  110.70   1.1130 ! variable 6 
      IC  C     N      CR    HR3      1.3418  122.57  300.00  110.70   1.1130 ! variable 7 
 
 
 END {NMA }
 !-----------------------------------------------------------
 
 RESIDUE ALAD ! Alanine dipeptide, Louis Kuchnir 
 
 GROUP  
       ATOM  CL    TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  HL1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HL2   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HL3   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  CLP   TYPE=C    CHARGE=     .5100  END  
       ATOM  OL    TYPE=O    CHARGE=    -.5100  END  
 GROUP  
       ATOM  NL    TYPE=NH1  CHARGE=    -.4700  END  
       ATOM  HL    TYPE=H    CHARGE=     .3100  END  
       ATOM  CA    TYPE=CT1  CHARGE=     .0700  END  
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END  
 GROUP  
       ATOM  CB    TYPE=CT3  CHARGE=    -.2700  END !     HL1     OL           OR           HR1 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !       \     ||   HL  HA  ||   HR      / 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !        \    ||   |   |   ||   |      / 
       ATOM  HB3   TYPE=HA   CHARGE=     .0900  END !   HL2---CL--CLP--NL--CA--CRP--NR---CR---HR2 
 GROUP						    !        /             |               \ 
       ATOM  CRP   TYPE=C    CHARGE=     .5100  END !       /         HB1--CB--HB3          \ 
       ATOM  OR    TYPE=O    CHARGE=    -.5100  END !     HL3              |                HR3 
 GROUP						    !                     HB2 
       ATOM  NR    TYPE=NH1  CHARGE=    -.4700  END  
       ATOM  HR    TYPE=H    CHARGE=     .3100  END  
       ATOM  CR    TYPE=CT3  CHARGE=    -.1100  END  
       ATOM  HR1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HR2   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HR3   TYPE=HA   CHARGE=     .0900  END  
 !END GROUP
 
      BOND  CL    CLP    
      BOND  CLP   NL     
      BOND  NL    CA     
      BOND  CA    CRP    
      BOND  CRP   NR     
      BOND  NR    CR     
      BOND  CLP   OL     
      BOND  CRP   OR     
      BOND  NL    HL     
      BOND  NR    HR     
      BOND  CA    HA     
      BOND  CA    CB     
      BOND  CL    HL1    
      BOND  CL    HL2    
      BOND  CL    HL3    
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CB    HB3    
      BOND  CR    HR1    
      BOND  CR    HR2    
      BOND  CR    HR3    

      DIHEDRAL  CL    CLP   NL    CA    ! multiple dihedral 
      DIHEDRAL  CL    CLP   NL    CA    ! multiple dihedral 
      DIHEDRAL  CA    CRP   NR    CR    ! multiple dihedral 
      DIHEDRAL  CA    CRP   NR    CR    ! multiple dihedral 

      IMPROPER  CLP   CL    NL    OL     
      IMPROPER  NL    CLP   CA    HL     
      IMPROPER  CRP   CA    NR    OR     
      IMPROPER  NR    CRP   CR    HR     
 
      IC  CLP   NL     CA    CRP       .0000     .00  180.00     .00    .0000 ! Phi 
      IC  CA    CLP    *NL   HL        .0000     .00  180.00     .00    .0000  
      IC  HL    NL     CA    CRP       .0000     .00     .00     .00    .0000  
      IC  NL    CA     CRP   NR        .0000     .00  180.00     .00    .0000 ! Psi 
      IC  CA    NR     *CRP  OR        .0000     .00  180.00     .00    .0000  
      IC  NL    CA     CRP   OR        .0000     .00     .00     .00    .0000  
      IC  CL    CLP    NL    CA        .0000     .00  180.00     .00    .0000 ! Omega Left 
      IC  NL    CL     *CLP  OL        .0000     .00  180.00     .00    .0000  
      IC  OL    CLP    NL    CA        .0000     .00     .00     .00    .0000  
      IC  CA    CRP    NR    CR        .0000     .00  180.00     .00    .0000 ! Omega Right 
      IC  CRP   CR     *NR   HR        .0000     .00  180.00     .00    .0000  
      IC  CA    CRP    NR    HR        .0000     .00  180.00     .00    .0000  
      IC  NL    CRP    *CA   HA        .0000     .00  240.00     .00    .0000  
      IC  NL    CRP    *CA   CB        .0000     .00  120.00     .00    .0000  
      IC  HL1   CL     CLP   NL        .0000     .00  180.00     .00    .0000  
      IC  HL2   CL     CLP   NL        .0000     .00   60.00     .00    .0000  
      IC  HL3   CL     CLP   OL        .0000     .00  120.00     .00    .0000  
      IC  HA    CA     CB    HB1       .0000     .00  180.00     .00    .0000  
      IC  NL    CA     CB    HB2       .0000     .00  180.00     .00    .0000  
      IC  CRP   CA     CB    HB3       .0000     .00  180.00     .00    .0000  
      IC  CRP   NR     CR    HR1       .0000     .00  180.00     .00    .0000  
      IC  CRP   NR     CR    HR2       .0000     .00   60.00     .00    .0000  
      IC  HR    NR     CR    HR3       .0000     .00  120.00     .00    .0000  
      IC  CA    CLP    *NL   HL        .0000     .00  180.00     .00    .0000  
      IC  CA    NR     *CRP  OR        .0000     .00  180.00     .00    .0000  
      IC  HB1   HB2    *CB   HB3       .0000     .00  120.00     .00    .0000  
      IC  HL1   HL2    *CL   HL3       .0000     .00  240.00     .00    .0000  
      IC  HR1   HR2    *CR   HR3       .0000     .00  240.00     .00    .0000  
      IC  HA    CA     NL    HL        .0000     .00  240.00     .00    .0000  
 
 
 END {ALAD}
 !-----------------------------------------------------------
 
 RESIDUE MESH ! methanethiol, DZUNG NGUYEN 
 
 GROUP  
       ATOM  H1    TYPE=HA   CHARGE=     .0900  END !   H1 
       ATOM  H2    TYPE=HA   CHARGE=     .0900  END !     \ 
       ATOM  H3    TYPE=HA   CHARGE=     .0900  END ! H2--CM--S 
       ATOM  CM    TYPE=CT3  CHARGE=    -.2000  END !     /    \ 
       ATOM  S     TYPE=S    CHARGE=    -.2300  END !   H3     H4 
       ATOM  H4    TYPE=HS   CHARGE=     .1600  END ! 
 !END GROUP
 
      BOND  CM    H1     
      BOND  CM    H2     
      BOND  CM    H3     
      BOND  CM    S      
      BOND  S     H4     
 
      IC  H1    CM     S     H4        .0000     .00     .00     .00    .0000  
      IC  H2    CM     S     H4        .0000     .00  120.00     .00    .0000  
      IC  H3    CM     S     H4        .0000     .00  240.00     .00    .0000  
      IC  CM    S      H4    H1        .0000     .00     .00     .00    .0000  
 
 
 END {MESH}
 !-----------------------------------------------------------
 
 RESIDUE MES1 ! methylthiolate, adm jr. 
 
 GROUP  
       ATOM  S     TYPE=SS   CHARGE=    -.8000  END !  H11 
       ATOM  C1    TYPE=CS   CHARGE=    -.4700  END !     \ 
       ATOM  H11   TYPE=HA   CHARGE=     .0900  END ! H12--C1--S 
       ATOM  H12   TYPE=HA   CHARGE=     .0900  END !     / 
       ATOM  H13   TYPE=HA   CHARGE=     .0900  END !  H13 
      ! 
 !END GROUP
 
      BOND  S     C1     
      BOND  C1    H11    
      BOND  C1    H12    
      BOND  C1    H13    
 
      IC  BLNK  H11    C1    S         .0000     .00   60.00     .00    .0000  
      IC  S     H11    *C1   H12       .0000     .00  120.00     .00    .0000  
      IC  S     H11    *C1   H13       .0000     .00 -120.00     .00    .0000  
      IC  H12   H13    *C1   H11       .0000     .00  120.00     .00    .0000  
 
 
 END {MES1}
 !-----------------------------------------------------------
 
 RESIDUE ETSH ! ethanethiol, Dzung Nguyen 
 
 GROUP  
       ATOM  H1    TYPE=HA   CHARGE=     .0900  END !   H1  H4   H5 
       ATOM  H2    TYPE=HA   CHARGE=     .0900  END !     \   \  / 
       ATOM  H3    TYPE=HA   CHARGE=     .0900  END !  H2-CM1--CM2 
       ATOM  CM1   TYPE=CT3  CHARGE=    -.2700  END !     /      \ 
 GROUP						    !   H3        S3--H6 
       ATOM  H4    TYPE=HA   CHARGE=     .0900  END  
       ATOM  H5    TYPE=HA   CHARGE=     .0900  END  
       ATOM  CM2   TYPE=CT2  CHARGE=    -.1100  END  
       ATOM  S3    TYPE=S    CHARGE=    -.2300  END  
       ATOM  H6    TYPE=HS   CHARGE=     .1600  END  
 !END GROUP
 
      BOND  CM1   H1     
      BOND  CM1   H2     
      BOND  CM1   H3     
      BOND  CM1   CM2    
      BOND  CM2   H4     
      BOND  CM2   H5     
      BOND  CM2   S3     
      BOND  S3    H6     

      DIHEDRAL  CM1  CM2  S3  H6 ! multiple dihedral
      DIHEDRAL  CM1  CM2  S3  H6 ! multiple dihedral
      DIHEDRAL  CM1  CM2  S3  H6 ! multiple dihedral
 
      IC  H1    CM1    CM2   S3        .0000     .00   60.00     .00    .0000  
      IC  H2    CM1    CM2   S3        .0000     .00  180.00     .00    .0000  
      IC  H3    CM1    CM2   S3        .0000     .00  300.00     .00    .0000  
      IC  H1    CM1    CM2   H4        .0000     .00  180.00     .00    .0000  
      IC  H2    CM1    CM2   H4        .0000     .00  300.00     .00    .0000  
      IC  H3    CM1    CM2   H4        .0000     .00   60.00     .00    .0000  
      IC  H1    CM1    CM2   H5        .0000     .00  300.00     .00    .0000  
      IC  H2    CM1    CM2   H5        .0000     .00   60.00     .00    .0000  
      IC  H3    CM1    CM2   H5        .0000     .00  180.00     .00    .0000  
      IC  CM1   CM2    S3    H6        .0000     .00   60.00     .00    .0000  
      IC  H4    CM2    S3    H6        .0000     .00  180.00     .00    .0000  
      IC  H5    CM2    S3    H6        .0000     .00  300.00     .00    .0000  
 
 
 END {ETSH}
 !-----------------------------------------------------------
 
 RESIDUE ES1  ! ethylthiolate, adm jr. 
 
 GROUP  
       ATOM  S     TYPE=SS   CHARGE=    -.8000  END !  H21  H11 H12 
       ATOM  C1    TYPE=CS   CHARGE=    -.3800  END !     \   \  / 
       ATOM  C2    TYPE=CT3  CHARGE=    -.2700  END ! H22--C2--C1 
       ATOM  H21   TYPE=HA   CHARGE=     .0900  END !     /      \ 
       ATOM  H22   TYPE=HA   CHARGE=     .0900  END !  H23        S 
       ATOM  H23   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H11   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H12   TYPE=HA   CHARGE=     .0900  END  
 !END GROUP
 
      BOND  S     C1     
      BOND  C1    H11    
      BOND  C1    H12    
      BOND  C1    C2     
      BOND  C2    H21    
      BOND  C2    H22    
      BOND  C2    H23    
 
      IC  BLNK  H11    C1    S         .0000     .00   60.00     .00    .0000  
      IC  S     H11    *C1   H12       .0000     .00  120.00     .00    .0000  
      IC  S     H11    *C1   C2        .0000     .00 -120.00     .00    .0000  
      IC  H12   C2     *C1   H11       .0000     .00  120.00     .00    .0000  
      IC  S     C1     C2    H21       .0000     .00  180.00     .00    .0000  
      IC  S     C1     C2    H22       .0000     .00   60.00     .00    .0000  
      IC  S     C1     C2    H23       .0000     .00  -60.00     .00    .0000  
 
 
 END {ES1 }
 !-----------------------------------------------------------
 
 RESIDUE DMDS ! dimethyldisulfide, Dzung Nguyen 
 
 GROUP  
       ATOM  H1    TYPE=HA   CHARGE=     .0900  END ! 
       ATOM  H2    TYPE=HA   CHARGE=     .0900  END !   H1 
       ATOM  H3    TYPE=HA   CHARGE=     .0900  END !    \ 
       ATOM  CM1   TYPE=CT3  CHARGE=    -.1900  END ! H2-CM1 
       ATOM  S2    TYPE=SM   CHARGE=    -.0800  END !    /  \ 
 GROUP						    !   H3   S2--S3    H4 
       ATOM  S3    TYPE=SM   CHARGE=    -.0800  END !              \  / 
       ATOM  CM4   TYPE=CT3  CHARGE=    -.1900  END !               CM4-H5 
       ATOM  H4    TYPE=HA   CHARGE=     .0900  END !                 \ 
       ATOM  H5    TYPE=HA   CHARGE=     .0900  END !                  H6 
       ATOM  H6    TYPE=HA   CHARGE=     .0900  END ! 
 !END GROUP
 
      BOND  H1    CM1    
      BOND  H2    CM1    
      BOND  H3    CM1    
      BOND  CM1   S2     
      BOND  S2    S3     
      BOND  S3    CM4    
      BOND  CM4   H4     
      BOND  CM4   H5     
      BOND  CM4   H6     
      
      DIHEDRAL CM1 S2 S3 CM4  ! required for multiple dihedrals
      DIHEDRAL CM1 S2 S3 CM4  ! required for multiple dihedrals
      DIHEDRAL CM1 S2 S3 CM4  ! required for multiple dihedrals

      IC  H1    CM1    S2    S3        .0000     .00   60.00     .00    .0000  
      IC  H2    CM1    S2    S3        .0000     .00  180.00     .00    .0000  
      IC  H3    CM1    S2    S3        .0000     .00  300.00     .00    .0000  
      IC  CM1   S2     S3    CM4       .0000     .00   90.00     .00    .0000  
      IC  S2    S3     CM4   H4        .0000     .00   60.00     .00    .0000  
      IC  S2    S3     CM4   H5        .0000     .00  180.00     .00    .0000  
      IC  S2    S3     CM4   H6        .0000     .00  300.00     .00    .0000  
 
 
 END {DMDS}
 !-----------------------------------------------------------
 
 RESIDUE EMS  ! ethylmethylsulfide, Dzung Nguyen 
 
 GROUP  
       ATOM  H1    TYPE=HA   CHARGE=     .0900  END ! 
       ATOM  H2    TYPE=HA   CHARGE=     .0900  END !   H1   H4 H5 
       ATOM  H3    TYPE=HA   CHARGE=     .0900  END !    \    \ / 
       ATOM  CM1   TYPE=CT3  CHARGE=    -.2700  END ! H2-CM1--CM2        H6 
 GROUP						    !            \      / 
       ATOM  H4    TYPE=HA   CHARGE=     .0900  END !             S3--CM4-H7 
       ATOM  H5    TYPE=HA   CHARGE=     .0900  END !                   \ 
       ATOM  CM2   TYPE=CT2  CHARGE=    -.1400  END !                    H8 
       ATOM  S3    TYPE=S    CHARGE=    -.0900  END ! 
       ATOM  CM4   TYPE=CT3  CHARGE=    -.2200  END ! 
       ATOM  H6    TYPE=HA   CHARGE=     .0900  END ! 
       ATOM  H7    TYPE=HA   CHARGE=     .0900  END ! 
       ATOM  H8    TYPE=HA   CHARGE=     .0900  END ! 
 !END GROUP
 
      BOND  CM1   H1     
      BOND  CM1   H2     
      BOND  CM1   H3     
      BOND  CM1   CM2    
      BOND  CM2   H4     
      BOND  CM2   H5     
      BOND  CM2   S3     
      BOND  S3    CM4    
      BOND  CM4   H6     
      BOND  CM4   H7     
      BOND  CM4   H8     

      DIHEDRAL  CM1  CM2  S3   CM4  ! multiple dihedral 
      DIHEDRAL  CM1  CM2  S3   CM4  ! multiple dihedral 

      IC  H1    CM1    CM2   S3        .0000     .00   60.00     .00    .0000  
      IC  H2    CM1    CM2   S3        .0000     .00  180.00     .00    .0000  
      IC  H3    CM1    CM2   S3        .0000     .00  300.00     .00    .0000  
      IC  H1    CM1    CM2   H4        .0000     .00  180.00     .00    .0000  
      IC  H2    CM1    CM2   H4        .0000     .00  300.00     .00    .0000  
      IC  H3    CM1    CM2   H4        .0000     .00   60.00     .00    .0000  
      IC  H1    CM1    CM2   H5        .0000     .00  300.00     .00    .0000  
      IC  H2    CM1    CM2   H5        .0000     .00   60.00     .00    .0000  
      IC  H3    CM1    CM2   H5        .0000     .00  180.00     .00    .0000  
      IC  CM1   CM2    S3    CM4       .0000     .00   60.00     .00    .0000  
      IC  H4    CM2    S3    CM4       .0000     .00  180.00     .00    .0000  
      IC  H5    CM2    S3    CM4       .0000     .00  300.00     .00    .0000  
      IC  CM2   S3     CM4   H6        .0000     .00   60.00     .00    .0000  
      IC  CM2   S3     CM4   H7        .0000     .00  180.00     .00    .0000  
      IC  CM2   S3     CM4   H8        .0000     .00  300.00     .00    .0000  
 
 
 END {EMS }
 !-----------------------------------------------------------
 
 RESIDUE IMIA ! Imidazole, adm jr. 
 
 GROUP  
       ATOM  CG    TYPE=CPH1 CHARGE=    -.0500  END !       HD1    HE1 
       ATOM  HG    TYPE=HR3  CHARGE=     .0900  END !        |     / 
       ATOM  CD2   TYPE=CPH1 CHARGE=     .2200  END !       ND1--CE1 
       ATOM  HD2   TYPE=HR3  CHARGE=     .1000  END !      /      | 
       ATOM  ND1   TYPE=NR1  CHARGE=    -.3600  END ! HG-CG       | 
       ATOM  HD1   TYPE=H    CHARGE=     .3200  END !      \      | 
       ATOM  CE1   TYPE=CPH2 CHARGE=     .2500  END !      CD2--NE2 
       ATOM  HE1   TYPE=HR1  CHARGE=     .1300  END !       | 
       ATOM  NE2   TYPE=NR2  CHARGE=    -.7000  END !      HD2 
      ! 
 !END GROUP
 
      BOND  NE2   CD2    
      BOND  ND1   CG     
      BOND  CD2   CG     
      BOND  CE1   ND1    
      BOND  NE2   CE1    
      BOND  CG    HG     
      BOND  ND1   HD1    
      BOND  CD2   HD2    
      BOND  CE1   HE1

      ! KEEPS HYDROGENS IN RING PLANE 
 
      IMPROPER  ND1   CG    CE1   HD1    
      IMPROPER  ND1   CE1   CG    HD1    
      IMPROPER  CD2   CG    NE2   HD2    
      IMPROPER  CD2   NE2   CG    HD2    
      IMPROPER  CE1   ND1   NE2   HE1    
      IMPROPER  CE1   NE2   ND1   HE1    
      IMPROPER  CG    CD2   ND1   HG     
      IMPROPER  CG    ND1   CD2   HG     
 
      DONOR  HD1   ND1    
      DONOR  HE1   CE1    
      DONOR  HG    CG     
 
      ACCEPTOR  NE2   NONE   
      ACCEPTOR  CE1   NONE   
      ACCEPTOR  ND1   NONE   
      ACCEPTOR  CD2   NONE   
 
      IC  HG    CG     ND1   CE1      1.5421  122.67 -173.67  109.79   1.2987  
      IC  CG    ND1    CE1   NE2      1.2854  109.79     .21  110.31   1.3071  
      IC  ND1   CE1    NE2   CD2      1.2987  110.31     .03  105.82   1.3165  
      IC  CE1   NE2    CD2   CG       1.3071  105.82    -.23  108.68   1.3758  
      IC  NE2   CD2    CG    ND1      1.3165  108.68     .35  105.39   1.2854  
      IC  NE2   CD2    CG    HG       1.3165  108.68  172.86  131.52   1.5421  
      IC  CD2   CG     ND1   CE1      1.3758  105.39    -.34  109.79   1.2987  
      IC  CD2   NE2    CE1   HE1      1.3165  105.82  149.51  119.57   1.0879  
      IC  NE2   CE1    ND1   HD1      1.3071  110.31  157.04  123.39    .9770  
      IC  HG    CG     CD2   HD2      1.5421  131.52  -48.16  118.30   1.0902  
 
 
 END {IMIA}
 !-----------------------------------------------------------
 
 RESIDUE MIMI ! 4-methylimidazole, adm jr. 
 
 GROUP  
       ATOM  ND1   TYPE=NR1  CHARGE=    -.3600  END !           HD1    HE1 
       ATOM  HD1   TYPE=H    CHARGE=     .3200  END !            |     / 
       ATOM  CG    TYPE=CPH1 CHARGE=    -.0500  END !    HB1    ND1--CE1 
       ATOM  CB    TYPE=CT2  CHARGE=    -.1800  END !     |    /      | 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END ! HB2-CB-CG       | 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !     |    \      | 
       ATOM  HB3   TYPE=HA   CHARGE=     .0900  END !    HB3    CD2--NE2 
 GROUP						    !            | 
       ATOM  NE2   TYPE=NR2  CHARGE=    -.7000  END !           HD2 
       ATOM  CD2   TYPE=CPH1 CHARGE=     .2200  END ! 
       ATOM  HD2   TYPE=HR3  CHARGE=     .1000  END ! 
       ATOM  CE1   TYPE=CPH2 CHARGE=     .2500  END ! 
       ATOM  HE1   TYPE=HR1  CHARGE=     .1300  END ! 
 !END GROUP
 
      BOND  NE2   CD2    
      BOND  ND1   CG     
      BOND  CD2   CG     
      BOND  CE1   ND1    
      BOND  NE2   CE1    
      BOND  CG    CB     
      BOND  ND1   HD1    
      BOND  CD2   HD2    
      BOND  CE1   HE1    
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CB    HB3

      ! KEEPS HYDROGENS IN RING PLANE 
 
      IMPROPER  ND1   CG    CE1   HD1    
      IMPROPER  CD2   CG    NE2   HD2    
      IMPROPER  CE1   ND1   NE2   HE1    
 
      DONOR  HD1   ND1    
 
      ACCEPTOR  NE2   NONE   
 
      IC  HG    CG     ND1   CE1      1.5421  122.67 -173.67  109.79   1.2987  
      IC  CG    ND1    CE1   NE2      1.2854  109.79     .21  110.31   1.3071  
      IC  ND1   CE1    NE2   CD2      1.2987  110.31     .03  105.82   1.3165  
      IC  CE1   NE2    CD2   CG       1.3071  105.82    -.23  108.68   1.3758  
      IC  NE2   CD2    CG    ND1      1.3165  108.68     .35  105.39   1.2854  
      IC  NE2   CD2    CG    HG       1.3165  108.68  172.86  131.52   1.5421  
      IC  CD2   CG     ND1   CE1      1.3758  105.39    -.34  109.79   1.2987  
      IC  CD2   NE2    CE1   HE1      1.3165  105.82  149.51  119.57   1.0879  
      IC  NE2   CE1    ND1   HD1      1.3071  110.31  157.04  123.39    .9770  
      IC  HG    CG     CD2   HD2      1.5421  131.52  -48.16  118.30   1.0902  
      IC  CD2   CG     CB    HB1       .0000     .00     .00     .00    .0000  
      IC  CD2   CG     CB    HB2       .0000     .00  120.00     .00    .0000  
      IC  CD2   CG     CB    HB3       .0000     .00  240.00     .00    .0000  
 
 
 END {MIMI}
 !-----------------------------------------------------------
 
 RESIDUE EIMI ! 4-ethylimidazole, adm jr. 
 
 GROUP  
       ATOM  CA    TYPE=CT3  CHARGE=    -.2700  END ! Optimized charges to minimize group size 
       ATOM  HA1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HA2   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HA3   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  ND1   TYPE=NR1  CHARGE=    -.3600  END !               HD1    HE1 
       ATOM  HD1   TYPE=H    CHARGE=     .3200  END !                |     / 
       ATOM  CG    TYPE=CPH1 CHARGE=    -.0500  END !  HA1   HB1    ND1--CE1 
       ATOM  CB    TYPE=CT2  CHARGE=    -.0900  END !    \    |    /      | 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END ! HA2-CA--CB-CG       | 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !    /    |    \      | 
 GROUP						    !  HA3   HB2    CD2--NE2 
       ATOM  NE2   TYPE=NR2  CHARGE=    -.7000  END !                | 
       ATOM  CD2   TYPE=CPH1 CHARGE=     .2200  END !               HD2 
       ATOM  HD2   TYPE=HR3  CHARGE=     .1000  END ! 
       ATOM  CE1   TYPE=CPH2 CHARGE=     .2500  END ! 
       ATOM  HE1   TYPE=HR1  CHARGE=     .1300  END ! 
 !END GROUP
 
      BOND  NE2   CD2    
      BOND  ND1   CG     
      BOND  CD2   CG     
      BOND  CE1   ND1    
      BOND  NE2   CE1    
      BOND  CG    CB     
      BOND  ND1   HD1    
      BOND  CD2   HD2    
      BOND  CE1   HE1    
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CB    CA     
      BOND  CA    HA1    
      BOND  CA    HA2    
      BOND  CA    HA3

      ! KEEPS HYDROGENS IN RING PLANE 
 
      IMPROPER  ND1   CG    CE1   HD1    
      IMPROPER  CD2   CG    NE2   HD2    
      IMPROPER  CE1   ND1   NE2   HE1    
 
      DONOR  HD1   ND1    
 
      ACCEPTOR  NE2   NONE   
 
      IC  HG    CG     ND1   CE1      1.5421  122.67 -173.67  109.79   1.2987  
      IC  CG    ND1    CE1   NE2      1.2854  109.79     .21  110.31   1.3071  
      IC  ND1   CE1    NE2   CD2      1.2987  110.31     .03  105.82   1.3165  
      IC  CE1   NE2    CD2   CG       1.3071  105.82    -.23  108.68   1.3758  
      IC  NE2   CD2    CG    ND1      1.3165  108.68     .35  105.39   1.2854  
      IC  NE2   CD2    CG    HG       1.3165  108.68  172.86  131.52   1.5421  
      IC  CD2   CG     ND1   CE1      1.3758  105.39    -.34  109.79   1.2987  
      IC  CD2   NE2    CE1   HE1      1.3165  105.82  149.51  119.57   1.0879  
      IC  NE2   CE1    ND1   HD1      1.3071  110.31  157.04  123.39    .9770  
      IC  HG    CG     CD2   HD2      1.5421  131.52  -48.16  118.30   1.0902  
      IC  CD2   CG     CB    HB1       .0000     .00     .00     .00    .0000  
      IC  CD2   CG     CB    HB2       .0000     .00  120.00     .00    .0000  
      IC  CD2   CG     CB    CA        .0000     .00  240.00     .00    .0000  
      IC  CG    CB     CA    HA1       .0000     .00 -180.00     .00    .0000  
      IC  CG    CB     CA    HA2       .0000     .00  -60.00     .00    .0000  
      IC  CG    CB     CA    HA3       .0000     .00   60.00     .00    .0000  
 
 
 END {EIMI}
 !-----------------------------------------------------------
 
 RESIDUE IMIM ! Imidazolium, adm jr. 
 
 GROUP  
       ATOM  CG    TYPE=CPH1 CHARGE=     .1900  END !       HD1    HE1 
       ATOM  HG    TYPE=HR1  CHARGE=     .1300  END !        |     / 
       ATOM  CD2   TYPE=CPH1 CHARGE=     .1900  END !       ND1--CE1 
       ATOM  HD2   TYPE=HR1  CHARGE=     .1300  END !      /      | 
 GROUP						    ! HG-CG       | 
       ATOM  ND1   TYPE=NR3  CHARGE=    -.5100  END !      \      | 
       ATOM  HD1   TYPE=H    CHARGE=     .4400  END !      CD2--NE2 
       ATOM  NE2   TYPE=NR3  CHARGE=    -.5100  END !       |      \ 
       ATOM  HE2   TYPE=H    CHARGE=     .4400  END !      HD2     HE2 
       ATOM  CE1   TYPE=CPH2 CHARGE=     .3200  END  
       ATOM  HE1   TYPE=HR2  CHARGE=     .1800  END  
 !END GROUP
 
      BOND  NE2   CD2    
      BOND  ND1   CG     
      BOND  CD2   CG     
      BOND  CE1   ND1    
      BOND  NE2   CE1    
      BOND  CG    HG     
      BOND  ND1   HD1    
      BOND  CD2   HD2    
      BOND  CE1   HE1    
      BOND  NE2   HE2

      ! KEEPS HYDROGENS IN RING PLANE 
 
      IMPROPER  HD1   CG    CE1   ND1    
      IMPROPER  HE2   CD2   CE1   NE2    
      IMPROPER  HD1   CE1   CG    ND1    
      IMPROPER  HE2   CE1   CD2   NE2    
 
      DONOR  HD1   ND1    
      DONOR  HE2   NE2    
 
      IC  HG    CG     ND1   CE1      1.5421  122.67 -173.67  109.79   1.2987  
      IC  CG    ND1    CE1   NE2      1.2854  109.79     .21  110.31   1.3071  
      IC  ND1   CE1    NE2   CD2      1.2987  110.31     .03  105.82   1.3165  
      IC  CE1   NE2    CD2   CG       1.3071  105.82    -.23  108.68   1.3758  
      IC  NE2   CD2    CG    ND1      1.3165  108.68     .35  105.39   1.2854  
      IC  NE2   CD2    CG    HG       1.3165  108.68  172.86  131.52   1.5421  
      IC  CD2   CG     ND1   CE1      1.3758  105.39    -.34  109.79   1.2987  
      IC  CD2   NE2    CE1   HE1      1.3165  105.82  149.51  119.57   1.0879  
      IC  NE2   CE1    ND1   HD1      1.3071  110.31  157.04  123.39    .9770  
      IC  HG    CG     CD2   HD2      1.5421  131.52  -48.16  118.30   1.0902  
      IC  HE1   CE1    NE2   HE2      1.0879  125.00     .00  125.00   1.0000  
 
 
 END {IMIM}
 !-----------------------------------------------------------
 
 RESIDUE EIMM ! Ethyl-Imidazolium, adm jr. 
 
 GROUP  
       ATOM  CA    TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  HA1   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HA2   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HA3   TYPE=HA   CHARGE=     .0900  END  
 GROUP ! 
       ATOM  CB    TYPE=CT2  CHARGE=    -.0500  END !               HD1    HE1 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !                |     / 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !  HA1   HB1    ND1--CE1 
       ATOM  CG    TYPE=CPH1 CHARGE=     .1900  END !    \    |    /      | 
       ATOM  CD2   TYPE=CPH1 CHARGE=     .1900  END ! HA2-CA--CB-CG       | 
       ATOM  HD2   TYPE=HR1  CHARGE=     .1300  END !    /    |    \      | 
 GROUP						    !  HA3   HB2    CD2--NE2 
       ATOM  ND1   TYPE=NR3  CHARGE=    -.5100  END !                |     \ 
       ATOM  HD1   TYPE=H    CHARGE=     .4400  END !               HD2    HE2 
       ATOM  NE2   TYPE=NR3  CHARGE=    -.5100  END ! 
       ATOM  HE2   TYPE=H    CHARGE=     .4400  END  
       ATOM  CE1   TYPE=CPH2 CHARGE=     .3200  END  
       ATOM  HE1   TYPE=HR2  CHARGE=     .1800  END  
 !END GROUP
 
      BOND  NE2   CD2    
      BOND  ND1   CG     
      BOND  CD2   CG     
      BOND  CE1   ND1    
      BOND  NE2   CE1    
      BOND  CG    CB     
      BOND  NE2   HE2    
      BOND  ND1   HD1    
      BOND  CD2   HD2    
      BOND  CE1   HE1    
      BOND  CB    HB1    
      BOND  CB    HB2    
      BOND  CB    CA     
      BOND  CA    HA1    
      BOND  CA    HA2    
      BOND  CA    HA3

      ! KEEPS HYDROGENS IN RING PLANE 
 
      IMPROPER  HD1   CG    CE1   ND1    
      IMPROPER  HE2   CD2   CE1   NE2    
      IMPROPER  HD1   CE1   CG    ND1    
      IMPROPER  HE2   CE1   CD2   NE2    
 
      DONOR  HD1   ND1    
 
      ACCEPTOR  NE2   NONE   
 
      IC  HG    CG     ND1   CE1      1.5421  122.67 -173.67  109.79   1.2987  
      IC  CG    ND1    CE1   NE2      1.2854  109.79     .21  110.31   1.3071  
      IC  ND1   CE1    NE2   CD2      1.2987  110.31     .03  105.82   1.3165  
      IC  CE1   NE2    CD2   CG       1.3071  105.82    -.23  108.68   1.3758  
      IC  NE2   CD2    CG    ND1      1.3165  108.68     .35  105.39   1.2854  
      IC  NE2   CD2    CG    HG       1.3165  108.68  172.86  131.52   1.5421  
      IC  CD2   CG     ND1   CE1      1.3758  105.39    -.34  109.79   1.2987  
      IC  CD2   NE2    CE1   HE1      1.3165  105.82  149.51  119.57   1.0879  
      IC  NE2   CE1    ND1   HD1      1.3071  110.31  157.04  123.39    .9770  
      IC  HG    CG     CD2   HD2      1.5421  131.52  -48.16  118.30   1.0902  
      IC  CD2   CG     CB    HB1       .0000     .00     .00     .00    .0000  
      IC  CD2   CG     CB    HB2       .0000     .00  120.00     .00    .0000  
      IC  CD2   CG     CB    CA        .0000     .00  240.00     .00    .0000  
      IC  CG    CB     CA    HA1       .0000     .00 -180.00     .00    .0000  
      IC  CG    CB     CA    HA2       .0000     .00  -60.00     .00    .0000  
      IC  CG    CB     CA    HA3       .0000     .00   60.00     .00    .0000  
 
 
 END {EIMM}
 !-----------------------------------------------------------
 
 RESIDUE BENZ ! benzene, adm jr. 
 
 GROUP  
       ATOM  CG    TYPE=CA   CHARGE=    -.1150  END ! 
       ATOM  HG    TYPE=HP   CHARGE=     .1150  END !      HD1  HE1 
 GROUP						    !       |    | 
       ATOM  CD1   TYPE=CA   CHARGE=    -.1150  END !      CD1--CE1 
       ATOM  HD1   TYPE=HP   CHARGE=     .1150  END !      /      \ 
 GROUP						    ! HG--CG      CZ--HZ 
       ATOM  CD2   TYPE=CA   CHARGE=    -.1150  END !      \      / 
       ATOM  HD2   TYPE=HP   CHARGE=     .1150  END !      CD2--CE2 
 GROUP						    !       |    | 
       ATOM  CE1   TYPE=CA   CHARGE=    -.1150  END !      HD2  HE2 
       ATOM  HE1   TYPE=HP   CHARGE=     .1150  END ! 
 GROUP  
       ATOM  CE2   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HE2   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CZ    TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HZ    TYPE=HP   CHARGE=     .1150  END  
 !END GROUP
 
      BOND  CD1   CG     
      BOND  CD2   CG     
      BOND  CE1   CD1    
      BOND  CE2   CD2    
      BOND  CZ    CE1    
      BOND  CZ    CE2    
      BOND  CG    HG     
      BOND  CD1   HD1    
      BOND  CD2   HD2    
      BOND  CE1   HE1    
      BOND  CE2   HE2    
      BOND  CZ    HZ     
 
      IC  CG    CD1    CE1   CZ        .0000     .00     .00     .00    .0000  
      IC  CD1   CE1    CZ    CE2       .0000     .00     .00     .00    .0000  
      IC  CE1   CZ     CE2   CD2       .0000     .00     .00     .00    .0000  
      IC  CD1   CD2    *CG   HG        .0000     .00  180.00     .00    .0000  
      IC  CE1   CG     *CD1  HD1       .0000     .00  180.00     .00    .0000  
      IC  CE2   CG     *CD2  HD2       .0000     .00  180.00     .00    .0000  
      IC  CZ    CD1    *CE1  HE1       .0000     .00  180.00     .00    .0000  
      IC  CZ    CD2    *CE2  HE2       .0000     .00  180.00     .00    .0000  
      IC  CE1   CE2    *CZ   HZ        .0000     .00  180.00     .00    .0000  
 
 
 END {BENZ}
 !-----------------------------------------------------------
 
 RESIDUE PHEN ! phenol, adm jr. 
 
 GROUP  
       ATOM  CG    TYPE=CA   CHARGE=    -.1150  END ! 
       ATOM  HG    TYPE=HP   CHARGE=     .1150  END !      HD1  HE1 
 GROUP						    !       |    | 
       ATOM  CD1   TYPE=CA   CHARGE=    -.1150  END !      CD1--CE1 
       ATOM  HD1   TYPE=HP   CHARGE=     .1150  END !      /      \ 
 GROUP						    ! HG--CG      CZ--OH 
       ATOM  CD2   TYPE=CA   CHARGE=    -.1150  END !      \      /     \ 
       ATOM  HD2   TYPE=HP   CHARGE=     .1150  END !      CD2--CE2      HH 
 GROUP						    !       |    | 
       ATOM  CE1   TYPE=CA   CHARGE=    -.1150  END !      HD2  HE2 
       ATOM  HE1   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CE2   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HE2   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CZ    TYPE=CA   CHARGE=     .1100  END  
       ATOM  OH    TYPE=OH1  CHARGE=    -.5400  END  
       ATOM  HH    TYPE=H    CHARGE=     .4300  END  
 !END GROUP
 
      BOND  CD1   CG     
      BOND  CD2   CG     
      BOND  CE1   CD1    
      BOND  CE2   CD2    
      BOND  CZ    CE1    
      BOND  CZ    CE2    
      BOND  CG    HG     
      BOND  CD1   HD1    
      BOND  CD2   HD2    
      BOND  CE1   HE1    
      BOND  CE2   HE2    
      BOND  CZ    OH     
      BOND  OH    HH     
 
      IC  CG    CD1    CE1   CZ        .0000     .00     .00     .00    .0000  
      IC  CD1   CE1    CZ    CE2       .0000     .00     .00     .00    .0000  
      IC  CE1   CZ     CE2   CD2       .0000     .00     .00     .00    .0000  
      IC  CD1   CD2    *CG   HG        .0000     .00  180.00     .00    .0000  
      IC  CE1   CG     *CD1  HD1       .0000     .00  180.00     .00    .0000  
      IC  CE2   CG     *CD2  HD2       .0000     .00  180.00     .00    .0000  
      IC  CZ    CD1    *CE1  HE1       .0000     .00  180.00     .00    .0000  
      IC  CZ    CD2    *CE2  HE2       .0000     .00  180.00     .00    .0000  
      IC  CE1   CE2    *CZ   OH        .0000     .00  180.00     .00    .0000  
      IC  CE1   CZ     OH    HH        .0000     .00  180.00     .00    .0000  
 
 
 END {PHEN}
 !-----------------------------------------------------------
 
 RESIDUE PHEO ! phenoxide, adm jr. 
 
 GROUP  
       ATOM  CG    TYPE=CA   CHARGE=    -.1150  END ! 
       ATOM  HG    TYPE=HP   CHARGE=     .1150  END !      HD1  HE1 
 GROUP						    !       |    | 
       ATOM  CD1   TYPE=CA   CHARGE=    -.1150  END !      CD1--CE1 
       ATOM  HD1   TYPE=HP   CHARGE=     .1150  END !      /      \ 
 GROUP						    ! HG--CG      CZ--OH 
       ATOM  CD2   TYPE=CA   CHARGE=    -.1150  END !      \      / 
       ATOM  HD2   TYPE=HP   CHARGE=     .1150  END !      CD2--CE2 
 GROUP						    !       |    | 
       ATOM  CE1   TYPE=CA   CHARGE=    -.6000  END !      HD2  HE2 
       ATOM  HE1   TYPE=HP   CHARGE=     .2800  END  
       ATOM  CE2   TYPE=CA   CHARGE=    -.6000  END  
       ATOM  HE2   TYPE=HP   CHARGE=     .2800  END  
       ATOM  CZ    TYPE=CA   CHARGE=     .4000  END  
       ATOM  OH    TYPE=OC   CHARGE=    -.7600  END  
 !END GROUP
 
      BOND  CD1   CG     
      BOND  CD2   CG     
      BOND  CE1   CD1    
      BOND  CE2   CD2    
      BOND  CZ    CE1    
      BOND  CZ    CE2    
      BOND  CG    HG     
      BOND  CD1   HD1    
      BOND  CD2   HD2    
      BOND  CE1   HE1    
      BOND  CE2   HE2    
      BOND  CZ    OH     
 
      IC  CG    CD1    CE1   CZ        .0000     .00     .00     .00    .0000  
      IC  CD1   CE1    CZ    CE2       .0000     .00     .00     .00    .0000  
      IC  CE1   CZ     CE2   CD2       .0000     .00     .00     .00    .0000  
      IC  CD1   CD2    *CG   HG        .0000     .00  180.00     .00    .0000  
      IC  CE1   CG     *CD1  HD1       .0000     .00  180.00     .00    .0000  
      IC  CE2   CG     *CD2  HD2       .0000     .00  180.00     .00    .0000  
      IC  CZ    CD1    *CE1  HE1       .0000     .00  180.00     .00    .0000  
      IC  CZ    CD2    *CE2  HE2       .0000     .00  180.00     .00    .0000  
      IC  CE1   CE2    *CZ   OH        .0000     .00  180.00     .00    .0000  
 
 
 END {PHEO}
 !-----------------------------------------------------------
 
 RESIDUE INDO ! indole, adm jr. 
 
 GROUP  
       ATOM  HG    TYPE=HP   CHARGE=     .1150  END !                   HE3 
       ATOM  CG    TYPE=CY   CHARGE=    -.1450  END !                    | 
       ATOM  CD2   TYPE=CPT  CHARGE=    -.0200  END !      HG           CE3 
       ATOM  CD1   TYPE=CA   CHARGE=     .0350  END !        \         /   \ 
       ATOM  HD1   TYPE=HP   CHARGE=     .1150  END !         CG-----CD2   CZ3-HZ3 
       ATOM  NE1   TYPE=NY   CHARGE=    -.6100  END !          |      |     | 
       ATOM  HE1   TYPE=H    CHARGE=     .3800  END !         CD1    CE2   CH2-HH2 
       ATOM  CE2   TYPE=CPT  CHARGE=     .1300  END !        /   \   / \   / 
 GROUP						    !      HD1    NE1   CZ2 
       ATOM  CE3   TYPE=CA   CHARGE=    -.1150  END !              |     | 
       ATOM  HE3   TYPE=HP   CHARGE=     .1150  END !             HE1   HZ2 
 GROUP  
       ATOM  CZ2   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HZ2   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CZ3   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HZ3   TYPE=HP   CHARGE=     .1150  END  
 GROUP  
       ATOM  CH2   TYPE=CA   CHARGE=    -.1150  END  
       ATOM  HH2   TYPE=HP   CHARGE=     .1150  END  
 !END GROUP
 
      BOND  CG    HG     
      BOND  CD1   CG     
      BOND  CD2   CG     
      BOND  NE1   CD1    
      BOND  CE2   CD2    
      BOND  CZ2   CE2    
      BOND  CZ3   CE3    
      BOND  CH2   CZ2    
      BOND  CZ3   CH2    
      BOND  CD2   CE3    
      BOND  NE1   CE2    
      BOND  CD1   HD1    
      BOND  NE1   HE1    
      BOND  CE3   HE3    
      BOND  CZ2   HZ2    
      BOND  CZ3   HZ3    
      BOND  CH2   HH2    
 
      DONOR  HE1   NE1    
 
      IC  CG    CD1    NE1   CE2       .0000     .00     .00     .00    .0000  
      IC  CD2   CB     *CG   CD1       .0000     .00  180.00     .00    .0000  
      IC  CD1   CG     CD2   CE2       .0000     .00     .00     .00    .0000  
      IC  CD2   CG     CD1   NE1       .0000     .00     .00     .00    .0000  
      IC  CE2   CG     *CD2  CE3       .0000     .00  180.00     .00    .0000  
      IC  CE2   CD2    CE3   CZ3       .0000     .00     .00     .00    .0000  
      IC  CD2   CE3    CZ3   CH2       .0000     .00     .00     .00    .0000  
      IC  CE3   CZ3    CH2   CZ2       .0000     .00     .00     .00    .0000  
      IC  CZ3   CD2    *CE3  HE3       .0000     .00  180.00     .00    .0000  
      IC  CH2   CE3    *CZ3  HZ3       .0000     .00  180.00     .00    .0000  
      IC  CZ2   CZ3    *CH2  HH2       .0000     .00  180.00     .00    .0000  
      IC  CE2   CH2    *CZ2  HZ2       .0000     .00  180.00     .00    .0000  
      IC  CD1   CE2    *NE1  HE1       .0000     .00  180.00     .00    .0000  
      IC  CG    NE1    *CD1  HD1       .0000     .00  180.00     .00    .0000  
      IC  CD1   CD2    *CG   HG        .0000     .00  180.00     .00    .0000  
 
 
 END {INDO}
 !-----------------------------------------------------------
 
 RESIDUE AP2  ! acetyl-prolineamide, R. Dunbrack 
 
 GROUP        					    !HY1 HY2 HY3 
       ATOM  N     TYPE=N    CHARGE=    -.2900  END !   \ | / 
       ATOM  CA    TYPE=CP1  CHARGE=     .0200  END !    CAY 
       ATOM  CB    TYPE=CP2  CHARGE=    -.1800  END !     | 
       ATOM  CG    TYPE=CP2  CHARGE=    -.1800  END !  OY=CY  HD1 HD2 
       ATOM  CD    TYPE=CP3  CHARGE=     .0000  END !      \    \ / 
       ATOM  CY    TYPE=C    CHARGE=     .5100  END !       N---CD   HG1 
       ATOM  OY    TYPE=O    CHARGE=    -.5100  END !       |     \  / 
       ATOM  C     TYPE=CC   CHARGE=     .5100  END !       |      CG 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !       |     /  \ 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !    HA-CA--CB   HG2 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !       |   / \ 
       ATOM  HG1   TYPE=HA   CHARGE=     .0900  END !       | HB1 HB2 
       ATOM  HG2   TYPE=HA   CHARGE=     .0900  END !     O=C 
       ATOM  HD1   TYPE=HA   CHARGE=     .0900  END !       | 
       ATOM  HD2   TYPE=HA   CHARGE=     .0900  END !       NT 
       ATOM  CAY   TYPE=CT3  CHARGE=    -.2700  END !      / \ 
       ATOM  HY1   TYPE=HA   CHARGE=     .0900  END !   HT1   HT2 
       ATOM  HY2   TYPE=HA   CHARGE=     .0900  END  
       ATOM  HY3   TYPE=HA   CHARGE=     .0900  END  
       ATOM  O     TYPE=O    CHARGE=    -.5100  END  
       ATOM  NT    TYPE=NH2  CHARGE=    -.6200  END  
       ATOM  HT1   TYPE=H    CHARGE=     .3100  END  
       ATOM  HT2   TYPE=H    CHARGE=     .3100  END  
 !END GROUP
 
      BOND  CY    N      
      BOND  CAY   HY1    
      BOND  CAY   HY2    
      BOND  CAY   HY3    
      BOND  CY    CAY    
      BOND  CY    OY     
      BOND  N     CA     
      BOND  CA    CB     
      BOND  CB    CG     
      BOND  CG    CD     
      BOND  CD    N      
      BOND  HA    CA     
      BOND  HG1   CG     
      BOND  HG2   CG     
      BOND  HD1   CD     
      BOND  HD2   CD     
      BOND  HB1   CB     
      BOND  HB2   CB     
      BOND  C     O      
      BOND  C     CA     
      BOND  C     NT     
      BOND  NT    HT1    
      BOND  NT    HT2    
 
      DIHEDRAL CB  CA  C   NT  ! multiple dihedral
      DIHEDRAL CB  CA  C   NT  ! multiple dihedral
      DIHEDRAL HA  CA  C   NT  ! multiple dihedral
      DIHEDRAL HA  CA  C   NT  ! multiple dihedral
      DIHEDRAL N   CA  C   NT  ! multiple dihedral
      DIHEDRAL N   CA  C   NT  ! multiple dihedral
      DIHEDRAL CB  CA  C   O   ! multiple dihedral
      DIHEDRAL CB  CA  C   O   ! multiple dihedral
      DIHEDRAL HA  CA  C   O   ! multiple dihedral
      DIHEDRAL HA  CA  C   O   ! multiple dihedral
      DIHEDRAL  CAY   CY    N     CA  ! multiple dihedral
      DIHEDRAL  CAY   CY    N     CA  ! multiple dihedral
      DIHEDRAL  CAY   CY    N     CD  ! multiple dihedral
      DIHEDRAL  CAY   CY    N     CD  ! multiple dihedral
      DIHEDRAL  OY    CY    N     CA  ! multiple dihedral
      DIHEDRAL  OY    CY    N     CA  ! multiple dihedral
      DIHEDRAL  OY    CY    N     CD  ! multiple dihedral
      DIHEDRAL  OY    CY    N     CD  ! multiple dihedral

      IMPROPER  CY    CAY   N     OY     
      IMPROPER  C     CA    NT    O      
      IMPROPER  C     NT    CA    O      
      IMPROPER  N     CY    CA    CD     
      IMPROPER  NT    C     HT2   HT1    
 
      DONOR  HT1   NT     
      DONOR  HT2   NT     
 
      ACCEPTOR  OY    CY     
      ACCEPTOR  O     C

      ! Improper ICs 
 
      IC  CY    CA     *N    CD        .0000     .00  180.00     .00    .0000  
      IC  N     C      *CA   CB        .0000     .00  120.00     .00    .0000  
      IC  N     C      *CA   HA        .0000     .00 -120.00     .00    .0000  
      IC  N     CAY    *CY   OY        .0000     .00  180.00     .00    .0000
      ! Backbone + Peptide bond IC's 
      IC  OY    CY     N     CA        .0000     .00     .00     .00    .0000 !Omega 
      IC  CAY   CY     N     CA        .0000     .00  180.00     .00    .0000 !Omega 
      IC  CY    N      CA    C         .0000  120.64  -60.00     .00    .0000 !Psi 
      IC  N     CA     C     NT        .0000     .00  180.00     .00    .0000 !Phi 
      IC  CAY   CY     N     CD        .0000     .00     .00  120.64    .0000  
      IC  CY    N      CD    CG        .0000     .00  168.60  103.28    .0000  
      IC  CY    N      CA    CB        .0000     .00  168.64  103.34    .0000
      ! Ring IC's 
      IC  N     CA     CB    CG        .0000  103.34   29.25  103.67    .0000  
      IC  CA    CB     CG    CD        .0000  103.67  -36.72  103.63    .0000  
      IC  CB    CG     CD    N         .0000  103.63   29.50  103.28    .0000  
      IC  CG    CD     N     CA        .0000  103.28  -11.61  112.90    .0000  
      IC  CD    N      CA    CB        .0000  112.90  -11.17  103.34    .0000
      ! Carbonyl IC's 
      IC  CD    N      CA    C         .0000  112.90  120.00  108.00    .0000
      ! Hydrogen IC's 
      IC  CD    N      CA    HA        .0000     .00 -120.00  108.00    .0000  
      IC  N     CA     CB    HB1       .0000     .00  120.00  108.00    .0000  
      IC  N     CA     CB    HB2       .0000     .00 -120.00  108.00    .0000  
      IC  CA    CB     CG    HG1       .0000     .00  120.00  108.00    .0000  
      IC  CA    CB     CG    HG2       .0000     .00 -120.00  108.00    .0000  
      IC  CB    CG     CD    HD1       .0000     .00  120.00  108.00    .0000  
      IC  CB    CG     CD    HD2       .0000     .00 -120.00  108.00    .0000  
      IC  NT    CA     *C    O         .0000     .00  180.00     .00    .0000  
      IC  CA    C      NT    HT2       .0000     .00  180.00     .00    .0000  
      IC  C     HT2    *NT   HT1       .0000     .00  180.00     .00    .0000  
 
 
 END {AP2 }
 !-----------------------------------------------------------
 
 RESIDUE TP2  ! prolineamide, R. Dunbrack 
 
 GROUP  
       ATOM  N     TYPE=NP   CHARGE=    -.0700  END ! 
       ATOM  HN1   TYPE=HC   CHARGE=     .2400  END ! 
       ATOM  HN2   TYPE=HC   CHARGE=     .2400  END ! 
       ATOM  CD    TYPE=CP3  CHARGE=     .1600  END ! 
       ATOM  CB    TYPE=CP2  CHARGE=    -.1800  END !    HN1   HD1 HD2 
       ATOM  CG    TYPE=CP2  CHARGE=    -.1800  END !      \    \ / 
       ATOM  CA    TYPE=CP1  CHARGE=     .1600  END !  HN2--N---CD   HG1 
       ATOM  C     TYPE=CC   CHARGE=     .5100  END !       |     \  / 
       ATOM  O     TYPE=O    CHARGE=    -.5100  END !       |      CG 
       ATOM  HA    TYPE=HB   CHARGE=     .0900  END !       |     /  \ 
       ATOM  HB1   TYPE=HA   CHARGE=     .0900  END !    HA-CA--CB   HG2 
       ATOM  HB2   TYPE=HA   CHARGE=     .0900  END !       |   / \ 
       ATOM  HG1   TYPE=HA   CHARGE=     .0900  END !       | HB1 HB2 
       ATOM  HG2   TYPE=HA   CHARGE=     .0900  END !     O=C 
       ATOM  HD1   TYPE=HA   CHARGE=     .0900  END !       | 
       ATOM  HD2   TYPE=HA   CHARGE=     .0900  END !       NT 
       ATOM  NT    TYPE=NH2  CHARGE=    -.6200  END !      / \ 
       ATOM  HT1   TYPE=H    CHARGE=     .3100  END !   HT1   HT2 
       ATOM  HT2   TYPE=H    CHARGE=     .3100  END ! 
 !END GROUP
 
      BOND  HN1   N      
      BOND  HN2   N      
      BOND  N     CA     
      BOND  CA    CB     
      BOND  CB    CG     
      BOND  CG    CD     
      BOND  CD    N      
      BOND  HA    CA     
      BOND  HG1   CG     
      BOND  HG2   CG     
      BOND  HD1   CD     
      BOND  HD2   CD     
      BOND  HB1   CB     
      BOND  HB2   CB     
      BOND  C     O      
      BOND  C     CA     
      BOND  C     NT     
      BOND  NT    HT1    
      BOND  NT    HT2    

      DIHEDRAL CB  CA  C   NT  ! multiple dihedral
      DIHEDRAL CB  CA  C   NT  ! multiple dihedral
      DIHEDRAL HA  CA  C   NT  ! multiple dihedral
      DIHEDRAL HA  CA  C   NT  ! multiple dihedral
      DIHEDRAL CB  CA  C   O   ! multiple dihedral
      DIHEDRAL CB  CA  C   O   ! multiple dihedral
      DIHEDRAL HA  CA  C   O   ! multiple dihedral
      DIHEDRAL HA  CA  C   O   ! multiple dihedral

      IMPROPER  C     CA    NT    O      
      IMPROPER  C     NT    CA    O      
      IMPROPER  NT    C     HT2   HT1    
 
      DONOR  HT1   NT     
      DONOR  HT2   NT     
      DONOR  HN1   N      
      DONOR  HN2   N      
 
      ACCEPTOR  O     C

      ! Improper ICs 
 
      IC  HN1   CA     *N    CD        .0000     .00  120.00     .00    .0000  
      IC  HN2   CA     *N    HN1       .0000     .00  120.00     .00    .0000  
      IC  N     C      *CA   CB        .0000     .00  120.00     .00    .0000  
      IC  N     C      *CA   HA        .0000     .00 -120.00     .00    .0000
      ! Backbone + Peptide bond IC's 
      IC  N     CA     C     NT        .0000     .00  180.00     .00    .0000 !Phi 
      ! Ring IC's 
      IC  N     CA     CB    CG        .0000  103.34   29.25  103.67    .0000  
      IC  CA    CB     CG    CD        .0000  103.67  -36.72  103.63    .0000  
      IC  CB    CG     CD    N         .0000  103.63   29.50  103.28    .0000  
      IC  CG    CD     N     CA        .0000  103.28  -11.61  112.90    .0000  
      IC  CD    N      CA    CB        .0000  112.90  -11.17  103.34    .0000
      ! Carbonyl IC's 
      IC  CD    N      CA    C         .0000  112.90  120.00  108.00    .0000
      ! Hydrogen IC's 
      IC  CD    N      CA    HA        .0000     .00 -120.00  108.00    .0000  
      IC  N     CA     CB    HB1       .0000     .00  120.00  108.00    .0000  
      IC  N     CA     CB    HB2       .0000     .00 -120.00  108.00    .0000  
      IC  CA    CB     CG    HG1       .0000     .00  120.00  108.00    .0000  
      IC  CA    CB     CG    HG2       .0000     .00 -120.00  108.00    .0000  
      IC  CB    CG     CD    HD1       .0000     .00  120.00  108.00    .0000  
      IC  CB    CG     CD    HD2       .0000     .00 -120.00  108.00    .0000  
      IC  NT    CA     *C    O         .0000     .00  180.00     .00    .0000  
      IC  CA    C      NT    HT2       .0000     .00  180.00     .00    .0000  
      IC  C     HT2    *NT   HT1       .0000     .00  180.00     .00    .0000  
 
 
 END {TP2 }
 !-----------------------------------------------------------
 
 RESIDUE ETHA ! ethane, S. Fischer 
 
 GROUP  
       ATOM  H1    TYPE=HA   CHARGE=     .0900  END !   H1      H4 
       ATOM  H2    TYPE=HA   CHARGE=     .0900  END !    \      / 
       ATOM  H3    TYPE=HA   CHARGE=     .0900  END !  H2-C1--C2-H5 
       ATOM  C1    TYPE=CT3  CHARGE=    -.2700  END !    /      \ 
 GROUP						    !   H3       H6 
       ATOM  H4    TYPE=HA   CHARGE=     .0900  END  
       ATOM  H5    TYPE=HA   CHARGE=     .0900  END  
       ATOM  H6    TYPE=HA   CHARGE=     .0900  END  
       ATOM  C2    TYPE=CT3  CHARGE=    -.2700  END  
 !END GROUP
 
      BOND  C1    H1     
      BOND  C1    H2     
      BOND  C1    H3     
      BOND  C1    C2     
      BOND  C2    H4     
      BOND  C2    H5     
      BOND  C2    H6     
 
      IC  H1    C1     C2    H4        .0000     .00     .00     .00    .0000  
      IC  H1    C1     C2    H5        .0000     .00  120.00     .00    .0000  
      IC  H1    C1     C2    H6        .0000     .00  240.00     .00    .0000  
      IC  H2    C1     C2    H6        .0000     .00  120.00     .00    .0000  
      IC  H3    C1     C2    H6        .0000     .00  240.00     .00    .0000  
 
 
 END {ETHA}
 !-----------------------------------------------------------
 
 RESIDUE PROP ! propane, adm jr. 
 
 GROUP ! 
       ATOM  H11   TYPE=HA   CHARGE=     .0900  END !  H11   H21    H31 
       ATOM  H12   TYPE=HA   CHARGE=     .0900  END !    \    |     / 
       ATOM  H13   TYPE=HA   CHARGE=     .0900  END ! H12-C1--C2--C3-H32 
       ATOM  C1    TYPE=CT3  CHARGE=    -.2700  END !    /    |     \ 
 GROUP						    !  H13   H22    H33 
       ATOM  C2    TYPE=CT2  CHARGE=    -.1800  END ! 
       ATOM  H21   TYPE=HA   CHARGE=     .0900  END ! 
       ATOM  H22   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  H31   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H32   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H33   TYPE=HA   CHARGE=     .0900  END  
       ATOM  C3    TYPE=CT3  CHARGE=    -.2700  END  
 !END GROUP
 
      BOND  C1    H11    
      BOND  C1    H12    
      BOND  C1    H13    
      BOND  C1    C2     
      BOND  C2    H21    
      BOND  C2    H22    
      BOND  C2    C3     
      BOND  C3    H31    
      BOND  C3    H32    
      BOND  C3    H33    
 
      IC  H11   C1     C2    C3        .0000     .00     .00     .00    .0000  
      IC  H11   C1     C2    H21       .0000     .00  120.00     .00    .0000  
      IC  H11   C1     C2    H22       .0000     .00  240.00     .00    .0000  
      IC  C1    C2     C3    H31       .0000     .00     .00     .00    .0000  
      IC  C1    C2     C3    H32       .0000     .00  120.00     .00    .0000  
      IC  C1    C2     C3    H33       .0000     .00  240.00     .00    .0000  
      IC  C3    C2     C1    H12       .0000     .00  120.00     .00    .0000  
      IC  C3    C2     C1    H13       .0000     .00  120.00     .00    .0000  
 
 
 END {PROP}
 !-----------------------------------------------------------
 
 RESIDUE BUTA ! butane, S. Fischer 
 
 GROUP  
       ATOM  H11   TYPE=HA   CHARGE=     .0900  END !  H11   H21 H31    H41 
       ATOM  H12   TYPE=HA   CHARGE=     .0900  END !    \    |   |     / 
       ATOM  H13   TYPE=HA   CHARGE=     .0900  END ! H12-C1--C2--C3--C4-H42 
       ATOM  C1    TYPE=CT3  CHARGE=    -.2700  END !    /    |   |     \ 
 GROUP						    !  H13   H22 H33    H43 
       ATOM  H21   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H22   TYPE=HA   CHARGE=     .0900  END  
       ATOM  C2    TYPE=CT2  CHARGE=    -.1800  END  
 GROUP  
       ATOM  H31   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H32   TYPE=HA   CHARGE=     .0900  END  
       ATOM  C3    TYPE=CT2  CHARGE=    -.1800  END  
 GROUP  
       ATOM  H41   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H42   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H43   TYPE=HA   CHARGE=     .0900  END  
       ATOM  C4    TYPE=CT3  CHARGE=    -.2700  END  
 !END GROUP
 
      BOND  H11   C1     
      BOND  H12   C1     
      BOND  H13   C1     
      BOND  C1    C2     
      BOND  H21   C2     
      BOND  H22   C2     
      BOND  C2    C3     
      BOND  H31   C3     
      BOND  H32   C3     
      BOND  C3    C4     
      BOND  H41   C4     
      BOND  H42   C4     
      BOND  H43   C4     
 
      IC  H11   C1     C2    C3        .0000     .00     .00     .00    .0000  
      IC  H11   C1     C2    H21       .0000     .00  120.00     .00    .0000  
      IC  H11   C1     C2    H22       .0000     .00  240.00     .00    .0000  
      IC  H12   C1     C2    C3        .0000     .00  120.00     .00    .0000  
      IC  H13   C1     C2    C3        .0000     .00  240.00     .00    .0000  
      IC  C1    C2     C3    C4        .0000     .00     .00     .00    .0000  
      IC  C1    C2     C3    H31       .0000     .00  120.00     .00    .0000  
      IC  C1    C2     C3    H32       .0000     .00  240.00     .00    .0000  
      IC  H21   C2     C3    C4        .0000     .00  120.00     .00    .0000  
      IC  H22   C2     C3    C4        .0000     .00  240.00     .00    .0000  
      IC  C2    C3     C4    H41       .0000     .00     .00     .00    .0000  
      IC  C2    C3     C4    H42       .0000     .00  120.00     .00    .0000  
      IC  C2    C3     C4    H43       .0000     .00  240.00     .00    .0000  
      IC  H31   C3     C4    H43       .0000     .00  120.00     .00    .0000  
      IC  H32   C3     C4    H43       .0000     .00  240.00     .00    .0000  
 
 
 END {BUTA}
 !-----------------------------------------------------------
 
 RESIDUE IBUT ! Iso-butane, S. Fischer 
 
 GROUP  
       ATOM  CT    TYPE=CT1  CHARGE=    -.0900  END !          H12 
       ATOM  HT    TYPE=HA   CHARGE=     .0900  END !           | 
 GROUP						    !       H11-C1-H13 
       ATOM  C1    TYPE=CT3  CHARGE=    -.2700  END !           | 
       ATOM  H11   TYPE=HA   CHARGE=     .0900  END !           CT-HT 
       ATOM  H12   TYPE=HA   CHARGE=     .0900  END !          / \ 
       ATOM  H13   TYPE=HA   CHARGE=     .0900  END !         /   \ 
 GROUP						    !   H21-C2     C3-H31 
       ATOM  C2    TYPE=CT3  CHARGE=    -.2700  END !      / |     | \ 
       ATOM  H21   TYPE=HA   CHARGE=     .0900  END !   H22 H23   H33 H32 
       ATOM  H22   TYPE=HA   CHARGE=     .0900  END ! 
       ATOM  H23   TYPE=HA   CHARGE=     .0900  END ! 
 GROUP  
       ATOM  C3    TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  H31   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H32   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H33   TYPE=HA   CHARGE=     .0900  END  
 !END GROUP
 
      BOND  CT    C1     
      BOND  CT    C2     
      BOND  CT    C3     
      BOND  CT    HT     
      BOND  C1    H11    
      BOND  C1    H12    
      BOND  C1    H13    
      BOND  C2    H21    
      BOND  C2    H22    
      BOND  C2    H23    
      BOND  C3    H31    
      BOND  C3    H32    
      BOND  C3    H33    
 
      IC  HT    CT     C1    H11       .0000     .00   60.00     .00    .0000  
      IC  CT    C1     H11   H12       .0000     .00     .00     .00    .0000  
      IC  HT    CT     C1    H12       .0000     .00  180.00     .00    .0000  
      IC  HT    CT     C1    H13       .0000     .00  300.00     .00    .0000  
      IC  H12   C1     CT    C2        .0000     .00   60.00     .00    .0000  
      IC  H12   C1     CT    C3        .0000     .00  -60.00     .00    .0000  
      IC  HT    CT     C2    H21       .0000     .00   60.00     .00    .0000  
      IC  HT    CT     C2    H22       .0000     .00  180.00     .00    .0000  
      IC  HT    CT     C2    H23       .0000     .00  300.00     .00    .0000  
      IC  HT    CT     C3    H31       .0000     .00   60.00     .00    .0000  
      IC  HT    CT     C3    H32       .0000     .00  180.00     .00    .0000  
      IC  HT    CT     C3    H33       .0000     .00  300.00     .00    .0000  
 
 
 END {IBUT}
 !-----------------------------------------------------------
 
 RESIDUE PENT ! pentane, adm jr. 
 
 GROUP  
       ATOM  C1    TYPE=CT3  CHARGE=    -.2700  END !  H11   H21 H31 H41    H51 
       ATOM  H11   TYPE=HA   CHARGE=     .0900  END !    \    |   |   |     / 
       ATOM  H12   TYPE=HA   CHARGE=     .0900  END ! H12-C1--C2--C3--C4--C5-H52 
       ATOM  H13   TYPE=HA   CHARGE=     .0900  END !    /    |   |   |     \ 
 GROUP						    !  H13   H22 H33 H42    H53 
       ATOM  C2    TYPE=CT2  CHARGE=    -.1800  END  
       ATOM  H21   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H22   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  C3    TYPE=CT2  CHARGE=    -.1800  END  
       ATOM  H31   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H32   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  C4    TYPE=CT2  CHARGE=    -.1800  END  
       ATOM  H41   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H42   TYPE=HA   CHARGE=     .0900  END  
 GROUP  
       ATOM  C5    TYPE=CT3  CHARGE=    -.2700  END  
       ATOM  H51   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H52   TYPE=HA   CHARGE=     .0900  END  
       ATOM  H53   TYPE=HA   CHARGE=     .0900  END  
 !END GROUP
 
      BOND  C1    C2     
      BOND  C2    C3     
      BOND  C3    C4     
      BOND  C4    C5     
      BOND  C1    H11    
      BOND  C1    H12    
      BOND  C1    H13    
      BOND  C2    H21    
      BOND  C2    H22    
      BOND  C3    H31    
      BOND  C3    H32    
      BOND  C4    H41    
      BOND  C4    H42    
      BOND  C5    H51    
      BOND  C5    H52    
      BOND  C5    H53    
 
      IC  H11   C1     C2    C3        .0000     .00  180.00     .00    .0000  
      IC  C1    C2     C3    C4        .0000     .00  180.00     .00    .0000  
      IC  C2    C3     C4    C5        .0000     .00  180.00     .00    .0000  
      IC  C3    C4     C5    H51       .0000     .00  180.00     .00    .0000  
      IC  H11   C2     *C1   H12       .0000     .00  120.00     .00    .0000  
      IC  H11   C2     *C1   H13       .0000     .00 -120.00     .00    .0000  
      IC  C1    C3     *C2   H21       .0000     .00  120.00     .00    .0000  
      IC  C1    C3     *C2   H22       .0000     .00 -120.00     .00    .0000  
      IC  C2    C4     *C3   H31       .0000     .00  120.00     .00    .0000  
      IC  C2    C4     *C3   H32       .0000     .00 -120.00     .00    .0000  
      IC  C3    C5     *C4   H41       .0000     .00  120.00     .00    .0000  
      IC  C3    C5     *C4   H42       .0000     .00 -120.00     .00    .0000  
      IC  C4    H51    *C5   H52       .0000     .00  120.00     .00    .0000  
      IC  C4    H51    *C5   H53       .0000     .00 -120.00     .00    .0000  
 
 
 END {PENT}
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
 
 PRESIDUE FHEM ! fix up the heme by deleting unwanted autogenerated angles 
	       ! unliganded heme patch, K. Kuczera 
 
 
 DELE ANGLE  1NA   1FE   1NC    
 DELE ANGLE  1NB   1FE   1ND    
 
 
 END {FHEM}
 !-----------------------------------------------------------
 
 RESIDUE TIP3 ! TIPS3P WATER MODEL, adm jr. 
 
 GROUP ! Generate noangle nodihedral 
       ATOM  OH2   TYPE=OT   CHARGE=    -.8340  END  
       ATOM  H1    TYPE=HT   CHARGE=     .4170  END  
       ATOM  H2    TYPE=HT   CHARGE=     .4170  END  
 !END GROUP
 
      BOND  OH2   H1     
      BOND  OH2   H2     
 
      ANGLE  H1    OH2   H2     
 
      ACCEPTOR  OH2   NONE   
 
 
 END {TIP3}
 
 SET ECHO=TRUE END 
