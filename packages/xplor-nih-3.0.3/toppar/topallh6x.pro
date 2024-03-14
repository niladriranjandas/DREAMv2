remarks  Residue Topology File for Proteins Using All Atom Hydrogens (ALLH)
remarks  Cambridge notation for atom names
remarks  atom type 3    ha  added for aliphatic hydrogens
remarks  atom type 16   ct  added for tetrahedral carbons
remarks  Impropers rearranged to match rtoph7 for gaunido & sidechain amide n
remarks  2nd donor antecedants converted to CA for peptides CD for arg HE
remarks  Amide charges from Hayes and Kollman JACS 98:3335,7811 (1796)
remarks  HA charges from Wiberg and Wendoloski J. Comp. Chem. 2:53 (1981)
remarks  ...
remarks  modified for NMR refinement:
remarks  improper term included to force all amino acids in L form
set echo=false end
autogenerate angles=true end

!
! The charges in this file were set up using the experimental formamide
! and N-methyl acetamide dipole moments of 3.7 D through the following
! proceedure.
! 	1 - the C=O and N-H2 of formamide were required to be neutral
! 	2 - the HA charge was arbitrarily set to zero (various Mulliken
! 		population analyses give it small magnitude and variable
! 		sign).
! 	3 - The above conditions and the direction of the formamide dipole
! 		competely determine the formamide charges
! 	4 - The CA peptide charge was arbitrarily set to 0.1 (various
! 		Mulliken population analysis give it 0.0 to 0.12 charges).
! 	5 - The peptide charges were obtained by keeping the same C=O
! 		and HN charges.  This requires a readjustment of the
! 		N charge to maintain neutrality.
! 		
! The resulting set of charges gives close to the same dipole moment
! when applied to N-methyl acetamide, but a somewhat different
! direction.
!
MASS     H      1.00800
MASS     HC     1.00800
MASS     HA     1.00800
MASS     LP     0.0                          
MASS     C     12.01100
MASS     CT    12.01100
MASS     CM    12.01100
MASS     N     14.00700
MASS     NR    14.00700
MASS     NP    14.00700
MASS     NH1   14.00700
MASS     NH2   14.00700
MASS     NH3   14.00700
MASS     NC2   14.00700
MASS     O     15.99900
MASS     OC    15.99900
MASS     OM    15.99940
MASS     OH1   15.99900
MASS     OH2   15.99900
MASS     OH1E  17.00700
MASS     OH2E  18.01500
MASS     S     32.06000
MASS     SH1E  33.06800
MASS     FE    55.84700



RESIdue ACE
GROUP
 ATOM CA   TYPE=CT  CHARge=-0.30000 END
 ATOM HA1  TYPE=HA  CHARge=0.1000   END
 ATOM HA2  TYPE=HA  CHARge=0.1000   END
 ATOM HA3  TYPE=HA  CHARge=0.1000   END
GROUP
 ATOM C    TYPE=C   CHARge=0.48000  END
 ATOM O    TYPE=O   CHARge=-0.480   END
BOND C CA     BOND O C    BOND CA HA1     BOND CA HA2    BOND CA HA3
DIHE O    C    CA   HA3 
ACCE O    C
END
!======================================================================

RESI ALA 
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360  END
 ATOM HN   TYPE=H   CHARge=0.260   END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000   END
 ATOM HA   TYPE=HA  CHARge=0.100   END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.300  END
 ATOM HB1  TYPE=HA  CHARge=0.100   END
 ATOM HB2  TYPE=HA  CHARge=0.100   END
 ATOM HB3  TYPE=HA  CHARge=0.100   END
GROUP
 ATOM C    TYPE=C   CHARge=0.480   END
 ATOM O    TYPE=O   CHARge=-0.480  END
BOND CB CA   BOND N HN    BOND N CA    BOND O C
BOND C CA    BOND CA HA   BOND CB HB1  BOND CB HB2   BOND CB HB3
DIHE HB3  CB   CA   N     DIHE HB2  CB   CA   N      DIHE HB1  CB   CA   N
IMPRoper CA   N    C    CB  !tetrahedral CA
DONO HN   N
ACCE O    C
END
!======================================================================


RESI ARG           ! Arginine with idealized charges on the guanido
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.200 END
 ATOM HB1  TYPE=HA  CHARge=0.100  END
 ATOM HB2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CG   TYPE=CT  CHARge=-0.200 END
 ATOM HG1  TYPE=HA  CHARge=0.100  END
 ATOM HG2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CD   TYPE=CT  CHARge=-0.100 END
 ATOM HD1  TYPE=HA  CHARge=0.150  END
 ATOM HD2  TYPE=HA  CHARge=0.150  END
GROUP
 ATOM NE   TYPE=NH1 CHARge=-0.600 END
 ATOM HE   TYPE=H   CHARge=0.400  END
 ATOM CZ   TYPE=C   CHARge=0.600  END
GROUP
 ATOM NH1  TYPE=NC2 CHARge=-0.600 END
 ATOM HH11 TYPE=HC  CHARge=0.400  END
 ATOM HH12 TYPE=HC  CHARge=0.400  END
GROUP
 ATOM NH2  TYPE=NC2 CHARge=-0.600 END
 ATOM HH21 TYPE=HC  CHARge=0.400  END
 ATOM HH22 TYPE=HC  CHARge=0.400  END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND CB CA     BOND CG CB     BOND CD CG    BOND NE CD     BOND CZ NE
BOND NH1 CZ    BOND NH2 CZ    BOND N HN     BOND N CA
BOND O C       BOND C CA      BOND CA HA    BOND CB HB1
BOND CB HB2    BOND CG HG1    BOND CG HG2   BOND CD HD1    BOND CD HD2
BOND NE HE     BOND NH1 HH11  BOND NH1 HH12 BOND NH2 HH21  BOND NH2 HH22
DIHE N    CA   CB   CG     DIHE CD   CG   CB   CA     DIHE NE   CD   CG   CB
DIHE CZ   NE   CD   CG     DIHE NH1  CZ   NH2  HH21   DIHE NE   CZ   NH2  HH22
DIHE NH2  CZ   NH1  HH11   DIHE NE   CZ   NH1  HH12
IMPR NH1  HH12 HH11 CZ     IMPR NH2  HH22 HH21 CZ
IMPR NE   CD   CZ   HE     IMPR CZ   NH1  NH2  NE
IMPRoper CA   N    C    CB  !tetrahedral CA
DONO HN   N
DONO HE   NE
DONO HH11 NH1
DONO HH12 NH1
DONO HH21 NH2
DONO HH22 NH2
ACCE O    C
END
!======================================================================


RESI ASN
GROUP
 ATOM N    TYPE=NH1  CHARge=-0.360 END
 ATOM HN   TYPE=H    CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT   CHARge=0.000  END
 ATOM HA   TYPE=HA   CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT   CHARge=-0.200 END
 ATOM HB1  TYPE=HA   CHARge=0.100  END
 ATOM HB2  TYPE=HA   CHARge=0.100  END
GROUP
 ATOM CG   TYPE=C    CHARge=0.480  END
 ATOM OD1  TYPE=O    CHARge=-0.480 END
GROUP
 ATOM ND2  TYPE=NH2  CHARge=-0.520 END
 ATOM HD21 TYPE=H    CHARge=0.260  END
 ATOM HD22 TYPE=H    CHARge=0.260  END
GROUP
 ATOM C    TYPE=C    CHARge=0.480  END
 ATOM O    TYPE=O    CHARge=-0.480 END
BOND CB CA    BOND CG CB    BOND OD1 CG      BOND ND2 CG       
BOND N HN     BOND N CA     BOND O C         BOND C CA  
BOND CA HA    BOND CB HB1   BOND CB HB2      BOND ND2 HD21   BOND ND2 HD22
DIHE CG   CB   CA   N         DIHE OD1  CG   CB   CA
DIHE OD1  CG   ND2  HD22      DIHE CB   CG   ND2  HD21  !! INCLUDED, MN
IMPRoper CA   N    C    CB  !tetrahedral CA
IMPR CG   CB   OD1  ND2       IMPR ND2  HD21 HD22 CG
DONO HN   N
DONO HD21 ND2
DONO HD22 ND2
ACCE OD1  CG
ACCE O    C
END
!======================================================================


RESI ASP
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.450 END
 ATOM HB1  TYPE=HA  CHARge=0.100  END
 ATOM HB2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CG   TYPE=C   CHARge=0.490  END
 ATOM OD1  TYPE=OC  CHARge=-0.620 END
 ATOM OD2  TYPE=OC  CHARge=-0.620 END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND CB CA     BOND CG CB    BOND OD1 CG     BOND OD2 CG 
BOND N HN      BOND N CA     BOND O C        BOND C CA   
BOND CA HA     BOND CB HB1   BOND CB HB2
DIHE CG   CB   CA   N        DIHE OD1  CG   CB   CA
IMPRoper CA   N    C    CB  !tetrahedral CA
IMPR CG   CB   OD1  OD2      
DONO HN   N
ACCE OD1  CG
ACCE OD2  CG
ACCE O    C
END
!======================================================================


RESI CYS
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.200 END
 ATOM HB1  TYPE=HA  CHARge=0.100  END
 ATOM HB2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM SG   TYPE=S   CHARge=-0.050 END
 ATOM HG1  TYPE=H   CHARge=0.050  END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND CB CA        BOND SG CB     BOND N HN    BOND N CA
BOND O C          BOND C CA      BOND CA HA   BOND CB HB1
BOND CB HB2       BOND SG HG1
DIHE SG   CB   CA   N         DIHE HG1  SG   CB   CA
IMPRoper CA   N    C    CB  !tetrahedral CA
DONO HN   N
DONO HG1  SG
!!! ACCE SG    !REMOVED, ATB
ACCE O    C
END
!======================================================================


RESI GLU
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.200 END
 ATOM HB1  TYPE=HA  CHARge=0.100  END
 ATOM HB2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CG   TYPE=CT  CHARge=-0.450 END
 ATOM HG1  TYPE=HA  CHARge=0.100  END
 ATOM HG2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CD   TYPE=C   CHARge=0.490  END
 ATOM OE1  TYPE=OC  CHARge=-0.620 END
 ATOM OE2  TYPE=OC  CHARge=-0.620 END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND CB CA      BOND CG CB     BOND CD CG     BOND OE1 CD     BOND OE2 CD
BOND N HN       BOND N CA      BOND O C       BOND C CA
BOND CA HA      BOND CB HB1    BOND CB HB2    BOND CG HG1
BOND CG HG2
DIHE CG   CB   CA   N     DIHE CD   CG   CB   CA     DIHE OE1  CD   CG   CB
IMPR CD   CG   OE1  OE2 
IMPRoper CA   N    C    CB  !tetrahedral CA
DONO HN   N
ACCE OE1  CD
ACCE OE2  CD
ACCE O    C
END
!======================================================================


RESI GLN
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.200 END
 ATOM HB1  TYPE=HA  CHARge=0.100  END
 ATOM HB2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CG   TYPE=CT  CHARge=-0.200 END
 ATOM HG1  TYPE=HA  CHARge=0.100  END
 ATOM HG2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CD   TYPE=C   CHARge=0.480  END
 ATOM OE1  TYPE=O   CHARge=-0.480 END
GROUP
 ATOM NE2  TYPE=NH2 CHARge=-0.520 END
 ATOM HE21 TYPE=H   CHARge=0.260  END
 ATOM HE22 TYPE=H   CHARge=0.260  END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND CB CA    BOND CG CB     BOND CD CG    BOND OE1 CD     BOND NE2 CD
BOND N HN     BOND N CA      BOND O C      BOND C CA
BOND CA HA    BOND CB HB1    BOND CB HB2   BOND CG HG1
BOND CG HG2   BOND NE2 HE21  BOND NE2 HE22
DIHE CG   CB   CA   N      DIHE CD   CG   CB   CA      DIHE OE1  CD   CG   CB
DIHE OE1  CD   NE2  HE22   DIHE CG   CD   NE2  HE21
IMPR CD   CG   OE1  NE2
IMPR NE2  HE21 HE22 CD
IMPRoper CA   N    C    CB  !tetrahedral CA
DONO HN   N
DONO HE21 NE2
DONO HE22 NE2
ACCE OE1  CD
ACCE O    C
END
!======================================================================


RESI GLY  
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=-0.100 END
 ATOM HA1  TYPE=HA  CHARge=0.100  END
 ATOM HA2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND N HN   BOND N CA  BOND O C   BOND C CA   BOND CA HA1   BOND CA HA2
DONO HN   N
ACCE O    C
END
!======================================================================


RESI HIS     ! Histidine with idealized charges on the ring
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.200 END
 ATOM HB1  TYPE=HA  CHARge=0.100  END
 ATOM HB2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CG   TYPE=C   CHARge=0.050  END
GROUP
 ATOM ND1  TYPE=NH1 CHARge=-0.400 END
 ATOM HD1  TYPE=H   CHARge=0.400  END
GROUP
 ATOM CD2  TYPE=C   CHARge=-0.140 END
 ATOM HD2  TYPE=HA  CHARge=0.140  END
GROUP
 ATOM CE1  TYPE=C   CHARge=-0.140 END
 ATOM HE1  TYPE=HA  CHARge=0.140  END
GROUP
 ATOM NE2  TYPE=NR  CHARge=-0.050 END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND CB CA    BOND CG CB    BOND ND1 CG     BOND CD2 CG    BOND CE1 ND1
BOND NE2 CD2  BOND N HN     BOND N CA       BOND O C
BOND C CA     BOND NE2 CE1  BOND CA HA      BOND CB HB1
BOND CB HB2   BOND ND1 HD1  BOND CD2 HD2    BOND CE1 HE1
DIHE CG   CB   CA   N       DIHE ND1  CG   CB   CA
IMPR CE1  ND1  CG   CD2     IMPR NE2  CD2  CG   ND1    IMPR NE2  CE1  ND1  CG
IMPR CE1  NE2  CD2  CG      IMPR CD2  NE2  CE1  ND1    
IMPR CG   CB   ND1  CD2     IMPR ND1  CG   CE1  HD1    IMPR CD2  CG   NE2  HD2
IMPR CE1  ND1  NE2  HE1
IMPRoper CA   N    C    CB  !tetrahedral CA
DONO HN   N
DONO HD1  ND1
ACCE NE2  " "
ACCE O    C
END
!======================================================================


RESI ILE 
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.100 END
 ATOM HB   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CG1  TYPE=CT  CHARge=-0.200 END
 ATOM HG11 TYPE=HA  CHARge=0.100  END
 ATOM HG12 TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CG2  TYPE=CT  CHARge=-0.300 END
 ATOM HG21 TYPE=HA  CHARge=0.100  END
 ATOM HG22 TYPE=HA  CHARge=0.100  END
 ATOM HG23 TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CD   TYPE=CT  CHARge=-0.300 END
 ATOM HD1  TYPE=HA  CHARge=0.100  END
 ATOM HD2  TYPE=HA  CHARge=0.100  END
 ATOM HD3  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND CB CA     BOND CG1 CB     BOND CG2 CB      BOND CD CG1
BOND N HN      BOND N CA      BOND O C         BOND C CA        
BOND CA HA     BOND CB HB     BOND CG1 HG11    BOND CG1 HG12    BOND CG2 HG21
BOND CG2 HG22  BOND CG2 HG23  BOND CD HD1      BOND CD  HD2     BOND CD HD3
IMPRoper CA   N    C    CB  !tetrahedral CA
DIHE CG1  CB   CA   N     DIHE CD   CG1  CB   CA
DIHE HG23 CG2  CB   CA    DIHE HG22 CG2  CB   CA      DIHE HG21 CG2  CB   CA
DIHE HD3  CD   CG1  CB    DIHE HD2  CD   CG1  CB      DIHE HD1  CD   CG1  CB
DONO HN   N
ACCE O    C
END
!======================================================================


RESI LEU
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.200 END
 ATOM HB1  TYPE=HA  CHARge=0.100  END
 ATOM HB2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CG   TYPE=CT  CHARge=-0.100 END
 ATOM HG   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CD1  TYPE=CT  CHARge=-0.300 END
 ATOM HD11 TYPE=HA  CHARge=0.100  END
 ATOM HD12 TYPE=HA  CHARge=0.100  END
 ATOM HD13 TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CD2  TYPE=CT  CHARge=-0.300 END
 ATOM HD21 TYPE=HA  CHARge=0.100  END
 ATOM HD22 TYPE=HA  CHARge=0.100  END
 ATOM HD23 TYPE=HA  CHARge=0.100  END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND CB CA     BOND CG CB    BOND CD1 CG    BOND CD2 CG 
BOND N HN      BOND N CA     BOND O C       BOND C CA 
BOND CA HA     BOND CB HB1   BOND CB HB2    BOND CG HG       BOND CD1 HD11
BOND CD1 HD12  BOND CD1 HD13 BOND CD2 HD21  BOND CD2 HD22    BOND CD2 HD23
IMPRoper CA   N    C    CB  !tetrahedral CA
DIHE CG   CB   CA   N     DIHE CD1  CG   CB   CA
DIHE HD13 CD1  CG   CB    DIHE HD12 CD1  CG   CB     DIHE HD11 CD1  CG   CB
DIHE HD23 CD2  CG   CB    DIHE HD22 CD2  CG   CB     DIHE HD21 CD2  CG   CB
DONO HN   N
ACCE O    C
END
!======================================================================


RESI LYS
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.200 END
 ATOM HB1  TYPE=HA  CHARge=0.100  END
 ATOM HB2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CG   TYPE=CT  CHARge=-0.200 END
 ATOM HG1  TYPE=HA  CHARge=0.100  END
 ATOM HG2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CD   TYPE=CT  CHARge=-0.200 END
 ATOM HD1  TYPE=HA  CHARge=0.100  END
 ATOM HD2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CE   TYPE=CT  CHARge=0.305  END
 ATOM HE1  TYPE=HA  CHARge=0.100  END
 ATOM HE2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM NZ   TYPE=NH3 CHARge=-0.810 END
 ATOM HZ1  TYPE=HC  CHARge=0.435  END
 ATOM HZ2  TYPE=HC  CHARge=0.435  END
 ATOM HZ3  TYPE=HC  CHARge=0.435  END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND CB CA    BOND CG CB    BOND CD CG    BOND CE CD    BOND NZ   CE
BOND N HN     BOND N CA     BOND O C      BOND C CA
BOND CA HA    BOND CB HB1   BOND CB HB2   BOND CG HG1
BOND CG HG2   BOND CD HD1   BOND CD HD2   BOND CE HE1   BOND CE HE2
BOND NZ HZ1   BOND NZ HZ2   BOND NZ HZ3
IMPRoper CA   N    C    CB  !tetrahedral CA
DIHE CG   CB   CA   N      DIHE CD   CG   CB   CA     DIHE CE   CD   CG   CB
DIHE NZ   CE   CD   CG
DIHE HZ3  NZ   CE   CD     DIHE HZ2  NZ   CE   CD     DIHE HZ1  NZ   CE   CD
DONO HN   N
DONO HZ1  NZ
DONO HZ2  NZ
DONO HZ3  NZ
ACCE O    C
END
!==========================================================================

RESI MET
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.200 END
 ATOM HB1  TYPE=HA  CHARge=0.100  END
 ATOM HB2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CG   TYPE=CT  CHARge=-0.115 END
 ATOM HG1  TYPE=HA  CHARge=0.100  END
 ATOM HG2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM SD   TYPE=S   CHARge=-0.170 END
GROUP 
 ATOM CE   TYPE=CT  CHARge=-0.215 END
 ATOM HE1  TYPE=HA  CHARge=0.100  END
 ATOM HE2  TYPE=HA  CHARge=0.100  END
 ATOM HE3  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND CB CA    BOND CG CB     BOND SD CG    BOND CE SD 
BOND N HN     BOND N CA      BOND O C      BOND C CA 
BOND CA HA    BOND CB HB1    BOND CB HB2   BOND CG HG1    BOND CG HG2
BOND CE HE1   BOND CE HE2    BOND CE HE3
IMPRoper CA   N    C    CB  !tetrahedral CA
DIHE CG   CB   CA   N     DIHE SD   CG   CB   CA      DIHE CE   SD   CG   CB
DIHE HE3  CE   SD   CG    DIHE HE2  CE   SD   CG      DIHE HE1  CE   SD   CG
DONO HN   N
ACCE O    C
END
!======================================================================


RESI PHE
GROUP
 ATOM N    TYPE=NH1    CHARge=-0.360 END
 ATOM HN   TYPE=H       CHARge=0.260 END
GROUP
 ATOM CA   TYPE=CT      CHARge=0.000 END
 ATOM HA   TYPE=HA      CHARge=0.100 END
GROUP
 ATOM CB   TYPE=CT     CHARge=-0.160 END
 ATOM HB1  TYPE=HA      CHARge=0.100 END
 ATOM HB2  TYPE=HA      CHARge=0.100 END
GROUP
 ATOM CG   TYPE=C       CHARge=0.030 EXCL=( CZ )  END
GROUP
 ATOM CD1  TYPE=C      CHARge=-0.160 EXCL=( CE2 ) END
 ATOM HD1  TYPE=HA      CHARge=0.140 END
GROUP
 ATOM CD2  TYPE=C      CHARge=-0.160 EXCL=( CE1 ) END
 ATOM HD2  TYPE=HA      CHARge=0.140 END
GROUP
 ATOM CE1  TYPE=C      CHARge=-0.150 END
 ATOM HE1  TYPE=HA      CHARge=0.140 END
GROUP
 ATOM CE2  TYPE=C      CHARge=-0.150 END
 ATOM HE2  TYPE=HA      CHARge=0.140 END
GROUP
 ATOM CZ   TYPE=C      CHARge=-0.150 END
 ATOM HZ   TYPE=HA      CHARge=0.140 END
GROUP
 ATOM C    TYPE=C       CHARge=0.480 END
 ATOM O    TYPE=O      CHARge=-0.480 END
BOND CB CA     BOND CG CB    BOND CD1 CG    BOND CD2 CG    BOND CE1 CD1
BOND CE2 CD2   BOND CZ CE1   BOND CZ CE2    BOND N HN
BOND N CA      BOND O C      BOND C CA      BOND CA HA
BOND CB HB1    BOND CB HB2   BOND CD1 HD1   BOND CD2 HD2   BOND CE1 HE1
BOND CE2 HE2   BOND CZ HZ
DIHE CG   CB   CA   N         DIHE CD1  CG   CB   CA
IMPRoper CA   N    C    CB  !tetrahedral CA
IMPR CE1  CD1  CG   CD2    IMPR CE2  CD2  CG   CD1     IMPR CZ   CE1  CD1  CG
IMPR CZ   CE2  CD2  CG     IMPR CE2  CZ   CE1  CD1     
IMPR CG   CB   CD1  CD2    IMPR CD1  CG   CE1  HD1     IMPR CD2  CG   CE2  HD2
IMPR CE1  CZ   CE2  CD2    IMPR CE1  CD1  CZ   HE1     IMPR CE2  CD2  CZ   HE2
IMPR CZ   CE1  CE2  HZ 
DONO HN   N
ACCE O    C
END
!======================================================================


RESI PRO
GROUP
 ATOM N    TYPE=N   CHARge=-0.360 END
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.200 END
 ATOM HB1  TYPE=HA  CHARge=0.100  END
 ATOM HB2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CG   TYPE=CT  CHARge=-0.200 END
 ATOM HG1  TYPE=HA  CHARge=0.100  END
 ATOM HG2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CD   TYPE=CT  CHARge=0.060  END
 ATOM HD1  TYPE=HA  CHARge=0.100  END
 ATOM HD2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND HA CA    BOND CB CA     BOND HB1 CB    BOND HB2 CB    BOND CG CB
BOND HG1 CG   BOND HG2 CG    BOND CD CG     BOND HD1 CD    BOND HD2 CD
BOND N CA     BOND O C       BOND C CA      BOND N CD
IMPRoper CA   N    C    CB  !tetrahedral CA
DIHE N    CA   CB   CG     DIHE CA   CB   CG   CD     DIHE N    CD   CG   CB
ACCE O    C
END
!======================================================================


RESI SER
GROUP
 ATOM N    TYPE=NH1  CHARge=-0.360 END
 ATOM HN   TYPE=H    CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT   CHARge=0.000  END
 ATOM HA   TYPE=HA   CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT   CHARge=0.080  END
 ATOM HB1  TYPE=HA   CHARge=0.100  END
 ATOM HB2  TYPE=HA   CHARge=0.100  END
GROUP
 ATOM OG   TYPE=OH1  CHARge=-0.680 END
 ATOM HG1  TYPE=H    CHARge=0.400  END
GROUP
 ATOM C    TYPE=C    CHARge=0.480  END
 ATOM O    TYPE=O    CHARge=-0.480 END
BOND CB CA      BOND OG CB        BOND N HN     BOND N CA
BOND O C        BOND C CA         BOND CA HA    BOND CB HB1
BOND CB HB2     BOND OG HG1
DIHE OG   CB   CA   N         DIHE HG1  OG   CB   CA
IMPRoper CA   N    C    CB  !tetrahedral CA
DONO HN   N
DONO HG1  OG
ACCE OG  " "
ACCE O    C
END
!======================================================================


RESI THR
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=0.180  END
 ATOM HB   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM OG1  TYPE=OH1 CHARge=-0.680 END
 ATOM HG1  TYPE=H   CHARge=0.400  END
GROUP
 ATOM CG2  TYPE=CT  CHARge=-0.300 END
 ATOM HG21 TYPE=HA  CHARge=0.100  END
 ATOM HG22 TYPE=HA  CHARge=0.100  END
 ATOM HG23 TYPE=HA  CHARge=0.100  END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND CB CA    BOND OG1 CB    BOND CG2 CB     BOND N HN
BOND N CA     BOND O C       BOND C CA       BOND CA HA
BOND CB HB    BOND OG1 HG1   BOND CG2 HG21   BOND CG2 HG22    BOND CG2 HG23
DIHE OG1  CB   CA   N     DIHE HG1  OG1  CB   CA
DIHE HG23 CG2  CB   CA    DIHE HG22 CG2  CB   CA       DIHE HG21 CG2  CB   CA
IMPRoper CA   N    C    CB  !tetrahedral CA
DONO HN   N
DONO HG1  OG1
ACCE OG1 " "
ACCE O    C
END
!======================================================================


RESI TYR
GROUP
 ATOM N    TYPE=NH1    CHARge=-0.360 END
 ATOM HN   TYPE=H       CHARge=0.260 END
GROUP
 ATOM CA   TYPE=CT      CHARge=0.000 END
 ATOM HA   TYPE=HA      CHARge=0.100 END
GROUP
 ATOM CB   TYPE=CT     CHARge=-0.200 END
 ATOM HB1  TYPE=HA      CHARge=0.100 END
 ATOM HB2  TYPE=HA      CHARge=0.100 END
GROUP
 ATOM CG   TYPE=C       CHARge=0.000 EXCL=( CZ ) END
GROUP
 ATOM CD1  TYPE=C      CHARge=-0.140 EXCL=( CE2 ) END
 ATOM HD1  TYPE=HA      CHARge=0.140 END
GROUP
 ATOM CD2  TYPE=C      CHARge=-0.140 EXCL=( CE1 ) END
 ATOM HD2  TYPE=HA      CHARge=0.140 END
GROUP
 ATOM CE1  TYPE=C      CHARge=-0.140 END
 ATOM HE1  TYPE=HA      CHARge=0.140 END
GROUP
 ATOM CE2  TYPE=C      CHARge=-0.140 END
 ATOM HE2  TYPE=HA      CHARge=0.140 END
GROUP
 ATOM CZ   TYPE=C       CHARge=0.200 END
 ATOM OH   TYPE=OH1    CHARge=-0.600 END
 ATOM HH   TYPE=H       CHARge=0.400 END
GROUP
 ATOM C    TYPE=C       CHARge=0.480 END
 ATOM O    TYPE=O      CHARge=-0.480 END
BOND CB CA    BOND CG CB   BOND CD1 CG    BOND CD2 CG    BOND CE1 CD1
BOND CE2 CD2  BOND CZ CE1  BOND CZ CE2    BOND OH CZ        
BOND N HN     BOND N CA    BOND O C       BOND C CA      
BOND CA HA    BOND CB HB1  BOND CB HB2    BOND CD1 HD1   BOND CD2 HD2
BOND CE1 HE1  BOND CE2 HE2 BOND OH HH
DIHE CG   CB   CA   N     DIHE CD1  CG   CB   CA    DIHE HH   OH   CZ   CE1
IMPR CE1  CD1  CG   CD2   IMPR CE2  CD2  CG   CD1   IMPR CZ   CE1  CD1  CG
IMPR CZ   CE2  CD2  CG    IMPR CE2  CZ   CE1  CD1   IMPR CE1  CZ   CE2  CD2
IMPR CG   CD1  CD2  CB    IMPR CD1  CG   CE1  HD1
IMPR CD2  CG   CE2  HD2   IMPR CE1  CD1  CZ   HE1   IMPR CE2  CD2  CZ   HE2
IMPR CZ   CE1  CE2  OH    
IMPRoper CA   N    C    CB  !tetrahedral CA
DONO HN   N
DONO HH   OH
ACCE OH  " "
ACCE O    C
END
!======================================================================


RESI VAL
GROUP
 ATOM N    TYPE=NH1  CHARge=-0.360 END
 ATOM HN   TYPE=H    CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT   CHARge=0.000  END
 ATOM HA   TYPE=HA   CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT   CHARge=-0.100 END
 ATOM HB   TYPE=HA   CHARge=0.100  END
GROUP
 ATOM CG1  TYPE=CT   CHARge=-0.300 END
 ATOM HG11 TYPE=HA   CHARge=0.100  END
 ATOM HG12 TYPE=HA   CHARge=0.100  END
 ATOM HG13 TYPE=HA   CHARge=0.100  END
GROUP
 ATOM CG2  TYPE=CT   CHARge=-0.300 END
 ATOM HG21 TYPE=HA   CHARge=0.100  END
 ATOM HG22 TYPE=HA   CHARge=0.100  END
 ATOM HG23 TYPE=HA   CHARge=0.100  END
GROUP
 ATOM C    TYPE=C    CHARge=0.480  END
 ATOM O    TYPE=O    CHARge=-0.480 END
BOND CB CA     BOND CG1 CB     BOND CG2 CB       BOND N HN
BOND N CA      BOND O C        BOND C CA         BOND CA HA
BOND CB HB     BOND CG1 HG11   BOND CG1 HG12     BOND CG1 HG13   BOND CG2 HG21
BOND CG2 HG22  BOND CG2 HG23
IMPRoper CA   N    C    CB  !tetrahedral CA
DIHE CG1  CB   CA   N
DIHE HG13 CG1  CB   CA    DIHE HG12 CG1  CB   CA     DIHE HG11 CG1  CB   CA
DIHE HG23 CG2  CB   CA    DIHE HG22 CG2  CB   CA     DIHE HG21 CG2  CB   CA
DONO HN   N
ACCE O    C
END
!======================================================================


RESI TRP
GROUP
 ATOM N    TYPE=NH1 CHARge=-0.360 END
 ATOM HN   TYPE=H   CHARge=0.260  END
GROUP
 ATOM CA   TYPE=CT  CHARge=0.000  END
 ATOM HA   TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CB   TYPE=CT  CHARge=-0.200 END
 ATOM HB1  TYPE=HA  CHARge=0.100  END
 ATOM HB2  TYPE=HA  CHARge=0.100  END
GROUP
 ATOM CG   TYPE=C   CHARge=-0.040 EXCL=( CD1  CD2  NE1  CE2  CE3 ) END
GROUP
 ATOM CD1  TYPE=C   CHARge=-0.010 EXCL=( CD2  NE1  CE2 ) END
 ATOM HD1  TYPE=HA  CHARge=0.140  END
GROUP
 ATOM CD2  TYPE=C   CHARge=-0.050 EXCL=( NE1  CE2  CE3  CZ2  CZ3 ) END
GROUP
 ATOM NE1  TYPE=NH1 CHARge=-0.490 EXCL=( CE2  CZ2 ) END
 ATOM HE1  TYPE=H   CHARge=0.400  END
GROUP
 ATOM CE2  TYPE=C   CHARge=0.130  EXCL=( CE3  CZ2  CH2 ) END
GROUP
 ATOM CE3  TYPE=C   CHARge=-0.160 EXCL=( CZ3  CH2 ) END
 ATOM HE3  TYPE=HA  CHARge=0.140  END
GROUP
 ATOM CZ2  TYPE=C   CHARge=-0.150 EXCL=( CZ3  CH2 ) END
 ATOM HZ2  TYPE=HA  CHARge=0.140  END
GROUP
 ATOM CZ3  TYPE=C   CHARge=-0.150 EXCL=( CH2 ) END
 ATOM HZ3  TYPE=HA  CHARge=0.140  END
GROUP
 ATOM CH2  TYPE=C   CHARge=-0.180 END
 ATOM HH2  TYPE=HA  CHARge=0.140  END
GROUP
 ATOM C    TYPE=C   CHARge=0.480  END
 ATOM O    TYPE=O   CHARge=-0.480 END
BOND CB CA     BOND CG CB    BOND CD1 CG    BOND CD2 CG    BOND NE1 CD1
BOND CE2 CD2   BOND CZ2 CE2  BOND CZ3 CE3   BOND CH2 CZ2       
BOND N HN      BOND N CA     BOND O C       BOND C CA  
BOND CZ3 CH2   BOND CD2 CE3  BOND NE1 CE2   BOND CA HA     BOND CB HB1
BOND CB HB2    BOND CD1 HD1  BOND NE1 HE1   BOND CE3 HE3   BOND CZ2 HZ2
BOND CZ3 HZ3   BOND CH2 HH2
OMIT ANGLe CG CD2 CE3
OMIT ANGLe NE1 CE2 CZ2
DIHE CG   CB   CA   N     DIHE CD1  CG   CB   CA
IMPRoper CA   N    C    CB  !tetrahedral CA
IMPR NE1  CD1  CG   CD2   IMPR CE2  CD2  CG   CD1   IMPR CE2  NE1  CD1  CG
IMPR NE1  CE2  CD2  CG    IMPR CH2  CZ2  CE2  CD2   IMPR CH2  CZ3  CE3  CD2
IMPR CZ3  CH2  CZ2  CE2   IMPR CE3  CZ3  CH2  CZ2   IMPR CE2  CD2  CE3  CZ3
IMPR CD1  NE1  CE2  CD2   IMPR CG   CD1  CD2  CB
IMPR CD1  CG   NE1  HD1   IMPR CD2  CG   CE2  CE3   IMPR NE1  CD1  CE2  HE1
IMPR CE2  CD2  CZ2  NE1   IMPR CE3  CZ3  CD2  HE3   IMPR CZ2  CE2  CH2  HZ2
IMPR CZ3  CE3  CH2  HZ3   IMPR CH2  CZ2  CZ3  HH2   
DONO HN   N
DONO HE1  NE1
ACCE O    C
END
!======================================================================


RESI FMA     ! FORMAMIDE WITH CHARGES FROM EXPERIMENTAL DIPOLE
GROUP
 ATOM HA   TYPE=HA   CHARge=0.000  END
 ATOM C    TYPE=C    CHARge=0.480  END
GROUP
 ATOM O    TYPE=O    CHARge=-0.480 END
 ATOM N    TYPE=NH2  CHARge=-0.520 END
 ATOM HT   TYPE=H    CHARge=0.260  END
 ATOM HC   TYPE=H    CHARge=0.260  END
BOND HA   C      BOND C O      BOND C  N    BOND N HT   BOND N    HC
DIHE O    C    N    HC        DIHE O    C    N    HT
IMPR C    N    O    HA        IMPR N    HC   HT   C
ACCE O    C
DONO HC   N
DONO HT   N
END
!======================================================================

!------------------------------------------------------------------

PRESidue NTER             { can be patched ( as NTER - * ... )
                            to any amino acid except PRO        }
  MODIfy ATOM +CA             CHARge=0.22  END
 GROUP
  MODIfy ATOM +N    TYPE=NH3  CHARge=-0.10 END
  DELETE ATOM +HN                          END
  ADD    ATOM +HT1  TYPE=HC   CHARge=0.26  END
  ADD    ATOM +HT2  TYPE=HC   CHARge=0.26  END
  ADD    ATOM +HT3  TYPE=HC   CHARge=0.26  END

 ADD BOND +HT1 +N
 ADD BOND +HT2 +N
 ADD BOND +HT3 +N

 ADD ANGLe +HT1  +N    +HT2
 ADD ANGLe +HT2  +N    +HT3
 ADD ANGLe +HT2  +N    +CA
 ADD ANGLe +HT1  +N    +HT3
 ADD ANGLe +HT1  +N    +CA
 ADD ANGLe +HT3  +N    +CA

 ADD DIHEdral +HT2  +N    +CA   +C
 ADD DIHEdral +HT1  +N    +CA   +C
 ADD DIHEdral +HT3  +N    +CA   +C

 ADD DONOr +HT1  +N
 ADD DONOr +HT2  +N
 ADD DONOr +HT3  +N

END {NTER}


PRESidue PROP                { this is the N-terminal for PROlines
                               PROP - PRO -...                      }
  MODIfy ATOM +CD             CHARge= 0.10   END
  MODIfy ATOM +CA             CHARge=0.10    END
 GROUp
  MODIfy ATOM +N    TYPE=NH3  CHARge=-0.02   END
  ADD    ATOM +HT1  TYPE=HC   CHARge= 0.26   END
  ADD    ATOM +HT2  TYPE=HC   CHARge= 0.26   END
!!  MODIfy ATOM +HA             CHARge=0.10  END

 ADD BOND +HT1  +N
 ADD BOND +HT2  +N

 ADD ANGLe +HT1  +N    +HT2
 ADD ANGLe +HT2  +N    +CA
 ADD ANGLe +HT1  +N    +CD
 ADD ANGLe +HT1  +N    +CA
 ADD ANGLe +CD   +N    +HT2

 ADD DIHEdral +HT2  +N    +CA   +C
 ADD DIHEdral +HT1  +N    +CA   +C

 ADD DONOr +HT1  +N
 ADD DONOr +HT2  +N

END {PROP}

!------------------------------------------------------------------

PRESidue CTER                { C-terminal for all amino acids
                                           ... * - CTER          }
 GROUp
  MODIfy ATOM -C             CHARge= 0.14  END
  ADD    ATOM -OT1  TYPE=OC  CHARge=-0.57  END
  ADD    ATOM -OT2  TYPE=OC  CHARge=-0.57  END
  DELETE ATOM -O                           END

 ADD BOND -C    -OT1
 ADD BOND -C    -OT2

 ADD ANGLe -CA   -C   -OT1
 ADD ANGLe -CA   -C   -OT2
 ADD ANGLe -OT1  -C   -OT2

 ADD DIHEdral -N    -CA    -C   -OT2

 ADD IMPRoper -C    -CA    -OT2 -OT1

 ADD ACCEptor -OT1 -C

 ADD ACCEptor -OT2 -C

END {CTER}


!-------------------------------------------------------------------
PRESidue CTN                { C-terminal for all amino acids
                                           ... * - CTN CONH2 at end         }
 GROUp
  MODIfy ATOM -C             CHARge= 0.48  END
 GROUP
  MODIFY ATOM -O             CHARge=-0.48  END
  ADD    ATOM -NT   TYPE=NH2 CHARge=-0.52  END
  ADD    ATOM  -H1   TYPE=H   CHARGE= 0.26  END
  ADD    ATOM  -H2   TYPE=H   CHARGE= 0.26  END
     
 ADD BOND -C    -NT
 ADD BOND -NT    -H1
 ADD BOND -NT    -H2
     
 ADD ANGLe -CA   -C   -NT
 ADD ANGLe -O    -C  -NT
 ADD ANGLe -CA   -C  -O     !added MN 23-jun-1988
 ADD ANGLe -C    -NT   -H1   !  "
 ADD ANGLe -C    -NT   -H2   !  "
 ADD ANGLe  -H1   -NT   -H2   !  "
     
 ADD DIHEdral -N   -CA    -C   -NT
 add DIHEdral -O    -C     -NT   -H1
 add DIHEdral -O    -C     -NT  -H2
     
 ADD IMPRoper -C    -CA    -NT  -O
 ADD IMPRoper -C    -NT    -O   -H1
 ADD IMPRoper -NT    -H1     -H2  -C
     
 ADD DONOR  -H1 -NT
 ADD DONOR  -H2 -NT
     
END {CTN}

!----------------------------------------------------------------------
     

PRESidue PEPT { PEPTide bond link, for all 
               amino acids ...*(-)     (+)*...
                                \ PEPT /

               except the  *(-) - (+)PRO link        }

 ADD BOND -C +N 

 ADD ANGLE -CA -C +N
 ADD ANGLE -O  -C +N
 ADD ANGLE -C  +N +CA
 ADD ANGLE -C  +N +HN

 ADD DIHEdral  -C +N +CA +C
 ADD DIHEdral  -N -CA -C +N
 ADD DIHEdral  -CA -C +N +CA

 ADD IMPRoper  -C -CA +N -O  {planar -C}
 ADD IMPRoper  +N -C +CA +HN  {planar +N}

END {PEPT}

!----------------------------------------------------------------------
PRESidue PEPP  { for  ...*(-) - (+)PRO  link
               same as PEPT except replacement H by CD
               and improper +N +CA +CD -C              }

 ADD BOND -C +N 

 ADD ANGLE -CA -C +N
 ADD ANGLE -O  -C +N
 ADD ANGLE -C  +N +CA
 ADD ANGLE -C  +N +CD

 ADD DIHEdral  -C +N +CD +CG    !!! ADDED
 ADD DIHEdral  -C +N +CA +C
 ADD DIHEdral  -N -CA -C +N
 ADD DIHEdral  -CA -C +N +CA

 ADD IMPRoper  -C -CA +N -O  {planar -C}
 ADD IMPRoper  +N -C +CA +CD  !!! +N  +CA +CD -C  {planar +N} MODIFIED

END {PEPP}

!------------------------------------------------------------------


RESI OH2  
GROUP
 ATOM OH2  TYPE=OH2  CHARGE=-0.50  END
 ATOM H1   TYPE=H    CHARGE=0.25   END
 ATOM H2   TYPE=H    CHARGE=0.25   END
BOND OH2  H1        BOND OH2  H2
DONO H1   OH2
DONO H2   OH2
ACCE OH2  " "
END

!------------------------------------------------------------------

PRESidue DISU       { disulfide bridge  ...CYS      CYS...
                                              \DISU/            }
DELETE ATOM 1HG1 END
DELETE ATOM 2HG1 END
 GROUP
  MODIfy ATOM 1CB           CHARge=-0.20  END
  MODIfy ATOM 1SG  TYPE=S   CHARge=-0.20  END
 GROUP
  MODIfy ATOM 2CB           CHARge=-0.20  END
  MODIfy ATOM 2SG  TYPE=S   CHARge=-0.20  END

 ADD BOND 1SG 2SG

 ADD ANGLe  1CB 1SG 2SG
 ADD ANGLe  1SG 2SG 2CB

 ADD DIHEdral   1CA 1CB 1SG 2SG
 ADD DIHEdral   1CB 1SG 2SG 2CB
 ADD DIHEdral   1SG 2SG 2CB 2CA

END {DISU}
!------------------------------------------------------------------------

set echo=true end
