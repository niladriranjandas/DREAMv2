REMARKS  TOPH22.pep -MACRO for protein sequence
SET ECHO=FALSE END

! this is a macro to define standard protein peptide bonds
! and termini to generate a protein sequence.
! 
! Author: Axel Brunger, 19-JAN-84
! To be used with the Charmm22 topology/parameter sets - PDA 19-JUL-93


LINK PEPP    HEAD - *     TAIL + PRO     END  { LINK to PRO }
LINK PPG1    HEAD - *     TAIL + GLY     END  { LINK to GLY }
LINK PPG2    HEAD - GLY   TAIL + *       END  { LINK from GLY }
LINK PPGP    HEAD - GLY   TAIL + PRO     END  { LINK from GLY to PRO }
LINK PPGG    HEAD - GLY   TAIL + GLY     END  { LINK from GLY to GLY }
LINK PEPT    HEAD - *     TAIL + *       END

FIRSt PROP                TAIL + PRO     END { Nter for PRO }
FIRSt GLYP                TAIL + GLY     END { Nter for GLY }
FIRSt NTER                TAIL + *       END

LAST  CTER   HEAD - *                    END

SET ECHO=TRUE END

