remarks  TOPH11.WAT
remarks  ==========
remarks  topology file for water
remarks  available: TIP3P model

set echo=false end

{* default masses *}
MASS   HT     1.00800! TIPS3P water hydrogen
MASS   OT    15.99940 ! TIPS3P water oxygen

AUTOGENERATE ANGLES=FALSE END

!------------------------------------------------------------------

RESIdue TIP3       { TIPS3P WATER MODEL }
 GROUp
  ATOM OH2  TYPE=OT   CHARge= -0.834  END
  ATOM H1   TYPE=HT   CHARge=  0.417  END
  ATOM H2   TYPE=HT   CHARge=  0.417  END

 BOND OH2  H1
 BOND OH2  H2
 ANGLE  H1 OH2 H2

END {TIP3P}

!------------------------------------------------------------------
set echo=true end
