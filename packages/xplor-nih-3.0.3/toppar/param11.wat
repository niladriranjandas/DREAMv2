remarks   PARAM11.WAT   (water parameters)
remarks   ===========
remarks   available: TIPS3P model
set echo=false end


{* TIP3P model *}
{* =========== *}
BOND HT   OT     450.0       0.9572    ! from TIPS3P geometry ( SHAKE w/PARAm)
BOND HT   HT       0.0       1.5139    ! from TIPS3P geometry ( SHAKE w/PARAm)
ANGLE HT   OT   HT      55.0     104.52    ! FROM TIPS3P geometry
NONBONDED OT        0.1521   3.1506       0.1521  3.1506  !TIPS3P water oxygen
NONBONDED HT        0.04598  0.4000       0.04598 0.4000  !TIPS3P water hydr.

set echo=true end
