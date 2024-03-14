remark   para_axis_3.pro
remark    NOTE: these parameters are DEPRECATED. Please instead use axis.par.
remark molecule for anisotropy

set message off echo off end


! BONDS

 bonds     OOO  XXX                  1000.0     3.00
 bonds     OOO  YYY                  1000.0     3.00
 bonds     OOO  ZZZ                  1000.0     3.00

! ANGLES
 angle     XXX OOO YYY             1000.0      90.0 
angle     XXX OOO ZZZ             1000.0      90.0
angle     YYY OOO ZZZ             1000.0      90.0


 NONBonded  XXX      0.145    3.2072      0.145    3.2072
NONBonded  YYY      0.145    3.2072      0.145    3.2072
NONBonded  ZZZ      0.145    3.2072      0.145    3.2072
 NONBonded  OOO      0.145    3.2072      0.145    3.2072

VDWOff OOO
VDWOff XXX
VDWOff YYY
VDWOff ZZZ

set message on echo on end
