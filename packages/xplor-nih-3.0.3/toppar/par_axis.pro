remark   para_axis.pro 
remark    NOTE: these parameters are DEPRECATED. Please instead use axis.par.
remark molecule for anisotropy

set message off echo off end

evaluate ($kbon = 1000)  ! kcal / mol-A^2
evaluate ($kang =  500)  ! kcal / mol-rad^2
evaluate ($kchi =  500)  ! kcal / mol-rad^2
evaluate ($kback = 500)
evaluate ($kssbon = 1000)
evaluate ($kssang = 500) 
evaluate ($kpla =  500)  ! kcal / mol-rad^2
evaluate ($kdih =    0)  ! kcal / mol-rad^2


! BONDS

 bonds     OOO  XXX                  $kbon     1.00
 bonds     OOO  YYY                  $kbon     1.00
 bonds     OOO  ZZZ                  $kbon     1.00

! ANGLES
 angle     XXX OOO YYY             $kang      90.0 
angle     XXX OOO ZZZ             $kang      90.0
angle     YYY OOO ZZZ             $kang      90.0


 NONBonded  XXX      0.145    3.2072      0.145    3.2072
NONBonded  YYY      0.145    3.2072      0.145    3.2072
NONBonded  ZZZ      0.145    3.2072      0.145    3.2072
 NONBonded  OOO      0.145    3.2072      0.145    3.2072

VDWOff OOO
VDWOff XXX
VDWOff YYY
VDWOff ZZZ

set message on echo on end