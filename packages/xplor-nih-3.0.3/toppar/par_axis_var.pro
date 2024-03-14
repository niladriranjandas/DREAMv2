remark   para_anis.pro
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
 bonds     OOO  PAR                  $kbon     1.00


! ANGLES
 angle     XXX OOO YYY             $kang      90.0
angle     XXX OOO ZZZ             $kang      90.0
angle     YYY OOO ZZZ             $kang      90.0


 NONBonded  XXX      0.145    0.001      0.145    0.001
NONBonded  YYY      0.145    0.001      0.145    0.001
NONBonded  ZZZ      0.145    0.001      0.145    0.001
 NONBonded  OOO      0.145    0.001      0.145    0.001
 NONBonded  PAR      0.145    0.001      0.145    0.001

VDWOff OOO
VDWOff XXX
VDWOff YYY
VDWOff ZZZ
VDWOff PAR

set message=$prev_messages echo=$prev_echo end

