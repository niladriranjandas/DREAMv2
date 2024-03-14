MASS OOO 1.0
MASS XXX 1.0
MASS YYY 1.0
MASS ZZZ 1.0

autogenerate angles=true end

residue ANI
  group
    atom X  type=XXX  charge=-0.0 end
    atom Y  type=YYY  charge= 0.0 end
    atom Z  type=ZZZ  charge= 0.0 end
    atom OO type=OOO  charge= 0.0 end
  angle X OO Y angle X OO Z angle Y OO Z
  bond OO X bond OO Y  bond OO Z

!  improper X Y OO Z  !stereo CB

end
