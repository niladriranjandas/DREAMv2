remark  file topallhdg.pro
remark   geometric energy function topology for distance geometry and
remark   simulated annealing.
remark 
remark   history:
remark
remark     XPLOR topology file TOPALLHSA.PRO, all hydrogens, 
remark     for simulated annealing, to be used with PARALLHSA.PRO
remark     based on TOPALLH6.PRO.
remark     internal coordinates from TOPH19.PRO.

! all dihedrals defining planarity have been converted to impropers,
! the only dihedrals left are around rotatable bonds.
! some additional impropers were included in the rings to define planarity
! better. patches are included to define chirality at stereospecifically 
! assigned prochiral centres and convert from trans to cis peptide bond.
! modifications by M. Nilges and P. Kraulis.
! added stereospecific impropers for all pro-chiral centers (ATB)


set message off echo off end

autogenerate 
  angles=true
  dihedrals=false
end

mass H    1.008
mass HC   1.008
mass HA   1.008
mass C   12.011
mass CT  12.011
mass CX  12.011
mass CW  12.011
mass CN  12.011
mass CA  12.011
mass CB  12.011
mass CC  12.011
mass CV  12.011
mass CR  12.011
mass CP  12.011
mass N   14.007
mass NA  14.007
mass NB  14.007
mass NH1 14.007
mass NH2 14.007
mass NH3 14.007
mass NC2 14.007
mass O   15.999
mass OC  15.999
mass OH  15.999
mass S   32.060


residue ALA
  group
    atom N   type=NH1 charge=-0.36 end
    atom HN  type=H   charge= 0.26 end
    atom CA  type=CT  charge= 0.00 end
    atom HA  type=HA  charge= 0.10 end
    atom CB  type=CT  charge=-0.30 end
    atom HB1 type=HA  charge= 0.10 end
    atom HB2 type=HA  charge= 0.10 end
    atom HB3 type=HA  charge= 0.10 end
    atom C   type=C   charge= 0.48 end
    atom O   type=O   charge=-0.48 end

  bond N  HN       
  bond N  CA    bond CA  HA   
  bond CA  CB   bond CB  HB1   bond CB  HB2     bond CB  HB3
  bond CA  C       
  bond C   O

  improper HA  N   C   CB  !stereo CA
  improper HB1 HB2 CA HB3  !stereo CB

  ic N C *CA CB  0.0 0.0  120.0  0.0 0.0
end
 

residue ARG
  group
    atom N    type=NH1  charge=-0.36 end
    atom HN   type=H    charge= 0.26 end
    atom CA   type=CT   charge= 0.00 end
    atom HA   type=HA   charge= 0.10 end
    atom CB   type=CT   charge=-0.20 end
    atom HB1  type=HA   charge= 0.10 end
    atom HB2  type=HA   charge= 0.10 end
    atom CG   type=CT   charge=-0.20 end
    atom HG1  type=HA   charge= 0.10 end
    atom HG2  type=HA   charge= 0.10 end
    atom CD   type=CT   charge=-0.10 end
    atom HD1  type=HA   charge= 0.15 end
    atom HD2  type=HA   charge= 0.15 end
    atom NE   type=NH1  charge=-0.60 end
    atom HE   type=H    charge= 0.40 end
    atom CZ   type=C    charge= 0.60 end
    atom NH1  type=NC2  charge=-0.60 end
    atom HH11 type=HC   charge= 0.40 end
    atom HH12 type=HC   charge= 0.40 end
    atom NH2   type=NC2 charge=-0.60 end
    atom HH21  type=HC  charge= 0.40 end
    atom HH22  type=HC  charge= 0.40 end
    atom C     type=C   charge= 0.48 end
    atom O     type=O   charge=-0.48 end
 
  bond N  HN
  bond N  CA     bond CA  HA
  bond CA CB     bond CB  HB1     bond CB  HB2
  bond CB CG     bond CG  HG1     bond CG  HG2
  bond CG CD     bond CD  HD1     bond CD  HD2
  bond CD NE     bond NE  HE
  bond NE CZ
  bond CZ NH1    bond NH1 HH11    bond NH1 HH12
  bond CZ NH2    bond NH2 HH21    bond NH2 HH22
  bond CA C
  bond C  O

  improper HA  N  C    CB  !stereo CA
  improper NE  CD CZ   HE  
  improper CZ  NE NH1  NH2
  improper NH1 CZ HH11 HH12
  improper NH2 CZ HH21 HH22
  improper NE  CZ NH1  HH11
  improper NE  CZ NH2  HH21
!  improper CZ  NH1 HE NE  ! planar HE, CZ
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB CD  !stereo CG
  improper HD1 HD2 CG NE  !stereo CD
  
  dihedral CG  CB  CA  N
  dihedral CD  CG  CB  CA
  dihedral NE  CD  CG  CB
  dihedral CZ  NE  CD  CG
  
  ic N   C  *CA CB   0.0 0.0  120.0  0.0 0.0
  ic N   CA CB  CG   0.0 0.0  180.0  0.0 0.0
  ic CA  CB CG  CD   0.0 0.0  180.0  0.0 0.0
  ic CB  CG CD  NE   0.0 0.0  180.0  0.0 0.0
  ic CG  CD NE  CZ   0.0 0.0  180.0  0.0 0.0
  ic CD  NE CZ  NH1  0.0 0.0  180.0  0.0 0.0
  ic NH1 NE *CZ NH2  0.0 0.0  180.0  0.0 0.0
  ic CD  CZ *NE HE   0.0 0.0  180.0  0.0 0.0
end

residue ASN
  group
    atom N    type=NH1 charge=-0.36 end
    atom HN   type=H   charge= 0.26 end
    atom CA   type=CT  charge= 0.00 end
    atom HA   type=HA  charge= 0.10 end
    atom CB   type=CT  charge=-0.20 end
    atom HB1  type=HA  charge= 0.10 end
    atom HB2  type=HA  charge= 0.10 end
    atom CG   type=C   charge= 0.48 end
    atom OD1  type=O   charge=-0.48 end
    atom ND2  type=NH2 charge=-0.52 end
    atom HD21 type=H   charge= 0.26 end
    atom HD22 type=H   charge= 0.26 end
    atom C    type=C   charge= 0.48 end
    atom O    type=O   charge=-0.48 end

  bond N  HN
  bond N  CA     bond CA  HA
  bond CA CB     bond CB  HB1     bond CB  HB2
  bond CB CG
  bond CG OD1
  bond CG ND2    bond ND2 HD21    bond ND2 HD22
  bond CA C
  bond C  O

  improper HA  N  C    CB  !stereo CA
  improper CG  CB OD1  ND2
  improper ND2 CG HD21 HD22
  improper CB  CG ND2  HD21
  improper HB1 HB2 CA CG  !stereo CB

  dihedral CG  CB  CA  N
  dihedral OD1 CG  CB  CA

  ic N   C  *CA CB   0.0 0.0  120.0  0.0 0.0
  ic N   CA CB  CG   0.0 0.0  180.0  0.0 0.0
  ic CA  CB CG  OD1  0.0 0.0  180.0  0.0 0.0
  ic OD1 CB *CG ND2  0.0 0.0  180.0  0.0 0.0
end


residue ASP
  group
    atom N   type=NH1 charge=-0.36 end
    atom HN  type=H   charge= 0.26 end
    atom CA  type=CT  charge= 0.00 end
    atom HA  type=HA  charge= 0.10 end
    atom CB  type=CT  charge=-0.45 end
    atom HB1 type=HA  charge= 0.10 end
    atom HB2 type=HA  charge= 0.10 end
    atom CG  type=C   charge= 0.49 end
    atom OD1 type=OC  charge=-0.62 end
    atom OD2 type=OC  charge=-0.62 end
    atom C   type=C   charge= 0.48 end
    atom O   type=O   charge=-0.48 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG
  bond CG OD1
  bond CG OD2
  bond CA C
  bond C  O

  improper HA N  C   CB   !stereo CA
  improper CG CB OD1 OD2
  improper HB1 HB2 CA CG  !stereo CB

  dihedral CG  CB  CA  N
  dihedral OD1 CG  CB  CA

  ic N   C  *CA CB   0.0 0.0  120.0  0.0 0.0
  ic N   CA CB  CG   0.0 0.0  180.0  0.0 0.0
  ic CA  CB CG  OD1  0.0 0.0  180.0  0.0 0.0
  ic OD1 CB *CG OD2  0.0 0.0  180.0  0.0 0.0
end


residue CYS
  group
    atom N   type=NH1 charge=-0.36 end
    atom HN  type=H   charge= 0.26 end
    atom CA  type=CT  charge= 0.00 end
    atom HA  type=HA  charge= 0.10 end
    atom CB  type=CT  charge=-0.20 end
    atom HB1 type=HA  charge= 0.10 end
    atom HB2 type=HA  charge= 0.10 end
    atom SG  type=S   charge=-0.05 end
    atom HG  type=H   charge= 0.05 end
    atom C   type=C   charge= 0.48 end
    atom O   type=O   charge=-0.48 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB SG     bond SG HG 
  bond CA C
  bond C  O

  improper HA N C CB   !stereo CA
  improper HB1 HB2 CA SG  !stereo CB

  dihedral SG  CB  CA  N

  ic N C  *CA CB  0.0 0.0  120.0  0.0 0.0
  ic N CA CB  SG  0.0 0.0  180.0  0.0 0.0
end


residue GLN
  group
    atom N    type=NH1 charge=-0.36 end
    atom HN   type=H   charge= 0.26 end
    atom CA   type=CT  charge= 0.00 end
    atom HA   type=HA  charge= 0.10 end
    atom CB   type=CT  charge=-0.20 end
    atom HB1  type=HA  charge= 0.10 end
    atom HB2  type=HA  charge= 0.10 end
    atom CG   type=CT  charge=-0.20 end
    atom HG1  type=HA  charge= 0.10 end
    atom HG2  type=HA  charge= 0.10 end
    atom CD   type=C   charge= 0.48 end
    atom OE1  type=O   charge=-0.48 end
    atom NE2  type=NH2 charge=-0.52 end
    atom HE21 type=H   charge= 0.26 end
    atom HE22 type=H   charge= 0.26 end
    atom C    type=C   charge= 0.48 end
    atom O    type=O   charge=-0.48 end

  bond N  HN
  bond N  CA     bond CA  HA
  bond CA CB     bond CB  HB1     bond CB  HB2
  bond CB CG     bond CG  HG1     bond CG  HG2
  bond CG CD
  bond CD OE1
  bond CD NE2    bond NE2 HE21    bond NE2 HE22
  bond CA C
  bond C  O

  improper HA  N  C    CB  !stereo CA
  improper CD  CG OE1  NE2
  improper NE2 CD HE21 HE22
  improper CG  CD NE2  HE21
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB CD  !stereo CG

  dihedral CG  CB  CA  N
  dihedral CD  CG  CB  CA
  dihedral OE1 CD  CG  CB

  ic N   C  *CA CB   0.0 0.0  120.0  0.0 0.0
  ic N   CA CB  CG   0.0 0.0  180.0  0.0 0.0
  ic CA  CB CG  CD   0.0 0.0  180.0  0.0 0.0
  ic CB  CG CD  OE1  0.0 0.0  180.0  0.0 0.0
  ic OE1 CG *CD NE2  0.0 0.0  180.0  0.0 0.0
end


residue GLU
  group
    atom N   type=NH1 charge=-0.36 end
    atom HN  type=H   charge= 0.26 end
    atom CA  type=CT  charge= 0.00 end
    atom HA  type=HA  charge= 0.10 end
    atom CB  type=CT  charge=-0.20 end
    atom HB1 type=HA  charge= 0.10 end
    atom HB2 type=HA  charge= 0.10 end
    atom CG  type=CT  charge=-0.45 end
    atom HG1 type=HA  charge= 0.10 end
    atom HG2 type=HA  charge= 0.10 end
    atom CD  type=C   charge= 0.49 end
    atom OE1 type=OC  charge=-0.62 end
    atom OE2 type=OC  charge=-0.62 end
    atom C   type=C   charge= 0.48 end
    atom O   type=O   charge=-0.48 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG HG1     bond CG HG2
  bond CG CD
  bond CD OE1
  bond CD OE2
  bond CA C
  bond C  O

  improper HA N  C   CB    !stereo CA
  improper CD CG OE1 OE2
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB CD  !stereo CG

  dihedral CG  CB  CA  N
  dihedral CD  CG  CB  CA
  dihedral OE1 CD  CG  CB

  ic N   C  *CA CB   0.0 0.0  120.0  0.0 0.0
  ic N   CA CB  CG   0.0 0.0  180.0  0.0 0.0
  ic CA  CB CG  CD   0.0 0.0  180.0  0.0 0.0
  ic CB  CG CD  OE1  0.0 0.0  180.0  0.0 0.0
  ic OE1 CG *CD OE2  0.0 0.0  180.0  0.0 0.0
end


residue GLY
  group
    atom N   type=NH1 charge=-0.36 end
    atom HN  type=H   charge= 0.26 end
    atom CA  type=CT  charge=-0.10 end
    atom HA1 type=HA  charge= 0.10 end
    atom HA2 type=HA  charge= 0.10 end
    atom C   type=C   charge= 0.48 end
    atom O   type=O   charge=-0.48 end

  bond N  HN
  bond N  CA     bond CA HA1     bond CA HA2
  bond CA C
  bond C  O

  improper HA1 HA2 N C  !stereo CA
end


residue HIS
  group
    atom N   type=NH1 charge=-0.36 end
    atom HN  type=H   charge= 0.26 end
    atom CA  type=CT  charge= 0.00 end
    atom HA  type=HA  charge= 0.10 end
    atom CB  type=CT  charge=-0.20 end
    atom HB1 type=HA  charge= 0.10 end
    atom HB2 type=HA  charge= 0.10 end
    atom CG  type=CC  charge= 0.05 end
    atom ND1 type=NA  charge=-0.40 end
    atom HD1 type=H   charge= 0.40 end
    atom CD2 type=CV  charge=-0.14 end
    atom HD2 type=HA  charge= 0.14 end
    atom CE1 type=CR  charge=-0.14 end
    atom HE1 type=HA  charge= 0.14 end
    atom NE2 type=NB  charge=-0.05 end
    atom C   type=C   charge= 0.48 end
    atom O   type=O   charge=-0.48 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB1     bond CB  HB2
  bond CB  CG
  bond CG  ND1    bond ND1 HD1
  bond ND1 CE1    bond CE1 HE1
  bond CG  CD2    bond CD2 HD2
  bond CD2 NE2
  bond CE1 NE2
  bond CA  C
  bond C   O

  improper HA  N   C   CB  !stereo CA
  improper CG  CB  ND1 CD2
  improper ND1 CE1 CG  HD1
  improper CD2 NE2 CG  HD2
  improper CE1 ND1 NE2 HE1
  improper CG  ND1 CE1 NE2
  improper ND1 CE1 NE2 CD2
  improper CE1 NE2 CD2 CG
  improper NE2 CD2 CG  ND1
  improper CD2 CG  ND1 CE1
  improper HB1 HB2 CA CG  !stereo CB

  dihedral CG  CB  CA  N
  dihedral ND1 CG  CB  CA

  ic N   C  *CA CB   0.0 0.0  120.0  0.0  0.0
  ic N   CA CB  CG   0.0 0.0  180.0  0.0  0.0
  ic CA  CB CG  CD2  0.0 0.0   90.0  0.0  0.0
  ic CD2 CB *CG ND1  0.0 0.0  180.0  0.0  0.0
  ic CD2 CG ND1 CE1  0.0 0.0    0.0  0.0  0.0
  ic ND1 CG CD2 NE2  0.0 0.0    0.0  0.0  0.0
end


residue ILE
  group
    atom N    type=NH1 charge=-0.36 end
    atom HN   type=H   charge= 0.26 end
    atom CA   type=CT  charge= 0.00 end
    atom HA   type=HA  charge= 0.10 end
    atom CB   type=CT  charge=-0.10 end
    atom HB   type=HA  charge= 0.10 end
    atom CG1  type=CT  charge=-0.20 end
    atom HG11 type=HA  charge= 0.10 excl = (HG21 HG22 HG23 HD11 HD12 HD13) end
    atom HG12 type=HA  charge= 0.10 excl = (HG21 HG22 HG23 HD11 HD12 HD13) end
    atom CG2  type=CT  charge=-0.30 end
    atom HG21 type=HA  charge= 0.10 excl = (HG11 HG12 HD11 HD12 HD13) end
    atom HG22 type=HA  charge= 0.10 excl = (HG11 HG12 HD11 HD12 HD13) end
    atom HG23 type=HA  charge= 0.10 excl = (HG11 HG12 HD11 HD12 HD13) end
    atom CD1  type=CT  charge=-0.30 end
    atom HD11 type=HA  charge= 0.10 excl = (HG21 HG22 HG23 HG11 HG12) end
    atom HD12 type=HA  charge= 0.10 excl = (HG21 HG22 HG23 HG11 HG12) end
    atom HD13 type=HA  charge= 0.10 excl = (HG21 HG22 HG23 HG11 HG12) end
    atom C    type=C   charge= 0.48 end
    atom O    type=O   charge=-0.48 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB
  bond CB  CG1    bond CG1 HG11    bond CG1 HG12
  bond CB  CG2    bond CG2 HG21    bond CG2 HG22    bond CG2 HG23
  bond CG1 CD1    bond CD1 HD11    bond CD1 HD12    bond CD1 HD13
  bond CA  C
  bond C   O

  improper HA N  C   CB   !stereo CA
  improper HB CA CG2 CG1  !stereo CB
  
  improper HG11 HG12 CB CD1  !stereo CG1
  improper HG21 HG22 CB HG23  !stereo CG2
  improper HD11 HD12 CG1 HD13  !stereo CD1

  dihedral CG1 CB  CA  N
  dihedral CD1 CG1 CB  CA

  ic N  C  *CA CB   0.0 0.0  120.0  0.0 0.0
  ic N  CA CB  CG1  0.0 0.0  180.0  0.0 0.0
  ic N  CA CB  CG2  0.0 0.0   60.0  0.0 0.0
  ic CA CB CG1 CD1   0.0 0.0  180.0  0.0 0.0
end


residue LEU
  group
    atom N    type=NH1 charge=-0.36 end
    atom HN   type=H   charge= 0.26 end
    atom CA   type=CT  charge= 0.00 end
    atom HA   type=HA  charge= 0.10 end
    atom CB   type=CT  charge=-0.20 end
    atom HB1  type=HA  charge= 0.10 end
    atom HB2  type=HA  charge= 0.10 end
    atom CG   type=CT  charge=-0.10 end
    atom HG   type=HA  charge= 0.10 end
    atom CD1  type=CT  charge=-0.30 end
    atom HD11 type=HA  charge= 0.10 excl = (HD21 HD22 HD23 HG) end
    atom HD12 type=HA  charge= 0.10 excl = (HD21 HD22 HD23 HG) end
    atom HD13 type=HA  charge= 0.10 excl = (HD21 HD22 HD23 HG) end
    atom CD2  type=CT  charge=-0.30 end
    atom HD21 type=HA  charge= 0.10 excl = (HD11 HD12 HD13 HG) end
    atom HD22 type=HA  charge= 0.10 excl = (HD11 HD12 HD13 HG) end
    atom HD23 type=HA  charge= 0.10 excl = (HD11 HD12 HD13 HG) end
    atom C    type=C   charge= 0.48 end
    atom O    type=O   charge=-0.48 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB1     bond CB  HB2
  bond CB  CG     bond CG  HG
  bond CG  CD1    bond CD1 HD11    bond CD1 HD12    bond CD1 HD13
  bond CG  CD2    bond CD2 HD21    bond CD2 HD22    bond CD2 HD23
  bond CA  C
  bond C   O

  improper HA   N    C   CB
  improper HG   CB   CD1 CD2   !stereo CG
  improper HB1 HB2 CA CG       !stereo CB
  improper HD11 HD12 CG HD13   !stereo CD1
  improper HD21 HD22 CG HD23   !stereo CD2
    
  dihedral CG  CB  CA  N
  dihedral CD1 CG  CB  CA

  ic N  C  *CA CB   0.0 0.0  120.0  0.0 0.0
  ic N  CA CB  CG   0.0 0.0  180.0  0.0 0.0
  ic CA CB CG  CD1  0.0 0.0  180.0  0.0 0.0
  ic CA CB CG  CD2  0.0 0.0  -60.0  0.0 0.0
end


residue LYS
  group
    atom N   type=NH1 charge=-0.360 end
    atom HN  type=H   charge= 0.260 end
    atom CA  type=CT  charge= 0.000 end
    atom HA  type=HA  charge= 0.100 end
    atom CB  type=CT  charge=-0.200 end
    atom HB1 type=HA  charge= 0.100 end
    atom HB2 type=HA  charge= 0.100 end
    atom CG  type=CT  charge=-0.200 end
    atom HG1 type=HA  charge= 0.100 end
    atom HG2 type=HA  charge= 0.100 end
    atom CD  type=CT  charge=-0.200 end
    atom HD1 type=HA  charge= 0.100 end
    atom HD2 type=HA  charge= 0.100 end
    atom CE  type=CT  charge= 0.305 end
    atom HE1 type=HA  charge= 0.100 end
    atom HE2 type=HA  charge= 0.100 end
    atom NZ  type=NH3 charge=-0.810 end
    atom HZ1 type=HC  charge= 0.435 end
    atom HZ2 type=HC  charge= 0.435 end
    atom HZ3 type=HC  charge= 0.435 end
    atom C   type=C   charge= 0.480 end
    atom O   type=O   charge=-0.480 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG HG1     bond CG HG2
  bond CG CD     bond CD HD1     bond CD HD2
  bond CD CE     bond CE HE1     bond CE HE2
  bond CE NZ     bond NZ HZ1     bond NZ HZ2     bond NZ HZ3
  bond CA C
  bond C  O

  improper HA N C CB      !stereo CA
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB CD  !stereo CG
  improper HD1 HD2 CG CE  !stereo CD
  improper HE1 HE2 CD NZ  !stereo CE
  improper HZ1 HZ2 CE HZ3 !stereo NZ

  dihedral CG  CB  CA  N
  dihedral CD  CG  CB  CA
  dihedral CE  CD  CG  CB
  dihedral NZ  CE  CD  CG

  ic N  C  *CA CB  0.0 0.0  120.0  0.0 0.0
  ic N  CA CB  CG  0.0 0.0  180.0  0.0 0.0
  ic CA CB CG  CD  0.0 0.0  180.0  0.0 0.0
  ic CB CG CD  CE  0.0 0.0  180.0  0.0 0.0
  ic CG CD CE  NZ  0.0 0.0  180.0  0.0 0.0
end

residue MET
  group
    atom N   type=NH1 charge=-0.360 end
    atom HN  type=H   charge= 0.260 end
    atom CA  type=CT  charge= 0.000 end
    atom HA  type=HA  charge= 0.100 end
    atom CB  type=CT  charge=-0.200 end
    atom HB1 type=HA  charge= 0.100 end
    atom HB2 type=HA  charge= 0.100 end
    atom CG  type=CT  charge=-0.115 end
    atom HG1 type=HA  charge= 0.100 end
    atom HG2 type=HA  charge= 0.100 end
    atom SD  type=S   charge=-0.170 end
    atom CE  type=CT  charge=-0.215 end
    atom HE1 type=HA  charge= 0.100 end
    atom HE2 type=HA  charge= 0.100 end
    atom HE3 type=HA  charge= 0.100 end
    atom C   type=C   charge= 0.480 end
    atom O   type=O   charge=-0.480 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG HG1     bond CG HG2
  bond CG SD
  bond SD CE     bond CE HE1     bond CE HE2     bond CE HE3
  bond CA C
  bond C  O

  improper HA N C CB      !stereo CA
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB SD  !stereo CG
  improper HE1 HE2 SD HE3 !stereo methyl

  dihedral CG  CB  CA  N
  dihedral SD  CG  CB  CA
  dihedral CE  SD  CG  CB

  ic N  C  *CA CB  0.0 0.0  120.0  0.0 0.0
  ic N  CA CB  CG  0.0 0.0  180.0  0.0 0.0
  ic CA CB CG  SD  0.0 0.0  180.0  0.0 0.0
  ic CB CG SD  CE  0.0 0.0  180.0  0.0 0.0
end


residue PHE
  group
    atom N   type=NH1 charge=-0.36 end
    atom HN  type=H   charge= 0.26 end
    atom CA  type=CT  charge= 0.00 end
    atom HA  type=HA  charge= 0.10 end
    atom CB  type=CT  charge=-0.16 end
    atom HB1 type=HA  charge= 0.10 end
    atom HB2 type=HA  charge= 0.10 end
    atom CG  type=CA  charge= 0.03 exclude=(CZ) end
    atom CD1 type=CA  charge=-0.16 exclude=(CE2) end
    atom HD1 type=HA  charge= 0.14 end
    atom CD2 type=CA  charge=-0.16 exclude=(CE1) end
    atom HD2 type=HA  charge= 0.14 end
    atom CE1 type=CA  charge=-0.15 exclude=(CD2) end
    atom HE1 type=HA  charge= 0.14 end
    atom CE2 type=CA  charge=-0.15 exclude=(CD1) end
    atom HE2 type=HA  charge= 0.14 end
    atom CZ  type=CA  charge=-0.15 exclude=(CG) end
    atom HZ  type=HA  charge= 0.14 end
    atom C   type=C   charge= 0.48 end
    atom O   type=O   charge=-0.48 end
  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB1     bond CB  HB2
  bond CB  CG
  bond CG  CD1    bond CD1 HD1
  bond CG  CD2    bond CD2 HD2
  bond CD1 CE1    bond CE1 HE1
  bond CD2 CE2    bond CE2 HE2
  bond CE1 CZ     bond CZ  HZ
  bond CE2 CZ
  bond CA  C
  bond C   O

  improper HA  N   C   CB !stereo CA
  improper HB1 HB2 CA CG  !stereo CB

! Hs and CB around the ring
  improper HD2 CD2 CE2 CZ
  improper HE2 CE2 CZ  CE1
  improper HZ  CZ  CE1 CD1
  improper HE1 CE1 CD1 CG
  improper HD1 CD1 CG  CD2
  improper CB  CG  CD2 CE2

! around the ring
  improper CG  CD1 CE1 CZ 
  improper CD1 CE1 CZ  CE2
  improper CE1 CZ  CE2 CD2
  improper CZ  CE2 CD2 CG
  improper CE2 CD2 CG  CD1
  improper CD2 CG  CD1 CE1

  dihedral CG  CB  CA  N
  dihedral CD1 CG  CB  CA

  ic N   C   *CA CB   0.0 0.0  120.0  0.0 0.0
  ic N   CA  CB  CG   0.0 0.0  180.0  0.0 0.0
  ic CA  CB  CG  CD1  0.0 0.0   90.0  0.0 0.0
  ic CD1 CB  *CG CD2  0.0 0.0  180.0  0.0 0.0
  ic CD1 CG  CD2 CE2  0.0 0.0    0.0  0.0 0.0
  ic CD2 CG  CD1 CE1  0.0 0.0    0.0  0.0 0.0
  ic CG  CD1 CE1 CZ   0.0 0.0    0.0  0.0 0.0
end


residue PRO
  group
    atom N   type=N  charge=-0.36 end
    atom CA  type=CT charge= 0.00 end
    atom HA  type=HA charge= 0.10 end
    atom CB  type=CT charge=-0.20 end
    atom HB1 type=HA charge= 0.10 end
    atom HB2 type=HA charge= 0.10 end
    atom CG  type=CP charge=-0.20 end
    atom HG1 type=HA charge= 0.10 end
    atom HG2 type=HA charge= 0.10 end
    atom CD  type=CP charge= 0.06 end
    atom HD2 type=HA charge= 0.10 end   ! the order of the two h's
    atom HD1 type=HA charge= 0.10 end   ! is inverted for the sake of
    atom C   type=C  charge= 0.48 end   ! hbuild
    atom O   type=O  charge=-0.48 end

  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG HG1     bond CG HG2
  bond CG CD     bond CD HD1     bond CD HD2
  bond CD N
  bond CA C
  bond C  O

  improper HA N C CB      !chiral CA
  improper N CA CB CG     !ring pucker
  improper HB1 HB2 CA CG  !stereo CB
  improper HG1 HG2 CB CD  !stereo CG
  improper HD1 HD2 CG N   !stereo CD

  ic C N  *CA CB  0.0 0.0 -120.0  0.0 0.0
  ic N CA CB  CG  0.0 0.0    0.0  0.0 0.0
end


residue SER
  group
    atom N   type=NH1 charge=-0.36 end
    atom HN  type=H   charge= 0.26 end
    atom CA  type=CT  charge= 0.00 end
    atom HA  type=HA  charge= 0.10 end
    atom CB  type=CT  charge= 0.08 end
    atom HB1 type=HA  charge= 0.10 end
    atom HB2 type=HA  charge= 0.10 end
    atom OG  type=OH  charge=-0.68 end
    atom HG  type=H   charge= 0.40 end
    atom C   type=C   charge= 0.48 end
    atom O   type=O   charge=-0.48 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB OG     bond OG HG
  bond O  C
  bond C  CA

  improper HA N C CB      !stereo CA
  improper HB1 HB2 CA OG  !stereo CB

  dihedral OG  CB  CA  N

  ic N C  *CA CB  0.0 0.0  120.0  0.0 0.0
  ic N CA CB  OG  0.0 0.0  180.0  0.0 0.0
end


residue THR
  group
    atom N    type=NH1 charge=-0.36 end
    atom HN   type=H   charge= 0.26 end
    atom CA   type=CT  charge= 0.00 end
    atom HA   type=HA  charge= 0.10 end
    atom CB   type=CT  charge= 0.18 end
    atom HB   type=HA  charge= 0.10 end
    atom OG1  type=OH  charge=-0.68 end
    atom HG1  type=H   charge= 0.40 end
    atom CG2  type=CT  charge=-0.30 end
    atom HG21 type=HA  charge= 0.10 end
    atom HG22 type=HA  charge= 0.10 end
    atom HG23 type=HA  charge= 0.10 end
    atom C    type=C   charge= 0.48 end
    atom O    type=O   charge=-0.48 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB
  bond CB OG1    bond OG1 HG1
  bond CB CG2    bond CG2 HG21    bond CG2 HG22    bond CG2 HG23
  bond CA C
  bond C  O

  improper HA N  C   CB        !stereo CA
  improper HB CA OG1 CG2       !stereo CB
  improper HG21 HG22 CB HG23  !stereo methyl

  dihedral OG1 CB  CA  N

  ic N C  *CA CB   0.0 0.0  120.0  0.0 0.0
  ic N CA CB  OG1  0.0 0.0  180.0  0.0 0.0
  ic N CA CB  CG2  0.0 0.0   60.0  0.0 0.0
end


residue TRP
  group
    atom N   type=NH1 charge=-0.36 end
    atom HN  type=H   charge= 0.26 end
    atom CA  type=CT  charge= 0.00 end
    atom HA  type=HA  charge= 0.10 end
    atom CB  type=CT  charge=-0.20 end
    atom HB1 type=HA  charge= 0.10 end
    atom HB2 type=HA  charge= 0.10 end
    atom CG  type=CX  charge=-0.04 end
    atom CD1 type=CW  charge=-0.01 end
    atom HD1 type=HA  charge= 0.14 end
    atom CD2 type=CB  charge=-0.05 exclude=(CH2) end
    atom NE1 type=NA charge=-0.49 end
    atom HE1 type=H   charge= 0.40 end
    atom CE2 type=CN  charge= 0.13 exclude=(CZ3) end
    atom CE3 type=CA  charge=-0.16 exclude=(CZ2) end
    atom HE3 type=HA  charge= 0.14 end
    atom CZ2 type=CA  charge=-0.15 exclude=(CE3) end
    atom HZ2 type=HA  charge= 0.14 end
    atom CZ3 type=CA  charge=-0.15 exclude=(CE2) end
    atom HZ3 type=HA  charge= 0.14 end
    atom CH2 type=CA  charge=-0.18 exclude=(CD2) end
    atom HH2 type=HA  charge= 0.14 end
    atom C   type=C   charge= 0.48 end
    atom O   type=O   charge=-0.48 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB1     bond CB  HB2
  bond CB  CG
  bond CG  CD1    bond CD1 HD1
  bond CG  CD2
  bond CD1 NE1    bond NE1 HE1
  bond NE1 CE2
  bond CD2 CE2
  bond CD2 CE3    bond CE3 HE3
  bond CE2 CZ2    bond CZ2 HZ2
  bond CE3 CZ3    bond CZ3 HZ3
  bond CZ2 CH2    bond CH2 HH2
  bond CZ3 CH2
  bond CA  C
  bond C   O

! chirality
  improper HA  N   C   CB !stereo CA
  improper HB1 HB2 CA CG  !stereo CB

! around the 6-ring
  improper CD2 CE2 CZ2 CH2
  improper CE2 CZ2 CH2 CZ3
  improper CZ2 CH2 CZ3 CE3
  improper CH2 CZ3 CE3 CD2
  improper CZ3 CE3 CD2 CE2
  improper CE3 CD2 CE2 CZ2

! link 5-ring to 6-ring
  improper CD1 NE1 CE2 CZ2
  improper CD1 CG  CD2 CE3
  improper NE1 CE2 CZ2 CH2
  improper NE1 CE2 CD2 CE3
  improper CG  CD2 CE3 CZ3
  improper CG  CD2 CE2 CZ2

! 6-ring hydrogens
  improper HZ2 CZ2 CH2 CZ3
  improper HH2 CH2 CZ3 CE3
  improper HZ3 CZ3 CH2 CZ2
  improper HE3 CE3 CZ3 CH2

! 5-ring hydrogens and CB
  improper HE1 NE1 CE2 CD2
  improper HD1 CD1 NE1 CE2
  improper CB  CG  CD2 CE2

  dihedral CG  CB  CA  N
  dihedral CD1 CG  CB  CA

  ic N   C   *CA  CB   0.0 0.0  120.0  0.0 0.0
  ic N   CA  CB   CG   0.0 0.0  180.0  0.0 0.0
  ic CA  CB  CG   CD2  0.0 0.0  180.0  0.0 0.0
  ic CD2 CB  *CG  CD1  0.0 0.0  180.0  0.0 0.0
  ic CD1 CG  CD2  CE2  0.0 0.0    0.0  0.0 0.0
  ic CD2 CG  CD1  NE1  0.0 0.0    0.0  0.0 0.0
  ic CE2 CG  *CD2 CE3  0.0 0.0  180.0  0.0 0.0
  ic CE2 CD2 CE3  CZ3  0.0 0.0    0.0  0.0 0.0
  ic CE3 CD2 CE2  CZ2  0.0 0.0    0.0  0.0 0.0
  ic CD2 CE2 CZ2  CH2  0.0 0.0    0.0  0.0 0.0
end

residue TYR
  group
    atom N   type=NH1 charge=-0.36 end
    atom HN  type=H   charge= 0.26 end
    atom CA  type=CT  charge= 0.00 end
    atom HA  type=HA  charge= 0.10 end
    atom CB  type=CT  charge=-0.20 end
    atom HB1 type=HA  charge= 0.10 end
    atom HB2 type=HA  charge= 0.10 end
    atom CG  type=CA  charge= 0.00 exclude=(CZ) end
    atom CD1 type=CA  charge=-0.14 exclude=(CE2) end
    atom HD1 type=HA  charge= 0.14 end
    atom CD2 type=CA  charge=-0.14 exclude=(CE1) end
    atom HD2 type=HA  charge= 0.14 end
    atom CE1 type=CA  charge=-0.14 exclude=(CD2) end
    atom HE1 type=HA  charge= 0.14 end
    atom CE2 type=CA  charge=-0.14 exclude=(CD1) end
    atom HE2 type=HA  charge= 0.14 end
    atom CZ  type=C   charge= 0.20 exclude=(CG) end
    atom OH  type=OH  charge=-0.60 end
    atom HH  type=H   charge= 0.40 end
    atom C   type=C   charge= 0.48 end
    atom O   type=O   charge=-0.48 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB1     bond CB  HB2
  bond CB  CG
  bond CG  CD1    bond CD1 HD1
  bond CG  CD2    bond CD2 HD2
  bond CD1 CE1    bond CE1 HE1
  bond CD2 CE2    bond CE2 HE2
  bond CE1 CZ
  bond CE2 CZ
  bond CZ  OH     bond OH  HH
  bond CA  C
  bond C   O

! chirality
  improper HA  N   C   CB !stereo CA
  improper HB1 HB2 CA CG  !stereo CB

! Hs, OH, and CB around the ring
  improper HD2 CD2 CE2 CZ
  improper HE2 CE2 CZ  CE1
  improper OH  CZ  CE1 CD1
  improper HE1 CE1 CD1 CG
  improper HD1 CD1 CG  CD2
  improper CB  CG  CD2 CE2

! around the ring
  improper CG  CD1 CE1 CZ 
  improper CD1 CE1 CZ  CE2
  improper CE1 CZ  CE2 CD2
  improper CZ  CE2 CD2 CG
  improper CE2 CD2 CG  CD1
  improper CD2 CG  CD1 CE1

  dihedral CG  CB  CA  N
  dihedral CD1 CG  CB  CA

  ic N   C   *CA CB   0.0 0.0  120.0  0.0 0.0
  ic N   CA  CB  CG   0.0 0.0  180.0  0.0 0.0
  ic CA  CB  CG  CD1  0.0 0.0   90.0  0.0 0.0
  ic CD1 CB  *CG CD2  0.0 0.0  180.0  0.0 0.0
  ic CD2 CG  CD1 CE1  0.0 0.0    0.0  0.0 0.0
  ic CD1 CG  CD2 CE2  0.0 0.0    0.0  0.0 0.0
  ic CG  CD1 CE1 CZ   0.0 0.0    0.0  0.0 0.0
  ic CE2 CE1 *CZ OH   0.0 0.0  180.0  0.0 0.0
end


residue VAL
  group
    atom N    type=NH1 charge=-0.36 end
    atom HN   type=H   charge= 0.26 end
    atom CA   type=CT  charge= 0.00 end
    atom HA   type=HA  charge= 0.10 end
    atom CB   type=CT  charge=-0.10 end
    atom HB   type=HA  charge= 0.10 end
    atom CG1  type=CT  charge=-0.30 end
    atom HG11 type=HA  charge= 0.10 excl = (HG21 HG22 HG23) end
    atom HG12 type=HA  charge= 0.10 excl = (HG21 HG22 HG23) end
    atom HG13 type=HA  charge= 0.10 excl = (HG21 HG22 HG23) end
    atom CG2  type=CT  charge=-0.30 end
    atom HG21 type=HA  charge= 0.10 excl = (HG11 HG12 HG13) end
    atom HG22 type=HA  charge= 0.10 excl = (HG11 HG12 HG13) end
    atom HG23 type=HA  charge= 0.10 excl = (HG11 HG12 HG13) end
    atom C    type=C   charge= 0.48 end
    atom O    type=O   charge=-0.48 end

  bond N   HN
  bond N   CA     bond CA  HA
  bond CA  CB     bond CB  HB
  bond CB  CG1    bond CG1 HG11    bond CG1 HG12    bond CG1 HG13
  bond CB  CG2    bond CG2 HG21    bond CG2 HG22    bond CG2 HG23
  bond CA  C
  bond C   O

  improper HA N C CB         !stereo CA
  improper HB   CA   CG1 CG2 !stereo CB
  improper HG11 HG12 CB HG13    !stereo G1
  improper HG21 HG22 CB HG23    !stereo G2

  dihedral CG1 CB  CA  N

  ic N C  *CA CB   0.0 0.0  120.0  0.0 0.0
  ic N CA CB  CG1  0.0 0.0  180.0  0.0 0.0
  ic N CA CB  CG2  0.0 0.0  -60.0  0.0 0.0
end




residue CHEX  !! ADDED BY MN
  group
    atom N   type=NH1 charge=-0.360 end
    atom HN  type=H   charge= 0.260 end
    atom CA  type=CT  charge= 0.000 end
    atom HA  type=HA  charge= 0.100 end
    atom CB  type=CT  charge=-0.200 end
    atom HB1 type=HA  charge= 0.100 end
    atom HB2 type=HA  charge= 0.100 end
    atom CG  type=CT  charge=-0.200 end
    atom HG  type=HA  charge= 0.100 end
    atom CD1  type=CT  charge=-0.200 end
    atom HD11 type=HA  charge= 0.100 end
    atom HD12 type=HA  charge= 0.100 end
    atom CD2  type=CT  charge=-0.200 end
    atom HD21 type=HA  charge= 0.100 end
    atom HD22 type=HA  charge= 0.100 end
    atom CE1  type=CT  charge=-0.200 end
    atom HE11 type=HA  charge= 0.100 end
    atom HE12 type=HA  charge= 0.100 end
    atom CE2  type=CT  charge=-0.200 end
    atom HE21 type=HA  charge= 0.100 end
    atom HE22 type=HA  charge= 0.100 end
    atom CZ  type=CT  charge=-0.200 end
    atom HZ1 type=HA  charge= 0.100 end
    atom HZ2 type=HA  charge= 0.100 end
    atom C   type=C   charge= 0.480 end
    atom O   type=O   charge=-0.480 end

  bond N  HN
  bond N  CA     bond CA HA
  bond CA CB     bond CB HB1     bond CB HB2
  bond CB CG     bond CG CD1     bond CG CD2     bond CG HG    
  bond CD1 CE1   bond CD1 HD11    bond CD1 HD12   
  bond CD2 CE2   bond CD2 HD21    bond CD2 HD22
  bond CE1 CZ    bond CE1 HE11    bond CE1 HE12
  bond CE2 CZ    bond CE2 HE21    bond CE2 HE22
  bond CZ HZ1    bond CZ  HZ2
  bond CA C
  bond C  O

  improper HA N C CB          !stereo CA
  improper HB1 HB2 CA CG  !stereo CB
  improper HG CB CD1 CD2      !stereo CG
  improper HD11 HD12 CG CE1      !stereo CD1
  improper HD21 HD22 CG CE2      !stereo CD2
  improper HE11 HE12 CE1 CZ      !stereo CE1
  improper HE11 HE22 CE2 CZ      !stereo CE2
  improper HZ1 HZ2 CE1 CE2  !stereo CZ

  dihedral CG  CB  CA  N
  dihedral CD1  CG  CB  CA

  ic N  C  *CA CB     0.0 0.0   120.0  0.0 0.0
  ic N  CA CB  CG     0.0 0.0   180.0  0.0 0.0
  ic CA CB CG  CD1    0.0 0.0   180.0  0.0 0.0
  ic CB CG CD1 CE1    0.0 0.0   180.0  0.0 0.0
  ic CG CD1 CE1 CZ    0.0 0.0   -60.0  0.0 0.0
  ic CD1 CE1 CZ CE2   0.0 0.0    60.0  0.0 0.0 
  ic CE1 CZ CE2 CD2   0.0 0.0   -60.0  0.0 0.0
end





residue ACE
  group
    atom CA  type=CT charge=-0.30 end
    atom HA1 type=HA charge= 0.10 end
    atom HA2 type=HA charge= 0.10 end
    atom HA3 type=HA charge= 0.10 end
    atom C   type=C  charge= 0.48 end
    atom O   type=O  charge=-0.48 end

  bond C  CA     bond CA HA1     bond CA HA2     bond CA HA3
  bond C  O
  improper HA1 HA2 C HA3  !stereo methyl
end




presidue NTER                      ! patch as "NTER - *" to any except PRO
  modify    atom +CA           charge= 0.22 end
  group
    modify    atom +N   type=NH3 charge=-0.10 end
    delete    atom +HN                        end
    add       atom +HT1 type=HC  charge= 0.26 end
    add       atom +HT2 type=HC  charge= 0.26 end
    add       atom +HT3 type=HC  charge= 0.26 end

  add bond +HT1 +N
  add bond +HT2 +N
  add bond +HT3 +N

  add angle +HT1 +N +HT2
  add angle +HT2 +N +HT3
  add angle +HT2 +N +CA
  add angle +HT1 +N +HT3
  add angle +HT1 +N +CA
  add angle +HT3 +N +CA
  add improper +HT1 +HT2 +CA +HT3  !stereo N ???

end


presidue PROP                        ! N-terminal for PRO: "PROP - PRO"
  modify    atom +CD           charge= 0.10 end
  modify    atom +CA           charge= 0.10 end
  modify    atom +N   type=NH3 charge=-0.02 end
  add       atom +HT1 type=HC  charge= 0.26 end
  add       atom +HT2 type=HC  charge= 0.26 end

  add bond +HT1 +N
  add bond +HT2 +N

  add angle +HT1 +N +HT2
  add angle +HT2 +N +CA
  add angle +HT1 +N +CD
  add angle +HT1 +N +CA
  add angle +CD  +N +HT2
  add improper +HT1 +HT2 +CA +CD  !stereo N 
end
 

presidue CTER               ! C-terminal for all amino acids "* - CTER"
  group
    modify    atom -C           charge= 0.14 end
    delete    atom -O                        end
    add       atom -OT1 type=OC charge=-0.57 end
    add       atom -OT2 type=OC charge=-0.57 end

  add bond -C -OT1
  add bond -C -OT2

  add angle -CA  -C -OT1
  add angle -CA  -C -OT2
  add angle -OT1 -C -OT2

  add improper -C -CA -OT2 -OT1
 
  add ic -N   -CA -C  -OT2  0.0 0.0  180.0  0.0 0.0
  add ic -OT2 -CA *-C -OT1  0.0 0.0  180.0  0.0 0.0
end


presidue CTN                  ! C-terminal for all, CONH2 at end "* - CTN"
  group
    modify    atom -C           charge= 0.48 end
    modify    atom -O           charge=-0.48 end
    add       atom -NT type=NH2 charge=-0.52 end
    add       atom -H1 type=H   charge= 0.26 end
    add       atom -H2 type=H   charge= 0.26 end

  add bond -C  -NT
  add bond -NT -H1
  add bond -NT -H2

  add angle -CA -C  -NT
  add angle -O  -C  -NT
  add angle -CA -C  -O
  add angle -C  -NT -H1
  add angle -C  -NT -H2
  add angle -H1 -NT -H2

  add improper -C  -CA -NT -O
  add improper -C  -NT -O  -H1
  add improper -NT -H1 -H2 -C
end
 

presidue PEPT     ! PEPTide bond link, for all except the  *(-) - (+)PRO link
                  ! "*(-) - PEPT - (+)*:
  add bond -C +N

  add angle -CA -C +N
  add angle -O  -C +N
  add angle -C  +N +CA
  add angle -C  +N +HN

  add improper -O  -C  +N  +CA                     ! planar -C    fixed for DG by JK
  add improper +HN +N  -C  -CA                     ! planar +N             
  add improper -CA -C  +N  +CA                     ! planar peptide        "

  add ic -C  +CA *+N +H   0.0 0.0  180.0  0.0 0.0  ! planar -C
  add ic -C  +N  +CA +C   0.0 0.0  180.0  0.0 0.0
  add ic -N  -CA -C  +N   0.0 0.0  180.0  0.0 0.0
  add ic -CA -C  +N  +CA  0.0 0.0  180.0  0.0 0.0
  add ic +N  -CA *-C -O   0.0 0.0  180.0  0.0 0.0  ! planar +N
end
 

presidue PEPP     ! for  ...*(-) - (+)PRO  link, same as PEPT except
                  ! replacement H by CD and improper +N +CA +CD -C
  add bond -C +N

  add angle -CA -C +N
  add angle -O  -C +N
  add angle -C  +N +CA
  add angle -C  +N +CD

  add improper -O  -C  +N  +CA                      ! planar -C
  add improper +CA +N  -C  +CD                      ! planar +N modified
  add improper -CA  -C +N +CA                       ! planar peptide

  add ic -C  +CA *+N +CD  0.0 0.0  180.0  0.0 0.0  ! planar +N
  add ic -C  +N  +CA +C   0.0 0.0  180.0  0.0 0.0
  add ic -N  -CA -C  +N   0.0 0.0  180.0  0.0 0.0
  add ic -CA -C  +N  +CA  0.0 0.0  180.0  0.0 0.0
  add ic +N  -CA *-C -O   0.0 0.0  180.0  0.0 0.0  ! planar -C
end


presidue DISU                ! disulfide bridge  ...CYS - DISU - CYS...
  group
    delete    atom 1HG               end
    modify    atom 1CB  charge= 0.20 end
    modify    atom 1SG  charge=-0.20 end
  group
    delete    atom 2HG               end
    modify    atom 2CB  charge= 0.20 end
    modify    atom 2SG  charge=-0.20 end

  add bond 1SG 2SG

  add angle 1CB 1SG 2SG
  add angle 1SG 2SG 2CB

  add ic 1CA 1CB 1SG 2SG  0.0 0.0  180.0  0.0 0.0
  add ic 1CB 1SG 2SG 2CB  0.0 0.0  180.0  0.0 0.0
  add ic 1SG 2SG 2CB 2CA  0.0 0.0  180.0  0.0 0.0
end


presidue DISN                ! disulfide bridge  ...CYS - DISU - CYS...
                             ! w/o the actual bond
  group                       
    delete    atom 1HG               end
    modify    atom 1CB  charge= 0.20 end
    modify    atom 1SG  charge=-0.20 end
  group
    delete    atom 2HG               end
    modify    atom 2CB  charge= 0.20 end
    modify    atom 2SG  charge=-0.20 end
end


presidue LTOD                        ! change from L to D amino acid
  delete improper HA N C CB
  add    improper HA C N CB
end

presidue CISP                        ! change from trans to cis peptide bond
  delete improper -CA -C  +N  +CA
  add    improper -C  -CA +N  +CA
  delete improper -O  -C  +N  +CA
  add    improper -C  -O  +N  +CA
  delete improper +HN +N  -C  -CA                     
  add    improper +N  +HN -C  -CA
end

presidue CIPP                        ! change from trans to cis peptide bond
  delete improper -CA -C  +N  +CA
  add    improper -C  -CA +N  +CA
  delete improper -O  -C  +N  +CA
  add    improper -C  -O  +N  +CA
  delete improper +CA +N  -C  +CD                     
  add    improper +CA +N  +CD -C
end

presidue HISE     !  change protonation of HIS
                  !  has to be patched as REFErence=NIL=<selection>
                       
  modify atom ND1 type=NA  charge=-0.05 end
  delete atom HD1                       end
  modify atom NE2 type=NB  charge=-0.40 end
  add    atom HE2 type=H   charge= 0.40 end

  add bond HE2 NE2
  add angle HE2 NE2 CD2
  add angle HE2 NE2 CE1
  add improper NE2 CE1 CD2 HE2

end {HISE}


set mess=$prev_messages echo=$prev_echo end
