 REMARKS FILENAME="parallh22x-eefx.pro"

 SET ECHO=FALSE END

eval ($protein_par_vers="charmm22-eef-1.0")

{ >>>>>>>>>> Developmental Parameter File for Proteins <<<<<<<<<<
 >>>>>>>>>>>>>>>>> Using All Hydrogens (ALLH) <<<<<<<<<<<<<<<<<<
 >>>>>>>>>>>>>>>>>>>>>>> May 1993 <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
 >>>>>>> Direct comments to Alexander D. MacKerell Jr. <<<<<<<<<
 >>>>>> 410-706-7442 or email: alex@mmiris.ab.umd.edu  <<<<<<<<<
 These files are a beta release; additional parameter development
 and testing may lead to alteration of the contents.}


 BOND  C    C        600.000 {SD=      .022}      1.335 ! ALLOW ARO HEM 
      ! Heme vinyl substituent (KK, from propene (JCS)) 
 BOND  CA   CA       305.000 {SD=      .031}      1.375 ! ALLOW   ARO 
      ! benzene, JES 8/25/89 
 BOND  CP1  C        250.000 {SD=      .034}      1.490 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  CP1  CC       250.000 {SD=      .034}      1.490 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  CP1  CD       200.000 {SD=      .038}      1.490 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  CP2  CP1      222.500 {SD=      .036}      1.527 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  CP2  CP2      222.500 {SD=      .036}      1.537 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  CP3  CP2      222.500 {SD=      .036}      1.537 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  CPB  C        450.000 {SD=      .026}      1.380 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 BOND  CPB  CPA      299.800 {SD=      .031}      1.443 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 BOND  CPB  CPB      340.700 {SD=      .029}      1.346 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 BOND  CPH1 CPH1     410.000 {SD=      .027}      1.360 ! ALLOW ARO 
      ! histidine, adm jr., 6/27/90 
 BOND  CPM  CPA      360.000 {SD=      .029}      1.372 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 BOND  CPT  CA       305.000 {SD=      .031}      1.368 ! ALLOW   ARO 
      ! adm jr., 12/30/91, for jwk 
 BOND  CPT  CPT      360.000 {SD=      .029}      1.400 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 BOND  CT1  C        250.000 {SD=      .034}      1.490 ! ALLOW   ALI PEP POL ARO 
      ! Ala Dipeptide ab initio calc's (LK) fixed from 10/90 (5/91) 
 BOND  CT1  CC       200.000 {SD=      .038}      1.522 ! ALLOW   POL 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 BOND  CT1  CD       200.000 {SD=      .038}      1.522 ! ALLOW   POL 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 BOND  CT1  CT1      222.500 {SD=      .036}      1.500 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 BOND  CT2  C        250.000 {SD=      .034}      1.490 ! ALLOW   ALI PEP POL ARO 
      ! Ala Dipeptide ab initio calc's (LK) fixed from 10/90 (5/91) 
 BOND  CT2  CA       230.000 {SD=      .036}      1.490 ! ALLOW   ALI ARO 
      ! phe,tyr, JES 8/25/89 
 BOND  CT2  CC       200.000 {SD=      .038}      1.522 ! ALLOW   POL 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 BOND  CT2  CD       200.000 {SD=      .038}      1.522 ! ALLOW   POL 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 BOND  CT2  CPB      230.000 {SD=      .036}      1.490 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 BOND  CT2  CPH1     229.630 {SD=      .036}      1.500 ! ALLOW ARO 
      ! his, adm jr., 7/22/89, FC from CT2CT, BL from crystals 
 BOND  CT2  CT1      222.500 {SD=      .036}      1.538 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 BOND  CT2  CT2      222.500 {SD=      .036}      1.530 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 BOND  CT3  C        250.000 {SD=      .034}      1.490 ! ALLOW   ALI PEP POL ARO 
      ! Ala Dipeptide ab initio calc's (LK) fixed from 10/90 (5/91) 
 BOND  CT3  CA       230.000 {SD=      .036}      1.490 ! ALLOW   ALI ARO 
      ! toluene, adm jr. 3/7/92 
 BOND  CT3  CC       200.000 {SD=      .038}      1.522 ! ALLOW   POL 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 BOND  CT3  CD       200.000 {SD=      .038}      1.522 ! ALLOW   POL 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 BOND  CT3  CPB      230.000 {SD=      .036}      1.490 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 BOND  CT3  CPH1     229.630 {SD=      .036}      1.500 ! ALLOW ARO 
      ! his, adm jr., 7/22/89, FC from CT2CT, BL from crystals 
 BOND  CT3  CS       190.000 {SD=      .039}      1.531 ! ALLOW   SUL 
      ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92 
 BOND  CT3  CT1      222.500 {SD=      .036}      1.538 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 BOND  CT3  CT2      222.500 {SD=      .036}      1.528 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 BOND  CT3  CT3      222.500 {SD=      .036}      1.530 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 BOND  CY   CA       350.000 {SD=      .029}      1.365 ! ALLOW   ARO 
      !adm jr., 5/08/91, indole CCDB structure search 
 BOND  CY   CPT      350.000 {SD=      .029}      1.440 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 BOND  CY   CT2      230.000 {SD=      .036}      1.510 ! ALLOW ARO 
      !JWK  Kb from alkane freq.. b0 from TRP crystal 
 BOND  FE   CM       258.000 {SD=      .034}      1.900 ! ALLOW HEM 
      ! Heme (6-liganded): CO ligand (KK 05/13/91) 
 BOND  FE   CPM         .000 {SD=999999.000}      3.381 ! ALLOW HEM 
      ! Heme (6-liganded): for "ic para" only (KK 05/13/91) 
 BOND  H    CD       330.000 {SD=      .030}      1.110 ! ALLOW   PEP POL ARO 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 BOND  HA   C        330.000 {SD=      .030}      1.100 ! ALLOW ARO HEM 
      ! Heme vinyl substituent (KK, from propene (JCS)) 
 BOND  HA   CA       340.000 {SD=      .030}      1.083 ! ALLOW ARO 
      ! trp, adm jr., 10/02/89 
 BOND  HA   CC       317.130 {SD=      .031}      1.100 ! ALLOW POL 
      ! adm jr., 5/13/91, formamide geometry and vibrations 
 BOND  HA   CP2      309.000 {SD=      .031}      1.111 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  HA   CP3      309.000 {SD=      .031}      1.111 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  HA   CPM      367.600 {SD=      .028}      1.090 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 BOND  HA   CS       300.000 {SD=      .031}      1.111 ! ALLOW   SUL 
      ! methylthiolate 6-31+G* geom/freq, adm jr., 6/1/92 
 BOND  HA   CT1      309.000 {SD=      .031}      1.111 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 BOND  HA   CT2      309.000 {SD=      .031}      1.111 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 BOND  HA   CT3      322.000 {SD=      .030}      1.111 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 BOND  HA   CY       330.000 {SD=      .030}      1.080 ! ALLOW   ARO 
      ! JWK 05/14/91 new r0 from indole 
 BOND  HB   CP1      330.000 {SD=      .030}      1.080 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  HB   CT1      330.000 {SD=      .030}      1.080 ! ALLOW   PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 BOND  HB   CT2      330.000 {SD=      .030}      1.080 ! ALLOW   PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 BOND  HB   CT3      330.000 {SD=      .030}      1.080 ! ALLOW   PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 BOND  HP   CA       340.000 {SD=      .030}      1.080 ! ALLOW   ARO 
      ! phe,tyr JES 8/25/89 
 BOND  HP   CY       350.000 {SD=      .029}      1.080 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 BOND  HR1  CPH1     375.000 {SD=      .028}      1.083 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 BOND  HR1  CPH2     340.000 {SD=      .030}      1.090 ! ALLOW ARO 
      ! his, adm jr., 6/28/29 
 BOND  HR2  CPH2     333.000 {SD=      .030}      1.070 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 BOND  HR3  CPH1     365.000 {SD=      .028}      1.083 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 BOND  HT   HT          .000 {SD=999999.000}      1.514 ! ALLOW WAT 
      ! FROM TIPS3P GEOMETRY (FOR SHAKE/W PARAM) 
 BOND  N    C        260.000 {SD=      .034}      1.300 ! ALLOW PEP POL ARO PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  N    CP1      320.000 {SD=      .030}      1.434 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  N    CP3      320.000 {SD=      .030}      1.455 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  NC2  C        463.000 {SD=      .025}      1.365 ! ALLOW   PEP POL ARO 
      ! 403.0->463.0, 1.305->1.365 guanidinium (KK) 
 BOND  NC2  CT2      261.000 {SD=      .034}      1.490 ! ALLOW   ALI POL 
      ! arg, (DS) 
 BOND  NC2  CT3      261.000 {SD=      .034}      1.490 ! ALLOW   ALI POL 
      ! methylguanidinium, adm jr., 3/26/92 
 BOND  NC2  HC       455.000 {SD=      .026}      1.000 ! ALLOW   POL 
      ! 405.0->455.0 GUANIDINIUM (KK) 
 BOND  NH1  C        370.000 {SD=      .028}      1.345 ! ALLOW   PEP POL ARO 
      ! Alanine Dipeptide ab initio calc's (LK) 
 BOND  NH1  CT1      320.000 {SD=      .030}      1.430 ! ALLOW   ALI PEP POL ARO 
      ! NMA Gas & Liquid Phase IR Spectra (LK) 
 BOND  NH1  CT2      320.000 {SD=      .030}      1.430 ! ALLOW   ALI PEP POL ARO 
      ! NMA Gas & Liquid Phase IR Spectra (LK) 
 BOND  NH1  CT3      320.000 {SD=      .030}      1.430 ! ALLOW   ALI PEP POL ARO 
      ! NMA Gas & Liquid Phase IR Spectra (LK) 
 BOND  NH1  H        440.000 {SD=      .026}       .997 ! ALLOW   PEP POL ARO 
      ! Alanine Dipeptide ab initio calc's (LK) 
 BOND  NH1  HC       405.000 {SD=      .027}       .980 ! ALLOW   PEP POL ARO 
      ! (DS) 
 BOND  NH2  CC       430.000 {SD=      .026}      1.360 ! ALLOW   PEP POL ARO 
      ! adm jr. 4/10/91, acetamide 
 BOND  NH2  CT3      240.000 {SD=      .035}      1.455 ! ALLOW   POL 
      ! methylamine geom/freq, adm jr., 6/2/92 
 BOND  NH2  H        480.000 {SD=      .025}      1.000 ! ALLOW   POL 
      ! adm jr. 8/13/90 acetamide geometry and vibrations 
 BOND  NH2  HC       460.000 {SD=      .025}      1.000 ! ALLOW   POL 
      ! methylamine geom/freq, adm jr., 6/2/92 
 BOND  NH3  CT1      200.000 {SD=      .038}      1.480 ! ALLOW   ALI POL 
      ! new stretch and bend; methylammonium (KK 03/10/92) 
 BOND  NH3  CT2      200.000 {SD=      .038}      1.480 ! ALLOW   ALI POL 
      ! new stretch and bend; methylammonium (KK 03/10/92) 
 BOND  NH3  CT3      200.000 {SD=      .038}      1.480 ! ALLOW   ALI POL 
      ! new stretch and bend; methylammonium (KK 03/10/92) 
 BOND  NH3  HC       403.000 {SD=      .027}      1.040 ! ALLOW   POL 
      ! new stretch and bend; methylammonium (KK 03/10/92) 
 BOND  NP   CP1      320.000 {SD=      .030}      1.485 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  NP   CP3      320.000 {SD=      .030}      1.502 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  NP   HC       460.000 {SD=      .025}      1.006 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 BOND  NPH  CPA      377.200 {SD=      .028}      1.376 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 BOND  NPH  FE       270.200 {SD=      .033}      1.958 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 BOND  NR1  CPH1     400.000 {SD=      .027}      1.380 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 BOND  NR1  CPH2     400.000 {SD=      .027}      1.360 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 BOND  NR1  H        466.000 {SD=      .025}      1.000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 BOND  NR2  CPH1     400.000 {SD=      .027}      1.380 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 BOND  NR2  CPH2     400.000 {SD=      .027}      1.320 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 BOND  NR2  FE        65.000 {SD=      .067}      2.200 ! ALLOW HEM 
      ! Heme (6-liganded): His ligand (KK 05/13/91) 
 BOND  NR3  CPH1     380.000 {SD=      .028}      1.370 ! ALLOW ARO 
      ! his, adm jr., 6/28/90 
 BOND  NR3  CPH2     380.000 {SD=      .028}      1.320 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 BOND  NR3  H        453.000 {SD=      .026}      1.000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 BOND  NY   CA       270.000 {SD=      .033}      1.370 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 BOND  NY   CPT      270.000 {SD=      .033}      1.375 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 BOND  NY   H        465.000 {SD=      .025}       .976 ! ALLOW   ARO 
      ! indole JWK   08/28/89 
 BOND  O    C        620.000 {SD=      .022}      1.230 ! ALLOW   PEP POL ARO 
      ! Peptide geometry, condensed phase (LK) 
 BOND  O    CC       650.000 {SD=      .021}      1.230 ! ALLOW   PEP POL ARO 
      ! adm jr. 4/10/91, acetamide 
 BOND  OB   CC       750.000 {SD=      .020}      1.220 ! ALLOW   PEP POL ARO 
      ! adm jr., 10/17/90, acetic acid vibrations and geom. 
 BOND  OB   CD       750.000 {SD=      .020}      1.220 ! ALLOW   PEP POL ARO 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 BOND  OC   CA       525.000 {SD=      .024}      1.260 ! ALLOW   PEP POL ARO ION 
      ! adm jr. 8/27/91, phenoxide 
 BOND  OC   CC       525.000 {SD=      .024}      1.260 ! ALLOW   PEP POL ARO ION 
      ! adm jr. 7/23/91, acetic acid 
 BOND  OC   CT2      450.000 {SD=      .026}      1.330 ! ALLOW   ALC 
      ! ethoxide 6-31+G* geom/freq, adm jr., 6/1/92 
 BOND  OC   CT3      450.000 {SD=      .026}      1.330 ! ALLOW   ALC 
      ! methoxide 6-31+G* geom/freq, adm jr., 6/1/92 
 BOND  OH1  CA       334.300 {SD=      .030}      1.411 ! ALLOW   ARO ALC 
      ! MeOH, EMB 10/10/89, 
 BOND  OH1  CD       230.000 {SD=      .036}      1.400 ! ALLOW   PEP POL ARO ALC 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 BOND  OH1  CT1      428.000 {SD=      .026}      1.420 ! ALLOW   ALI ALC ARO 
      ! methanol vib fit EMB 11/21/89 
 BOND  OH1  CT2      428.000 {SD=      .026}      1.420 ! ALLOW   ALI ALC ARO 
      ! methanol vib fit EMB 11/21/89 
 BOND  OH1  CT3      428.000 {SD=      .026}      1.420 ! ALLOW   ALI ALC ARO 
      ! methanol vib fit EMB 11/21/89 
 BOND  OH1  H        545.000 {SD=      .023}       .960 ! ALLOW   ALC ARO 
      ! EMB 11/21/89 methanol vib fit 
 BOND  OM   CM      1115.000 {SD=      .016}      1.128 ! ALLOW HEM 
      ! Heme (6-liganded): CO ligand (KK 05/13/91) 
 BOND  OM   FE       250.000 {SD=      .034}      1.800 ! ALLOW HEM 
      ! Heme (6-liganded): O2 ligand (KK 05/13/91) 
 BOND  OM   OM       600.000 {SD=      .022}      1.230 ! ALLOW HEM 
      ! Heme (6-liganded): O2 ligand (KK 05/13/91) 
 BOND  OS   CD       150.000 {SD=      .044}      1.334 ! ALLOW POL PEP 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 BOND  OS   CT3      340.000 {SD=      .030}      1.430 ! ALLOW POL PEP 
      ! adm jr., 4/05/91, for PRES CT1 from methylacetate 
 BOND  OT   HT       450.000 {SD=      .026}       .957 ! ALLOW   WAT 
      ! FROM TIPS3P GEOM 
 BOND  S    CT2      198.000 {SD=      .039}      1.818 ! ALLOW   ALI SUL ION 
      ! fitted to C-S s   9/26/92 (FL) 
 BOND  S    CT3      240.000 {SD=      .035}      1.816 ! ALLOW   ALI SUL ION 
      ! fitted to C-S s   9/26/92 (FL) 
 BOND  S    HS       275.000 {SD=      .033}      1.325 ! ALLOW   SUL ION 
      ! methanethiol pure solvent, adm jr., 6/22/92 
 BOND  SM   CT2      214.000 {SD=      .037}      1.816 ! ALLOW   SUL ION 
      ! improved CSSC dihedral in DMDS  5/15/92 (FL) 
 BOND  SM   CT3      214.000 {SD=      .037}      1.816 ! ALLOW   SUL ION 
      ! improved CSSC dihedral in DMDS  5/15/92 (FL) 
 BOND  SM   SM       173.000 {SD=      .041}      2.029 ! ALLOW   SUL ION 
      ! improved CSSC dihedral in DMDS  5/15/92 (FL) 
 BOND  SS   CS       205.000 {SD=      .038}      1.836 ! ALLOW   SUL 
      ! methylthiolate 6-31+G* geom/freq, adm jr., 6/1/92 

 ANGLE  CA   CA   CA        40.00 {SD=      .086}    120.0000    UB     35.000       2.416 ! ALLOW   ARO 
      ! JES 8/25/89 
 ANGLE  CP1  N    C         60.00 {SD=      .070}    117.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  CP2  CP1  C         52.00 {SD=      .075}    112.3000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  CP2  CP1  CC        52.00 {SD=      .075}    112.3000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  CP2  CP1  CD        50.00 {SD=      .077}    112.3000 ! ALLOW PRO PEP 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  CP2  CP2  CP1       70.00 {SD=      .065}    108.5000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  CP3  CP2  CP2       70.00 {SD=      .065}    107.5000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  CP3  N    C         60.00 {SD=      .070}    117.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  CP3  N    CP1      100.00 {SD=      .054}    114.2000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  CP3  NP   CP1      100.00 {SD=      .054}    111.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  CPA  CPB  C         70.00 {SD=      .065}    126.7400 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 ANGLE  CPA  CPM  CPA       94.20 {SD=      .056}    125.1200 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 ANGLE  CPA  NPH  CPA      139.30 {SD=      .046}    103.9000 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 ANGLE  CPB  C    C         70.00 {SD=      .065}    121.5000 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 ANGLE  CPB  CPB  C         70.00 {SD=      .065}    126.7500 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 ANGLE  CPB  CPB  CPA       30.80 {SD=      .098}    106.5100 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 ANGLE  CPH2 NR1  CPH1     130.00 {SD=      .048}    107.5000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 ANGLE  CPH2 NR2  CPH1     130.00 {SD=      .048}    104.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 ANGLE  CPH2 NR3  CPH1     145.00 {SD=      .045}    108.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 ANGLE  CPM  CPA  CPB       61.60 {SD=      .069}    124.0700 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 ANGLE  CPT  CA   CA        60.00 {SD=      .070}    118.0000 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  CPT  CPT  CA        60.00 {SD=      .070}    122.0000 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  CPT  CY   CA       120.00 {SD=      .050}    107.4000    UB     25.000       2.261 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  CPT  NY   CA       110.00 {SD=      .052}    108.0000 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  CT1  CT1  C         52.00 {SD=      .075}    108.0000 ! ALLOW   ALI PEP POL ARO 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  CT1  CT1  CC        52.00 {SD=      .075}    108.0000 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 ANGLE  CT1  CT1  CT1       53.35 {SD=      .074}    111.0000    UB      8.000       2.561 ! ALLOW ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  CT1  CT2  CA        51.80 {SD=      .076}    107.5000 ! ALLOW   ALI ARO 
      ! PARALLH19 (JES) 
 ANGLE  CT1  CT2  CC        52.00 {SD=      .075}    108.0000 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 ANGLE  CT1  CT2  CD        52.00 {SD=      .075}    108.0000 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  CT1  CT2  CPH1      58.35 {SD=      .071}    113.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, from CT2CT2CT, U-B omitted 
 ANGLE  CT1  CT2  CT1       58.35 {SD=      .071}    113.5000    UB     11.160       2.561 ! ALLOW   ALI 
      ! alkane frequencies (MJF), alkane geometries (SF) 
 ANGLE  CT1  NH1  C         50.00 {SD=      .077}    120.0000 ! ALLOW   ALI PEP POL ARO 
      ! NMA Vib Modes (LK) 
 ANGLE  CT2  CA   CA        45.80 {SD=      .080}    122.3000 ! ALLOW   ALI ARO 
      ! PARALLH19 (JES) 
 ANGLE  CT2  CPB  CPA       65.00 {SD=      .067}    126.7400 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 ANGLE  CT2  CPB  CPB       65.00 {SD=      .067}    126.7500 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 ANGLE  CT2  CPH1 CPH1      45.80 {SD=      .080}    130.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FC=>CT2CA CA,BA=> CRYSTALS 
 ANGLE  CT2  CT1  C         52.00 {SD=      .075}    108.0000 ! ALLOW   ALI PEP POL ARO 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  CT2  CT1  CC        52.00 {SD=      .075}    108.0000 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 ANGLE  CT2  CT1  CD        52.00 {SD=      .075}    108.0000 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  CT2  CT1  CT1       53.35 {SD=      .074}    111.0000    UB      8.000       2.561 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  CT2  CT2  C         52.00 {SD=      .075}    108.0000 ! ALLOW   ALI PEP POL ARO 
      ! from CT2  CT1  C, for lactams, adm jr. 
 ANGLE  CT2  CT2  CC        52.00 {SD=      .075}    108.0000 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 ANGLE  CT2  CT2  CD        52.00 {SD=      .075}    108.0000 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  CT2  CT2  CPB       70.00 {SD=      .065}    113.0000 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 ANGLE  CT2  CT2  CT1       58.35 {SD=      .071}    113.5000    UB     11.160       2.561 ! ALLOW   ALI 
      ! alkane frequencies (MJF), alkane geometries (SF) 
 ANGLE  CT2  CT2  CT2       58.35 {SD=      .071}    113.6000    UB     11.160       2.561 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  CT2  CT3  CT1       58.35 {SD=      .071}    113.5000    UB     11.160       2.561 ! ALLOW   ALI 
      ! alkane frequencies (MJF), alkane geometries (SF) 
 ANGLE  CT2  CY   CA        45.80 {SD=      .080}    129.4000 ! ALLOW ARO 
      !adm jr., 5/08/91,  indole CCDB structure search 
 ANGLE  CT2  CY   CPT       45.80 {SD=      .080}    124.0000 ! ALLOW ARO 
      !adm jr., 5/08/91,  indole CCDB structure search 
 ANGLE  CT2  NC2  C         62.30 {SD=      .069}    120.0000 ! ALLOW   ALI POL PEP ARO 
      ! 107.5->120.0 to make planar Arg (KK) 
 ANGLE  CT2  NH1  C         50.00 {SD=      .077}    120.0000 ! ALLOW   ALI PEP POL ARO 
      ! NMA Vib Modes (LK) 
 ANGLE  CT2  OS   CD        40.00 {SD=      .086}    109.6000    UB     30.000       2.265 ! ALLOW  POL PEP 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  CT3  CA   CA        45.80 {SD=      .080}    122.3000 ! ALLOW   ALI ARO 
      ! toluene, adm jr., 3/7/92 
 ANGLE  CT3  CPB  CPA       65.00 {SD=      .067}    126.7400 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 ANGLE  CT3  CPB  CPB       65.00 {SD=      .067}    126.7500 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 ANGLE  CT3  CPH1 CPH1      45.80 {SD=      .080}    130.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FC=>CT2CA CA,BA=> CRYSTALS 
 ANGLE  CT3  CT1  C         52.00 {SD=      .075}    108.0000 ! ALLOW   ALI PEP POL ARO 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  CT3  CT1  CC        52.00 {SD=      .075}    108.0000 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 4/09/92, for ALA cter 
 ANGLE  CT3  CT1  CT1       53.35 {SD=      .074}    108.5000    UB      8.000       2.561 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  CT3  CT1  CT2       53.35 {SD=      .074}    113.0000    UB      8.000       2.561 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  CT3  CT1  CT3       53.35 {SD=      .074}    113.0000    UB      8.000       2.561 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  CT3  CT2  CA        51.80 {SD=      .076}    107.5000 ! ALLOW   ALI ARO 
      ! ethylbenzene, adm jr., 3/7/92 
 ANGLE  CT3  CT2  CPH1      58.35 {SD=      .071}    113.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, from CT2CT2CT, U-B omitted 
 ANGLE  CT3  CT2  CT1       58.35 {SD=      .071}    113.5000    UB     11.160       2.561 ! ALLOW   ALI 
      ! alkane frequencies (MJF), alkane geometries (SF) 
 ANGLE  CT3  CT2  CT2       58.00 {SD=      .071}    115.0000    UB      8.000       2.561 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  CT3  CT2  CT3       53.35 {SD=      .074}    114.0000    UB      8.000       2.561 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  CT3  NC2  C         62.30 {SD=      .069}    120.0000 ! ALLOW   ALI POL PEP ARO 
      ! methylguanidinium, adm jr., 3/26/92 
 ANGLE  CT3  NH1  C         50.00 {SD=      .077}    120.0000 ! ALLOW   ALI PEP POL ARO 
      ! NMA Vib Modes (LK) 
 ANGLE  CT3  OS   CD        40.00 {SD=      .086}    109.6000    UB     30.000       2.265 ! ALLOW  POL PEP 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  CT3  S    CT2       34.00 {SD=      .093}     95.0000 ! ALLOW   ALI SUL ION 
      ! expt. MeEtS,    3/26/92 (FL) 
 ANGLE  CY   CPT  CA       160.00 {SD=      .043}    130.6000 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  CY   CPT  CPT      110.00 {SD=      .052}    107.4000 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  CY   CT2  CT1       58.35 {SD=      .071}    114.0000 ! ALLOW ARO 
      ! from TRP crystal, JWK 
 ANGLE  CY   CT2  CT3       58.35 {SD=      .071}    114.0000 ! ALLOW ARO 
      ! from TRP crystal, JWK 
 ANGLE  FE   NPH  CPA       96.15 {SD=      .055}    128.0500 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 ANGLE  FE   NR2  CPH1      30.00 {SD=      .099}    133.0000 ! ALLOW HEM 
      ! Heme (6-liganded): ligand links (KK 05/13/91) 
 ANGLE  FE   NR2  CPH2      30.00 {SD=      .099}    123.0000 ! ALLOW HEM 
      ! Heme (6-liganded): ligand links (KK 05/13/91) 
 ANGLE  H    NH1  C         34.00 {SD=      .093}    123.0000 ! ALLOW   PEP POL ARO 
      ! NMA Vib Modes (LK) 
 ANGLE  H    NH1  CT1       35.00 {SD=      .092}    117.0000 ! ALLOW   PEP POL ARO ALI 
      ! NMA Vibrational Modes (LK) 
 ANGLE  H    NH1  CT2       35.00 {SD=      .092}    117.0000 ! ALLOW   PEP POL ARO ALI 
      ! NMA Vibrational Modes (LK) 
 ANGLE  H    NH1  CT3       35.00 {SD=      .092}    117.0000 ! ALLOW   PEP POL ARO ALI 
      ! NMA Vibrational Modes (LK) 
 ANGLE  H    NH2  CC        50.00 {SD=      .077}    120.0000 ! ALLOW   POL PEP ARO 
      ! his, adm jr. 8/13/90 acetamide geometry and vibrations 
 ANGLE  H    NH2  H         23.00 {SD=      .113}    120.0000 ! ALLOW   POL 
      ! adm jr. 8/13/90 acetamide geometry and vibrations 
 ANGLE  H    NR1  CPH1      30.00 {SD=      .099}    125.5000    UB     20.000       2.150 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 ANGLE  H    NR1  CPH2      30.00 {SD=      .099}    127.0000    UB     20.000       2.140 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 ANGLE  H    NR3  CPH1      25.00 {SD=      .109}    126.0000    UB     15.000       2.130 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 ANGLE  H    NR3  CPH2      25.00 {SD=      .109}    126.0000    UB     15.000       2.090 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 ANGLE  H    NY   CA        28.00 {SD=      .103}    126.0000 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  H    NY   CPT       28.00 {SD=      .103}    126.0000 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  H    OH1  CA        65.00 {SD=      .067}    108.0000 ! ALLOW   ALC ARO 
      ! JES 8/25/89 phenol 
 ANGLE  H    OH1  CD        55.00 {SD=      .073}    115.0000 ! ALLOW   ALC ARO PEP POL 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  H    OH1  CT1       57.50 {SD=      .072}    106.0000 ! ALLOW   ALC ARO ALI 
      ! methanol vib fit EMB 11/21/89 
 ANGLE  H    OH1  CT2       57.50 {SD=      .072}    106.0000 ! ALLOW   ALC ARO ALI 
      ! methanol vib fit EMB 11/21/89 
 ANGLE  H    OH1  CT3       57.50 {SD=      .072}    106.0000 ! ALLOW   ALC ARO ALI 
      ! methanol vib fit EMB 11/21/89 
 ANGLE  HA   C    C         50.00 {SD=      .077}    120.5000 ! ALLOW   PEP POL ARO 
      ! Heme vinyl substituent (KK from propene (JCS)) 
 ANGLE  HA   C    CPB       50.00 {SD=      .077}    120.0000 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 ANGLE  HA   C    HA        50.00 {SD=      .077}    118.0000 ! ALLOW   PEP POL ARO 
      ! Heme vinyl substituent (KK from propene (JCS)) 
 ANGLE  HA   CA   CA        29.00 {SD=      .101}    120.0000    UB     25.000       2.152 ! ALLOW ARO 
      ! trp, adm jr., 10/02/89 
 ANGLE  HA   CA   CPT       41.00 {SD=      .085}    122.0000 ! ALLOW   ARO 
      !adm jr., 5/08/91, indole CCDB structure search 
 ANGLE  HA   CA   CY        32.00 {SD=      .096}    125.0000    UB     25.000       2.173 ! ALLOW   ARO 
      ! JWK 05/14/91 new theta0 and r0UB from indole 
 ANGLE  HA   CP2  CP1       33.43 {SD=      .094}    110.1000    UB     22.530       2.179 ! ALLOW ALI PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HA   CP2  CP2       26.50 {SD=      .106}    110.1000    UB     22.530       2.179 ! ALLOW ALI PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HA   CP2  CP3       26.50 {SD=      .106}    110.1000    UB     22.530       2.179 ! ALLOW ALI PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HA   CP2  HA        35.50 {SD=      .091}    109.0000    UB      5.400       1.802 ! ALLOW ALI PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HA   CP3  CP2       26.50 {SD=      .106}    110.1000    UB     22.530       2.179 ! ALLOW ALI PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HA   CP3  HA        35.50 {SD=      .091}    109.0000    UB      5.400       1.802 ! ALLOW ALI PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HA   CPM  CPA       12.70 {SD=      .153}    117.4400 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 ANGLE  HA   CPM  FE          .00 {SD=999999.000}    180.0000 ! ALLOW HEM 
      ! Heme (6-liganded): for "ic para" only (KK 05/13/91) 
 ANGLE  HA   CS   CT3       34.60 {SD=      .093}    110.1000    UB     22.530       2.179 ! ALLOW SUL 
      ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92 
 ANGLE  HA   CS   HA        35.50 {SD=      .091}    108.4000    UB     14.000       1.775 ! ALLOW SUL 
      ! methylthiolate 6-31+G* geom/freq, adm jr., 6/1/92 
 ANGLE  HA   CT1  C         33.00 {SD=      .095}    109.5000    UB     30.000       2.163 ! ALLOW   ALI PEP POL ARO 
      ! alanine dipeptide, LK, replaced, adm jr., 5/09/91 
 ANGLE  HA   CT1  CD        33.00 {SD=      .095}    109.5000    UB     30.000       2.163 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  HA   CT1  CT1       34.50 {SD=      .093}    110.1000    UB     22.530       2.179 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  HA   CT1  CT2       34.50 {SD=      .093}    109.1000    UB     22.530       2.179 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  HA   CT1  CT3       34.50 {SD=      .093}    109.1000    UB     22.530       2.179 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  HA   CT1  HA        35.50 {SD=      .091}    109.0000    UB      5.400       1.802 ! TEST for test cpd 
      ! based on HA   CT2  HA 
 ANGLE  HA   CT2  C         33.00 {SD=      .095}    109.5000    UB     30.000       2.163 ! ALLOW   ALI PEP POL ARO 
      ! alanine dipeptide, LK, replaced, adm jr., 5/09/91 
 ANGLE  HA   CT2  CA        49.30 {SD=      .077}    107.5000 ! ALLOW   ALI ARO 
      ! PARALLH19 (JES) 
 ANGLE  HA   CT2  CC        33.00 {SD=      .095}    109.5000    UB     30.000       2.163 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 ANGLE  HA   CT2  CD        33.00 {SD=      .095}    109.5000    UB     30.000       2.163 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  HA   CT2  CPB       50.00 {SD=      .077}    109.5000 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 ANGLE  HA   CT2  CPH1      33.43 {SD=      .094}    109.5000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, from CT2CT2HA, U-B OMITTED 
 ANGLE  HA   CT2  CT1       33.43 {SD=      .094}    110.1000    UB     22.530       2.179 ! ALLOW   ALI 
      ! alkane frequencies (MJF), alkane geometries (SF) 
 ANGLE  HA   CT2  CT2       26.50 {SD=      .106}    110.1000    UB     22.530       2.179 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  HA   CT2  CT3       34.60 {SD=      .093}    110.1000    UB     22.530       2.179 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  HA   CT2  CY        33.43 {SD=      .094}    109.5000 ! ALLOW ARO 
      ! ADM JR., 10/02/89, from CT2CT2HA (U-B OMITTED), FOR JOANNA 
 ANGLE  HA   CT2  HA        35.50 {SD=      .091}    109.0000    UB      5.400       1.802 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  HA   CT3  C         33.00 {SD=      .095}    109.5000    UB     30.000       2.163 ! ALLOW   ALI PEP POL ARO 
      ! alanine dipeptide, LK, replaced, adm jr., 5/09/91 
 ANGLE  HA   CT3  CA        49.30 {SD=      .077}    107.5000 ! ALLOW   ALI ARO 
      ! toluene, adm jr. 3/7/92 
 ANGLE  HA   CT3  CC        33.00 {SD=      .095}    109.5000    UB     30.000       2.163 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 ANGLE  HA   CT3  CD        33.00 {SD=      .095}    109.5000    UB     30.000       2.163 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  HA   CT3  CPB       50.00 {SD=      .077}    109.5000 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 ANGLE  HA   CT3  CPH1      33.43 {SD=      .094}    109.5000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, from CT2CT2HA, U-B OMITTED 
 ANGLE  HA   CT3  CS        34.60 {SD=      .093}    110.1000    UB     22.530       2.179 ! ALLOW SUL 
      ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92 
 ANGLE  HA   CT3  CT1       33.43 {SD=      .094}    110.1000    UB     22.530       2.179 ! ALLOW   ALI 
      ! alkane frequencies (MJF), alkane geometries (SF) 
 ANGLE  HA   CT3  CT2       34.60 {SD=      .093}    110.1000    UB     22.530       2.179 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  HA   CT3  CT3       37.50 {SD=      .089}    110.1000    UB     22.530       2.179 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  HA   CT3  HA        35.50 {SD=      .091}    108.4000    UB      5.400       1.802 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 ANGLE  HA   CY   CA        20.00 {SD=      .122}    126.4000    UB     25.000       2.186 ! ALLOW   ARO 
      ! JWK 05/14/91 new theta0 and r0UB from indole 
 ANGLE  HA   CY   CPT       32.00 {SD=      .096}    126.4000    UB     25.000       2.255 ! ALLOW   ARO 
      ! JWK 05/14/91 new theta0 and r0UB from indole 
 ANGLE  HB   CP1  C         50.00 {SD=      .077}    111.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HB   CP1  CC        50.00 {SD=      .077}    112.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HB   CP1  CD        50.00 {SD=      .077}    112.0000 ! ALLOW PEP POL PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HB   CP1  CP2       35.00 {SD=      .092}    115.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HB   CT1  C         50.00 {SD=      .077}    109.5000 ! ALLOW  PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  HB   CT1  CC        50.00 {SD=      .077}    109.5000 ! ALLOW  PEP POL 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 ANGLE  HB   CT1  CD        50.00 {SD=      .077}    109.5000 ! ALLOW  PEP POL 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  HB   CT1  CT1       35.00 {SD=      .092}    111.0000 ! ALLOW  PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  HB   CT1  CT2       35.00 {SD=      .092}    111.0000 ! ALLOW  PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  HB   CT1  CT3       35.00 {SD=      .092}    111.0000 ! ALLOW  PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  HB   CT2  C         50.00 {SD=      .077}    109.5000 ! ALLOW  PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  HB   CT2  CC        50.00 {SD=      .077}    109.5000 ! ALLOW  PEP POL 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 ANGLE  HB   CT2  CD        50.00 {SD=      .077}    109.5000 ! ALLOW  PEP POL 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  HB   CT2  HB        36.00 {SD=      .091}    115.0000 ! ALLOW   PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  HB   CT3  C         50.00 {SD=      .077}    109.5000 ! ALLOW  PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  HC   NC2  C         49.00 {SD=      .078}    120.0000 ! ALLOW   POL PEP ARO 
      ! 35.3->49.0 GUANIDINIUM (KK) 
 ANGLE  HC   NC2  CT2       40.40 {SD=      .086}    120.0000 ! ALLOW   POL ALI 
      ! 107.5->120.0 to make planar Arg (KK) 
 ANGLE  HC   NC2  CT3       40.40 {SD=      .086}    120.0000 ! ALLOW   POL ALI 
      ! methylguanidinium, adm jr., 3/26/92 
 ANGLE  HC   NC2  HC        25.00 {SD=      .109}    120.0000 ! ALLOW   POL 
      ! 40.0->25.0 GUANIDINIUM (KK) 
 ANGLE  HC   NH2  CT3       50.00 {SD=      .077}    111.0000 ! ALLOW   POL 
      ! methylamine geom/freq, adm jr., 6/2/92 
 ANGLE  HC   NH2  HC        39.00 {SD=      .087}    106.5000 ! ALLOW   POL 
      ! 40.0->25.0 GUANIDINIUM (KK) 
 ANGLE  HC   NH3  CT1       30.00 {SD=      .099}    109.5000    UB     20.000       2.074 ! ALLOW   POL ALI 
      ! new stretch and bend; methylammonium (KK 03/10/92) 
 ANGLE  HC   NH3  CT2       30.00 {SD=      .099}    109.5000    UB     20.000       2.074 ! ALLOW   POL ALI 
      ! new stretch and bend; methylammonium (KK 03/10/92) 
 ANGLE  HC   NH3  CT3       30.00 {SD=      .099}    109.5000    UB     20.000       2.074 ! ALLOW   POL ALI 
      ! new stretch and bend; methylammonium (KK 03/10/92) 
 ANGLE  HC   NH3  HC        44.00 {SD=      .082}    109.5000 ! ALLOW   POL 
      ! new stretch and bend; methylammonium (KK 03/10/92) 
 ANGLE  HC   NP   CP1       33.00 {SD=      .095}    109.5000    UB      4.000       2.056 ! ALLOW POL ALI PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HC   NP   CP3       33.00 {SD=      .095}    109.5000    UB      4.000       2.056 ! ALLOW POL ALI PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HC   NP   HC        51.00 {SD=      .076}    107.5000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  HP   CA   CA        30.00 {SD=      .099}    120.0000    UB     22.000       2.152 ! ALLOW   ARO 
      ! JES 8/25/89 benzene 
 ANGLE  HP   CA   CPT       30.00 {SD=      .099}    122.0000    UB     22.000       2.146 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  HP   CA   CY        32.00 {SD=      .096}    125.0000    UB     25.000       2.173 ! ALLOW   ARO 
      ! JWK 05/14/91 new theta0 and r0UB from indole 
 ANGLE  HP   CY   CA        32.00 {SD=      .096}    126.4000    UB     25.000       2.186 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  HP   CY   CPT       32.00 {SD=      .096}    126.4000    UB     25.000       2.255 ! ALLOW   ARO 
      ! JWK 05/14/91 new theta0 and r0UB from indole 
 ANGLE  HR1  CPH1 CPH1      22.00 {SD=      .116}    130.0000    UB     15.000       2.215 ! ALLOW ARO 
      ! adm jr., 6/27/90, his 
 ANGLE  HR3  CPH1 CPH1      25.00 {SD=      .109}    130.0000    UB     20.000       2.200 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 ANGLE  HS   S    CT2       38.80 {SD=      .087}     95.0000 ! ALLOW   SUL ION ALI 
      ! methanethiol pure solvent, adm jr., 6/22/92 
 ANGLE  HS   S    CT3       43.00 {SD=      .083}     95.0000 ! ALLOW   SUL ION ALI 
      ! methanethiol pure solvent, adm jr., 6/22/92 
 ANGLE  HT   OT   HT        55.00 {SD=      .073}    104.5200 ! ALLOW WAT 
      ! TIP3P GEOMETRY, ADM JR. 
 ANGLE  N    C    CP1       20.00 {SD=      .122}    111.5000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  N    C    CT1       20.00 {SD=      .122}    114.0000 ! ALLOW ALI PEP POL ARO PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  N    C    CT2       20.00 {SD=      .122}    112.5000 ! ALLOW ALI PEP POL ARO PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  N    C    CT3       20.00 {SD=      .122}    112.5000 ! ALLOW ALI PEP POL ARO PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  N    CP1  C         50.00 {SD=      .077}    107.2000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  N    CP1  CC        50.00 {SD=      .077}    108.2000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  N    CP1  CD        50.00 {SD=      .077}    108.2000 ! ALLOW PRO PEP 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  N    CP1  CP2       70.00 {SD=      .065}    107.8000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  N    CP1  HB        48.00 {SD=      .079}    111.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  N    CP3  CP2       70.00 {SD=      .065}    108.5000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  N    CP3  HA        48.00 {SD=      .079}    108.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  NC2  C    NC2       52.00 {SD=      .075}    120.0000    UB     90.000       2.364 ! ALLOW   POL PEP ARO 
      ! changed from 60.0/120.3 for guanidinium (KK) 
 ANGLE  NC2  CT2  CT2       67.70 {SD=      .066}    107.5000 ! ALLOW   ALI POL 
      ! arg, (DS) 
 ANGLE  NC2  CT2  HA        51.50 {SD=      .076}    107.5000 ! ALLOW   ALI POL 
      ! arg, (DS) 
 ANGLE  NC2  CT3  HA        51.50 {SD=      .076}    107.5000 ! ALLOW   ALI POL 
      ! methylguanidinium, adm jr., 3/26/92 
 ANGLE  NH1  C    CP1       80.00 {SD=      .061}    116.5000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  NH1  C    CT1       80.00 {SD=      .061}    116.5000 ! ALLOW   ALI PEP POL ARO 
      ! NMA Vib Modes (LK) 
 ANGLE  NH1  C    CT2       80.00 {SD=      .061}    116.5000 ! ALLOW   ALI PEP POL ARO 
      ! NMA Vib Modes (LK) 
 ANGLE  NH1  C    CT3       80.00 {SD=      .061}    116.5000 ! ALLOW   ALI PEP POL ARO 
      ! NMA Vib Modes (LK) 
 ANGLE  NH1  CT1  C         50.00 {SD=      .077}    107.0000 ! ALLOW   PEP POL ARO ALI 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  NH1  CT1  CC        50.00 {SD=      .077}    107.0000 ! ALLOW   PEP POL ARO ALI 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 ANGLE  NH1  CT1  CD        50.00 {SD=      .077}    107.0000 ! ALLOW   PEP POL ARO ALI 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  NH1  CT1  CT1       70.00 {SD=      .065}    113.5000 ! ALLOW   ALI PEP POL ARO 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  NH1  CT1  CT2       70.00 {SD=      .065}    113.5000 ! ALLOW   ALI PEP POL ARO 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  NH1  CT1  CT3       70.00 {SD=      .065}    113.5000 ! ALLOW   ALI PEP POL ARO 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  NH1  CT1  HB        48.00 {SD=      .079}    108.0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  NH1  CT2  C         50.00 {SD=      .077}    107.0000 ! ALLOW   PEP POL ARO ALI 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  NH1  CT2  CC        50.00 {SD=      .077}    107.0000 ! ALLOW   PEP POL ARO ALI 
      ! adm jr. 5/20/92, for asn,asp,gln,glu and cters 
 ANGLE  NH1  CT2  CD        50.00 {SD=      .077}    107.0000 ! ALLOW   PEP POL ARO ALI 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  NH1  CT2  CT2       70.00 {SD=      .065}    113.5000 ! ALLOW   ALI PEP POL ARO 
      ! from NH1  CT1  CT2, for lactams, adm jr. 
 ANGLE  NH1  CT2  HA        51.50 {SD=      .076}    109.5000 ! ALLOW   ALI PEP POL ARO 
      ! from NH1  CT3  HA, for lactams, adm jr. 
 ANGLE  NH1  CT2  HB        48.00 {SD=      .079}    108.0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  NH1  CT3  HA        51.50 {SD=      .076}    109.5000 ! ALLOW   ALI PEP POL ARO 
      ! NMA crystal (JCS) 
 ANGLE  NH2  CC   CP1       80.00 {SD=      .061}    112.5000 ! ALLOW ALI PEP POL ARO PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  NH2  CC   CT1       50.00 {SD=      .077}    116.5000    UB     50.000       2.450 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 8/13/90 acetamide geometry and vibrations 
 ANGLE  NH2  CC   CT2       50.00 {SD=      .077}    116.5000    UB     50.000       2.450 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 8/13/90 acetamide geometry and vibrations 
 ANGLE  NH2  CC   CT3       50.00 {SD=      .077}    116.5000    UB     50.000       2.450 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 8/13/90 acetamide geometry and vibrations 
 ANGLE  NH2  CC   HA        44.00 {SD=      .082}    111.0000    UB     50.000       1.980 ! ALLOW POL 
      ! adm jr., 5/13/91, formamide geometry and vibrations 
 ANGLE  NH2  CT3  HA        38.00 {SD=      .088}    109.5000    UB     50.000       2.140 ! ALLOW POL 
      ! methylamine geom/freq, adm jr., 6/2/92 
 ANGLE  NH3  CT1  C         43.70 {SD=      .082}    110.0000 ! ALLOW   PEP POL ARO ALI 
      ! new aliphatics, adm jr., 2/3/92 
 ANGLE  NH3  CT1  CC        43.70 {SD=      .082}    110.0000 ! ALLOW   PEP POL ARO ALI 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 ANGLE  NH3  CT1  CT1       67.70 {SD=      .066}    110.0000 ! ALLOW   ALI POL 
      ! new aliphatics, adm jr., 2/3/92 
 ANGLE  NH3  CT1  CT2       67.70 {SD=      .066}    110.0000 ! ALLOW   ALI POL 
      ! new aliphatics, adm jr., 2/3/92 
 ANGLE  NH3  CT1  CT3       67.70 {SD=      .066}    110.0000 ! ALLOW   ALI POL 
      ! new aliphatics, adm jr., 2/3/92 
 ANGLE  NH3  CT1  HB        51.50 {SD=      .076}    107.5000 ! ALLOW   ALI POL PEP 
      ! new aliphatics, adm jr., 2/3/92 
 ANGLE  NH3  CT2  C         43.70 {SD=      .082}    110.0000 ! ALLOW   PEP POL ARO ALI 
      ! alanine (JCS) 
 ANGLE  NH3  CT2  CC        43.70 {SD=      .082}    110.0000 ! ALLOW   PEP POL ARO ALI 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 ANGLE  NH3  CT2  CD        43.70 {SD=      .082}    110.0000 ! ALLOW   PEP POL ARO ALI 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  NH3  CT2  CT2       67.70 {SD=      .066}    110.0000 ! ALLOW   ALI POL 
      ! alanine (JCS) 
 ANGLE  NH3  CT2  HA        45.00 {SD=      .081}    107.5000    UB     35.000       2.101 ! ALLOW   ALI POL 
      ! new stretch and bend; methylammonium (KK 03/10/92) 
 ANGLE  NH3  CT2  HB        51.50 {SD=      .076}    107.5000 ! ALLOW   ALI POL PEP 
      ! for use on NTER -- from NH3 CT2HA (JCS) -- (LK) 
 ANGLE  NH3  CT3  HA        45.00 {SD=      .081}    107.5000    UB     35.000       2.101 ! ALLOW   ALI POL 
      ! new stretch and bend; methylammonium (KK 03/10/92) 
 ANGLE  NP   CP1  C         50.00 {SD=      .077}    106.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  NP   CP1  CC        50.00 {SD=      .077}    106.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  NP   CP1  CD        50.00 {SD=      .077}    106.0000 ! ALLOW PRO PEP 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  NP   CP1  CP2       70.00 {SD=      .065}    108.5000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  NP   CP1  HB        51.50 {SD=      .076}    107.5000 ! ALLOW ALI POL PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  NP   CP3  CP2       70.00 {SD=      .065}    108.5000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  NP   CP3  HA        51.50 {SD=      .076}    109.1500 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  NPH  CPA  CPB      122.00 {SD=      .049}    111.5400 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 ANGLE  NPH  CPA  CPM       88.00 {SD=      .058}    124.3900 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 ANGLE  NPH  FE   CM        50.00 {SD=      .077}     90.0000 ! ALLOW HEM 
      ! Heme (6-liganded): ligand links (KK 05/13/91) 
 ANGLE  NPH  FE   CPM         .00 {SD=999999.000}     45.0000 ! ALLOW HEM 
      ! Heme (6-liganded): for "ic para" only  (KK 05/13/91) 
 ANGLE  NPH  FE   NPH       14.39 {SD=      .143}     90.0000 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 ANGLE  NR1  CPH1 CPH1     130.00 {SD=      .048}    106.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 ANGLE  NR1  CPH1 CT2       45.80 {SD=      .080}    124.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FC FROM CA CT2CT 
 ANGLE  NR1  CPH1 CT3       45.80 {SD=      .080}    124.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FC FROM CA CT2CT 
 ANGLE  NR1  CPH1 HR3       25.00 {SD=      .109}    124.0000    UB     20.000       2.140 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 ANGLE  NR1  CPH2 HR1       25.00 {SD=      .109}    122.5000    UB     20.000       2.140 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 ANGLE  NR2  CPH1 CPH1     130.00 {SD=      .048}    110.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 ANGLE  NR2  CPH1 CT2       45.80 {SD=      .080}    120.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FC FROM CA CT2CT 
 ANGLE  NR2  CPH1 HR3       25.00 {SD=      .109}    120.0000    UB     20.000       2.140 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 ANGLE  NR2  CPH2 HR1       25.00 {SD=      .109}    125.0000    UB     20.000       2.120 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 ANGLE  NR2  CPH2 NR1      130.00 {SD=      .048}    112.5000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 ANGLE  NR2  FE   CM        50.00 {SD=      .077}    180.0000 ! ALLOW HEM 
      ! Heme (6-liganded): ligand links (KK 05/13/91) 
 ANGLE  NR2  FE   NPH       50.00 {SD=      .077}     90.0000 ! ALLOW HEM 
      ! Heme (6-liganded): ligand links (KK 05/13/91) 
 ANGLE  NR3  CPH1 CPH1     145.00 {SD=      .045}    108.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 ANGLE  NR3  CPH1 CT2       45.80 {SD=      .080}    122.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FC FROM CA CT2CT 
 ANGLE  NR3  CPH1 HR1       22.00 {SD=      .116}    122.0000    UB     15.000       2.180 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 ANGLE  NR3  CPH2 HR2       32.00 {SD=      .096}    126.0000    UB     25.000       2.140 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 ANGLE  NR3  CPH2 NR3      145.00 {SD=      .045}    108.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 ANGLE  NY   CA   CY       120.00 {SD=      .050}    110.0000    UB     25.000       2.240 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  NY   CA   HA        32.00 {SD=      .096}    125.0000    UB     25.000       2.177 ! ALLOW   ARO 
      ! JWK 05/14/91 new theta0 and r0UB from indole 
 ANGLE  NY   CA   HP        32.00 {SD=      .096}    125.0000    UB     25.000       2.177 ! ALLOW   ARO 
      ! JWK 05/14/91 new theta0 and r0UB from indole 
 ANGLE  NY   CPT  CA       160.00 {SD=      .043}    130.6000 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  NY   CPT  CPT      110.00 {SD=      .052}    107.4000 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 ANGLE  O    C    CP1       80.00 {SD=      .061}    118.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  O    C    CT1       80.00 {SD=      .061}    121.0000 ! ALLOW   ALI PEP POL ARO 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  O    C    CT2       80.00 {SD=      .061}    121.0000 ! ALLOW   ALI PEP POL ARO 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  O    C    CT3       80.00 {SD=      .061}    121.0000 ! ALLOW   ALI PEP POL ARO 
      ! Alanine Dipeptide ab initio calc's (LK) 
 ANGLE  O    C    H         50.00 {SD=      .077}    121.7000 ! ALLOW   PEP POL ARO 
      ! acetaldehyde (JCS) 
 ANGLE  O    C    N         80.00 {SD=      .061}    122.5000 ! ALLOW PRO PEP POL ARO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  O    C    NH1       80.00 {SD=      .061}    122.5000 ! ALLOW   PEP POL ARO 
      ! NMA Vib Modes (LK) 
 ANGLE  O    CC   CP1       80.00 {SD=      .061}    118.0000 ! ALLOW ALI PEP POL ARO PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  O    CC   CT1       15.00 {SD=      .140}    121.0000    UB     50.000       2.440 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 4/10/91, acetamide update 
 ANGLE  O    CC   CT2       15.00 {SD=      .140}    121.0000    UB     50.000       2.440 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 4/10/91, acetamide update 
 ANGLE  O    CC   CT3       15.00 {SD=      .140}    121.0000    UB     50.000       2.440 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 4/10/91, acetamide update 
 ANGLE  O    CC   HA        44.00 {SD=      .082}    122.0000 ! ALLOW POL 
      ! adm jr., 5/13/91, formamide geometry and vibrations 
 ANGLE  O    CC   NH2       75.00 {SD=      .063}    122.5000    UB     50.000       2.370 ! ALLOW   POL PEP ARO 
      ! adm jr. 4/10/91, acetamide update 
 ANGLE  OB   CD   CP1       70.00 {SD=      .065}    125.0000    UB     20.000       2.442 ! ALLOW ALI PEP POL ARO PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  OB   CD   CT1       70.00 {SD=      .065}    125.0000    UB     20.000       2.442 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  OB   CD   CT2       70.00 {SD=      .065}    125.0000    UB     20.000       2.442 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  OB   CD   CT3       70.00 {SD=      .065}    125.0000    UB     20.000       2.442 ! ALLOW   ALI PEP POL ARO 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 ANGLE  OC   CA   CA        40.00 {SD=      .086}    120.0000 ! ALLOW  POL ARO 
      ! adm jr. 8/27/91, phenoxide 
 ANGLE  OC   CC   CP1       40.00 {SD=      .086}    118.0000    UB     50.000       2.388 ! ALLOW ALI PEP POL ARO ION PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  OC   CC   CT1       40.00 {SD=      .086}    118.0000    UB     50.000       2.388 ! ALLOW   ALI PEP POL ARO ION 
      ! adm jr. 7/23/91, correction, ACETATE (KK) 
 ANGLE  OC   CC   CT2       40.00 {SD=      .086}    118.0000    UB     50.000       2.388 ! ALLOW   ALI PEP POL ARO ION 
      ! adm jr. 7/23/91, correction, ACETATE (KK) 
 ANGLE  OC   CC   CT3       40.00 {SD=      .086}    118.0000    UB     50.000       2.388 ! ALLOW   ALI PEP POL ARO ION 
      ! adm jr. 7/23/91, correction, ACETATE (KK) 
 ANGLE  OC   CC   OC       100.00 {SD=      .054}    124.0000    UB     70.000       2.225 ! ALLOW   POL ION PEP ARO 
      ! adm jr. 7/23/91, correction, ACETATE (KK) 
 ANGLE  OC   CT2  CT3       65.00 {SD=      .067}    122.0000 ! ALLOW  ALC 
      ! ethoxide 6-31+G* geom/freq, adm jr., 6/1/92 
 ANGLE  OC   CT2  HA        65.00 {SD=      .067}    118.3000 ! ALLOW  ALC 
      ! ethoxide 6-31+G* geom/freq, adm jr., 6/1/92 
 ANGLE  OC   CT3  HA        65.00 {SD=      .067}    118.3000 ! ALLOW  ALC 
      ! methoxide 6-31+G* geom/freq, adm jr., 6/1/92 
 ANGLE  OH1  CA   CA        45.20 {SD=      .081}    120.0000 ! ALLOW   ARO ALC 
      ! PARALLH19 WITH [122.3] (JES) 
 ANGLE  OH1  CD   CT2       55.00 {SD=      .073}    110.5000 ! ALLOW   ALI PEP POL ARO ALC 
      ! adm jr, 10/17/90, acetic acid vibrations 
 ANGLE  OH1  CD   CT3       55.00 {SD=      .073}    110.5000 ! ALLOW   ALI PEP POL ARO ALC 
      ! adm jr, 10/17/90, acetic acid vibrations 
 ANGLE  OH1  CD   OB        50.00 {SD=      .077}    123.0000    UB    210.000       2.262 ! ALLOW   PEP POL ARO ALC 
      ! adm jr, 10/17/90, acetic acid vibrations 
 ANGLE  OH1  CT1  CT1       75.70 {SD=      .063}    110.1000 ! ALLOW   ALI ALC ARO 
      ! MeOH, EMB, 10/10/89 
 ANGLE  OH1  CT1  CT3       75.70 {SD=      .063}    110.1000 ! ALLOW   ALI ALC ARO 
      ! MeOH, EMB, 10/10/89 
 ANGLE  OH1  CT1  HA        45.90 {SD=      .080}    108.8900 ! ALLOW   ALI ALC ARO 
      ! MeOH, EMB, 10/10/89 
 ANGLE  OH1  CT2  CT1       75.70 {SD=      .063}    110.1000 ! ALLOW   ALI ALC ARO 
      ! MeOH, EMB, 10/10/89 
 ANGLE  OH1  CT2  CT3       75.70 {SD=      .063}    110.1000 ! ALLOW   ALI ALC ARO 
      ! MeOH, EMB, 10/10/89 
 ANGLE  OH1  CT2  HA        45.90 {SD=      .080}    108.8900 ! ALLOW   ALI ALC ARO 
      ! MeOH, EMB, 10/10/89 
 ANGLE  OH1  CT3  HA        45.90 {SD=      .080}    108.8900 ! ALLOW   ALI ALC ARO 
      ! MeOH, EMB, 10/10/89 
 ANGLE  OM   CM   FE        35.00 {SD=      .092}    180.0000 ! ALLOW HEM 
      ! Heme (6-liganded): ligand links (KK 05/13/91) 
 ANGLE  OM   FE   NPH        5.00 {SD=      .243}     90.0000 ! ALLOW HEM 
      ! Heme (6-liganded): ligand links (KK 05/13/91) 
 ANGLE  OM   OM   FE          .00 {SD=999999.000}    180.0000 ! ALLOW HEM 
      ! Heme (6-liganded): ligand links (KK 05/13/91) 
 ANGLE  OS   CD   CP1       55.00 {SD=      .073}    109.0000    UB     20.000       2.326 ! ALLOW POL PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 ANGLE  OS   CD   CT1       55.00 {SD=      .073}    109.0000    UB     20.000       2.326 ! ALLOW POL PEP 
      ! adm jr., 4/05/91, for PRES CT1 from methylacetate 
 ANGLE  OS   CD   CT2       55.00 {SD=      .073}    109.0000    UB     20.000       2.326 ! ALLOW POL PEP 
      ! adm jr., 4/05/91, for PRES CT1 from methylacetate 
 ANGLE  OS   CD   CT3       55.00 {SD=      .073}    109.0000    UB     20.000       2.326 ! ALLOW POL PEP 
      ! adm jr., 4/05/91, for PRES CT1 from methylacetate 
 ANGLE  OS   CD   OB        90.00 {SD=      .057}    125.9000    UB    160.000       2.258 ! ALLOW  PEP POL 
      ! adm jr. 3/19/92, from lipid methyl acetate 
 ANGLE  OS   CT2  HA        60.00 {SD=      .070}    109.5000 ! ALLOW PEP POL 
      ! adm jr. 4/05/91, for PRES CT1 from methyl acetate 
 ANGLE  OS   CT3  HA        60.00 {SD=      .070}    109.5000 ! ALLOW PEP POL 
      ! adm jr. 4/05/91, for PRES CT1 from methyl acetate 
 ANGLE  S    CT2  CT1       58.00 {SD=      .071}    112.5000 ! ALLOW   ALI SUL ION 
      ! as in expt.MeEtS & DALC crystal,  5/15/92 
 ANGLE  S    CT2  CT2       58.00 {SD=      .071}    114.5000 ! ALLOW   ALI SUL ION 
      ! expt. MeEtS,     3/26/92 (FL) 
 ANGLE  S    CT2  CT3       58.00 {SD=      .071}    114.5000 ! ALLOW   ALI SUL ION 
      ! expt. MeEtS,     3/26/92 (FL) 
 ANGLE  S    CT2  HA        46.10 {SD=      .080}    111.3000 ! ALLOW   ALI SUL ION 
      ! vib. freq. and HF/6-31G* geo. (DTN) 8/24/90 
 ANGLE  S    CT3  HA        46.10 {SD=      .080}    111.3000 ! ALLOW   ALI SUL ION 
      ! vib. freq. and HF/6-31G* geo. (DTN) 8/24/90 
 ANGLE  SM   CT2  CT1       58.00 {SD=      .071}    112.5000 ! ALLOW   ALI SUL ION 
      ! as in expt.MeEtS & DALC crystal,  5/15/92 
 ANGLE  SM   CT2  HA        38.00 {SD=      .088}    111.0000 ! ALLOW   ALI SUL ION 
      ! new S-S atom type 8/24/90 
 ANGLE  SM   CT3  HA        38.00 {SD=      .088}    111.0000 ! ALLOW   ALI SUL ION 
      ! new S-S atom type 8/24/90 
 ANGLE  SM   SM   CT2       72.50 {SD=      .064}    103.3000 ! ALLOW   ALI SUL ION 
      ! expt. dimethyldisulfide,    3/26/92 (FL) 
 ANGLE  SM   SM   CT3       72.50 {SD=      .064}    103.3000 ! ALLOW   ALI SUL ION 
      ! expt. dimethyldisulfide,    3/26/92 (FL) 
 ANGLE  SS   CS   CT3       55.00 {SD=      .073}    118.0000 ! ALLOW SUL 
      ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92 
 ANGLE  SS   CS   HA        40.00 {SD=      .086}    112.3000 ! ALLOW SUL 
      ! methylthiolate 6-31+G* geom/freq, adm jr., 6/1/92 

!! Added for Gly terminal YT
 ANGLE  HB   CT1  HB        35.50 {SD=      .091}    109.0000    UB      5.400       1.802 ! TEST for test cpd 
      ! based on HA   CT2  HA 


 DIHEDRAL  C    CT1  NH1  C            .20 {SD=     1.217}    1    180.0000 ! ALLOW PEP 
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c 
 DIHEDRAL  C    CT2  NH1  C            .20 {SD=     1.217}    1    180.0000 ! ALLOW PEP 
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c 
 DIHEDRAL  C    N    CP1  C            .80 {SD=      .608}    3       .0000 ! ALLOW PRO PEP 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CA   CA   CA   CA          3.10 {SD=      .309}    2    180.0000 ! ALLOW   ARO 
      ! JES 8/25/89 
 DIHEDRAL  CA   CPT  CPT  CA          3.10 {SD=      .309}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  CA   CT2  CT1  C            .04 {SD=     2.721}    3       .0000 ! ALLOW   ARO 
      ! 2.7 kcal/mole CH3 rot in ethylbenzene, adm jr, 3/7/92 
 DIHEDRAL  CA   CY   CPT  CA          3.00 {SD=      .314}    2    180.0000 ! ALLOW   ARO 
      ! JWK 09/05/89 
 DIHEDRAL  CA   NY   CPT  CA          3.00 {SD=      .314}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  CC   CP1  N    C            .80 {SD=      .608}    3       .0000 ! ALLOW PRO PEP 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CC   CT1  CT2  CA           .04 {SD=     2.721}    3       .0000 ! ALLOW   ARO 
      ! 2.7 kcal/mole CH3 rot in ethylbenzene, adm jr, 3/7/92 
 DIHEDRAL  CC   CT1  NH1  C            .20 {SD=     1.217}    1    180.0000 ! ALLOW PEP POL 
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c 
 DIHEDRAL  CC   CT2  NH1  C            .20 {SD=     1.217}    1    180.0000 ! ALLOW PEP POL 
      ! Alanine dipeptide; NMA; acetate; etc. adm jr., 3/3/93c 
 DIHEDRAL  CD   CP1  N    C            .00 {SD=999999.000}    1    180.0000 ! ALLOW PRO PEP 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CD   CT1  NH1  C            .20 {SD=     1.217}    1    180.0000 ! ALLOW PEP POL 
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c 
 DIHEDRAL  CD   CT2  NH1  C            .20 {SD=     1.217}    1    180.0000 ! ALLOW PEP POL 
      ! Alanine dipeptide; NMA; acetate; etc. backbon adm jr., 3/3/93c 
 DIHEDRAL  CP1  C    N    CP1    MULTIPLE= 2       2.75     2    180.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .30     4       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CP2  CP1  N    C            .80 {SD=      .608}    3       .0000 ! ALLOW PRO PEP 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CP2  CP3  N    C            .00 {SD=999999.000}    3    180.0000 ! ALLOW PRO PEP 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CP2  CP3  N    CP1          .10 {SD=     1.721}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CP2  CP3  NP   CP1          .08 {SD=     1.924}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CP3  N    C    CP1    MULTIPLE= 2       2.75     2    180.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .30     4       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CP3  N    CP1  C            .10 {SD=     1.721}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CP3  N    CP1  CC           .10 {SD=     1.721}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CP3  N    CP1  CP2          .10 {SD=     1.721}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CP3  NP   CP1  C            .08 {SD=     1.924}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CP3  NP   CP1  CC           .08 {SD=     1.924}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CP3  NP   CP1  CD           .08 {SD=     1.924}    3       .0000 ! ALLOW PRO PEP 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CP3  NP   CP1  CP2          .08 {SD=     1.924}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CPH2 NR1  CPH1 CPH1       14.00 {SD=      .145}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 DIHEDRAL  CPH2 NR2  CPH1 CPH1       14.00 {SD=      .145}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 DIHEDRAL  CPH2 NR3  CPH1 CPH1       12.00 {SD=      .157}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 DIHEDRAL  CPT  CA   CA   CA          3.10 {SD=      .309}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  CPT  CPT  CA   CA          3.10 {SD=      .309}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  CPT  CPT  CY   CA          4.00 {SD=      .272}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  CPT  CPT  NY   CA          5.00 {SD=      .243}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  CT1  C    N    CP1    MULTIPLE= 2       2.75     2    180.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .30     4       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CT1  C    N    CP3    MULTIPLE= 2       2.75     2    180.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .30     4       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CT1  C    NH1  CT1    MULTIPLE= 2       1.60     1       .0000 !  ALLOW PEP 
      ! Revised to adjust NMA cis/trans energy difference. (LK) 
                                                   2.50     2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  CT1  CT1  NH1  C           1.80 {SD=      .406}    1       .0000 ! ALLOW PEP 
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c 
 DIHEDRAL  CT1  CT2  CA   CA           .23 {SD=     1.135}    2    180.0000 ! ALLOW   ARO ALI 
      ! ethylbenzene ethyl rotation, adm jr. 3/7/92 
 DIHEDRAL  CT1  CT2  CPH1 CPH1   MULTIPLE= 3        .20     1       .0000 ! ALLOW ARO 
      ! 4-ethylimidazole 4-21G rot bar, adm jr. 3/4/92 
                                                    .27     2       .0000 ! ALLOW ARO 
      ! 4-ethylimidazole 4-21G rot bar, adm jr. 3/4/92 
                                                    .00     3       .0000 ! ALLOW ARO 
      ! 4-ethylimidazole 4-21G rot bar, adm jr. 3/4/92 
 DIHEDRAL  CT1  CT2  CY   CA           .23 {SD=     1.135}    2    180.0000 ! ALLOW ARO 
      ! from ethylbenzene, adm jr., 3/7/92 
 DIHEDRAL  CT1  CT2  CY   CPT          .23 {SD=     1.135}    2    180.0000 ! ALLOW ARO 
      ! from ethylbenzene, adm jr., 3/7/92 
 DIHEDRAL  CT1  NH1  C    CP1    MULTIPLE= 2       1.60     1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                   2.50     2    180.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CT2  C    N    CP1    MULTIPLE= 2       2.75     2    180.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .30     4       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CT2  C    N    CP3    MULTIPLE= 2       2.75     2    180.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .30     4       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CT2  C    NH1  CT1    MULTIPLE= 2       1.60     1       .0000 !  ALLOW PEP 
      ! Revised to adjust NMA cis/trans energy difference. (LK) 
                                                   2.50     2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  CT2  C    NH1  CT2    MULTIPLE= 2       1.60     1       .0000 !  ALLOW PEP 
      ! Revised to adjust NMA cis/trans energy difference. (LK) 
                                                   2.50     2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  CT2  CA   CA   CA          3.10 {SD=      .309}    2    180.0000 ! ALLOW   ARO 
      ! JES 8/25/89 toluene and ethylbenzene 
 DIHEDRAL  CT2  CPH1 NR1  CPH2        3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FROM HA CPH1 NR1 CPH2 
 DIHEDRAL  CT2  CPH1 NR2  CPH2        3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FROM HA CPH1 NR2 CPH2 
 DIHEDRAL  CT2  CPH1 NR3  CPH2        2.50 {SD=      .344}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  CT2  CT1  NH1  C           1.80 {SD=      .406}    1       .0000 ! ALLOW PEP 
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c 
 DIHEDRAL  CT2  CT2  CPH1 CPH1         .40 {SD=      .860}    1       .0000 ! ALLOW ARO 
      ! 4-methylimidazole 4-21G//6-31G* rot bar. ADM JR., 9/4/89 
 DIHEDRAL  CT2  CT2  CT2  CT2          .15 {SD=     1.405}    1       .0000 ! ALLOW ALI 
      ! alkane update, adm jr., 3/2/92, butane trans/gauche 
 DIHEDRAL  CT2  CT2  NH1  C           1.80 {SD=      .406}    1       .0000 ! ALLOW PEP 
      ! from CT2  CT1  NH1  C, for lactams, adm jr. 
 DIHEDRAL  CT2  CY   CPT  CA          3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! JWK 
 DIHEDRAL  CT2  CY   CPT  CPT         3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      !JWK 
 DIHEDRAL  CT2  NH1  C    CP1    MULTIPLE= 2       1.60     1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                   2.50     2    180.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CT2  NH1  C    CT1    MULTIPLE= 2       1.60     1       .0000 !  ALLOW PEP 
      ! Revised to adjust NMA cis/trans energy difference. (LK) 
                                                   2.50     2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  CT2  SM   SM   CT2    MULTIPLE= 3       1.00     1       .0000 ! ALLOW   ALI SUL ION 
      ! improved CSSC dihedral in DMDS  5/15/92 (FL) 
                                                   4.10     2       .0000 ! ALLOW   ALI SUL ION 
      ! mp 6-311G** dimethyldisulfide,  3/26/92 (FL) 
                                                    .90     3       .0000 ! ALLOW   ALI SUL ION 
      ! improved CSSC dihedral in DMDS  5/15/92 (FL) 
 DIHEDRAL  CT3  C    N    CP1    MULTIPLE= 2       2.75     2    180.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .30     4       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CT3  C    N    CP3    MULTIPLE= 2       2.75     2    180.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .30     4       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CT3  C    NH1  CT1    MULTIPLE= 2       1.60     1       .0000 !  ALLOW PEP 
      ! Revised to adjust NMA cis/trans energy difference. (LK) 
                                                   2.50     2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  CT3  C    NH1  CT3    MULTIPLE= 2       1.60     1       .0000 !  ALLOW PEP 
      ! Revised to adjust NMA cis/trans energy difference. (LK) 
                                                   2.50     2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  CT3  CA   CA   CA          3.10 {SD=      .309}    2    180.0000 ! ALLOW   ARO 
      ! toluene, adm jr., 3/7/92 
 DIHEDRAL  CT3  CPH1 NR1  CPH2        3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FROM HA CPH1 NR1 CPH2 
 DIHEDRAL  CT3  CT1  NH1  C           1.80 {SD=      .406}    1       .0000 ! ALLOW PEP 
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c 
 DIHEDRAL  CT3  CT2  CA   CA           .23 {SD=     1.135}    2    180.0000 ! ALLOW   ARO ALI 
      ! ethylbenzene ethyl rotation, adm jr. 3/7/92 
 DIHEDRAL  CT3  CT2  CPH1 CPH1   MULTIPLE= 3        .20     1       .0000 ! ALLOW ARO 
      ! 4-ethylimidazole 4-21G rot bar, adm jr. 3/4/92 
                                                    .27     2       .0000 ! ALLOW ARO 
      ! 4-ethylimidazole 4-21G rot bar, adm jr. 3/4/92 
                                                    .00     3       .0000 ! ALLOW ARO 
      ! 4-ethylimidazole 4-21G rot bar, adm jr. 3/4/92 
 DIHEDRAL  CT3  CT2  CT2  CT2          .15 {SD=     1.405}    1       .0000 ! ALLOW ALI 
      ! alkane update, adm jr., 3/2/92, butane trans/gauche 
 DIHEDRAL  CT3  CT2  CT2  CT3          .15 {SD=     1.405}    1       .0000 ! ALLOW ALI 
      ! alkane update, adm jr., 3/2/92, butane trans/gauche 
 DIHEDRAL  CT3  CT2  CY   CA           .23 {SD=     1.135}    2    180.0000 ! ALLOW ARO 
      ! from ethylbenzene, adm jr., 3/7/92 
 DIHEDRAL  CT3  CT2  CY   CPT          .23 {SD=     1.135}    2    180.0000 ! ALLOW ARO 
      ! from ethylbenzene, adm jr., 3/7/92 
 DIHEDRAL  CT3  CT2  S    CT3    MULTIPLE= 2        .24     1    180.0000 ! ALOW    ALI SUL ION 
      ! expt. MeEtS,      3/26/92 (FL) 
                                                    .37     3       .0000 ! ALOW    ALI SUL ION 
      ! DTN 8/24/90 
 DIHEDRAL  CT3  NH1  C    CP1    MULTIPLE= 2       1.60     1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                   2.50     2    180.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  CT3  NH1  C    CT1    MULTIPLE= 2       1.60     1       .0000 !  ALLOW PEP 
      ! Revised to adjust NMA cis/trans energy difference. (LK) 
                                                   2.50     2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  CT3  S    CT2  CT2    MULTIPLE= 2        .24     1    180.0000 ! ALOW    ALI SUL ION 
      ! expt. MeEtS,      3/26/92 (FL) 
                                                    .37     3       .0000 ! ALOW    ALI SUL ION 
      ! expt. MeEtS,      3/26/92 (FL) 
 DIHEDRAL  CT3  SM   SM   CT3    MULTIPLE= 3       1.00     1       .0000 ! ALLOW   ALI SUL ION 
      ! improved CSSC dihedral in DMDS  5/15/92 (FL) 
                                                   4.10     2       .0000 ! ALLOW   ALI SUL ION 
      ! mp 6-311G** dimethyldisulfide,   3/26/92 (FL) 
                                                    .90     3       .0000 ! ALLOW   ALI SUL ION 
      ! improved CSSC dihedral in DMDS  5/15/92 (FL) 
 DIHEDRAL  CY   CA   NY   CPT         5.00 {SD=      .243}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  CY   CPT  CA   CA          3.00 {SD=      .314}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  CY   CPT  CPT  CA         10.00 {SD=      .172}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  H    NH1  C    CP1         2.50 {SD=      .344}    2    180.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  H    NH1  C    CT1         2.50 {SD=      .344}    2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  H    NH1  C    CT2         2.50 {SD=      .344}    2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  H    NH1  C    CT3         2.50 {SD=      .344}    2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  H    NH1  CT1  C            .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  H    NH1  CT1  CC           .00 {SD=999999.000}    1       .0000 ! ALLOW PEP POL 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 DIHEDRAL  H    NH1  CT1  CD           .00 {SD=999999.000}    1       .0000 ! ALLOW PEP POL 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 DIHEDRAL  H    NH1  CT1  CT1          .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  H    NH1  CT1  CT2          .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  H    NH1  CT1  CT3          .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  H    NH1  CT2  C            .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  H    NH1  CT2  CC           .00 {SD=999999.000}    1       .0000 ! ALLOW PEP POL 
      ! Alanine dipeptide; NMA; acetate; etc. backbone param. RLD 3/22/92 
 DIHEDRAL  H    NH1  CT2  CD           .00 {SD=999999.000}    1       .0000 ! ALLOW PEP POL 
      ! adm jr. 5/02/91, acetic acid pure solvent 
 DIHEDRAL  H    NH1  CT2  CT2          .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! from H    NH1  CT2  CT3, for lactams, adm jr. 
 DIHEDRAL  H    NH1  CT2  CT3          .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  H    NH2  CC   CT1         1.40 {SD=      .460}    2    180.0000 !  ALLOW   PEP POL ARO PRO 
      ! adm jr. 4/10/91, acetamide update 
 DIHEDRAL  H    NH2  CC   CT2         1.40 {SD=      .460}    2    180.0000 !  ALLOW   PEP POL ARO PRO 
      ! adm jr. 4/10/91, acetamide update 
 DIHEDRAL  H    NH2  CC   CT3         1.40 {SD=      .460}    2    180.0000 !  ALLOW   PEP POL ARO PRO 
      ! adm jr. 4/10/91, acetamide update 
 DIHEDRAL  H    NH2  CC   CP1         2.50 {SD=      .344}    2    180.0000 ! ALLOW PEP POL ARO PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  H    NR1  CPH1 CPH1        1.00 {SD=      .544}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 7/20/89 
 DIHEDRAL  H    NR1  CPH1 CT2         1.00 {SD=      .544}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 7/22/89, FROM HA CPH1 NR1 H 
 DIHEDRAL  H    NR1  CPH1 CT3         1.00 {SD=      .544}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 7/22/89, FROM HA CPH1 NR1 H 
 DIHEDRAL  H    NR3  CPH1 CPH1        1.40 {SD=      .460}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  H    NR3  CPH1 CT2         3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 7/22/89, FROM HC NR3 CPH1 HA 
 DIHEDRAL  H    NR3  CPH1 CT3         3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 7/22/89, FROM HC NR3 CPH1 HA 
 DIHEDRAL  H    NY   CA   CY           .80 {SD=      .608}    2    180.0000 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 DIHEDRAL  H    NY   CPT  CA           .80 {SD=      .608}    2    180.0000 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 DIHEDRAL  H    NY   CPT  CPT          .80 {SD=      .608}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  H    OH1  CA   CA           .99 {SD=      .547}    2    180.0000 ! ALLOW   ARO ALC 
      ! phenol OH rot bar, 3.37 kcal/mole, adm jr. 3/7/92 
 DIHEDRAL  H    OH1  CT1  CT1    MULTIPLE= 3       1.33     1       .0000 ! ALLOW ALC 
      ! 2-propanol OH hf/6-31g* torsional surface, adm jr., 3/2/93 
                                                    .18     2       .0000 ! ALLOW ALC 
      ! 2-propanol OH hf/6-31g* torsional surface, adm jr., 3/2/93 
                                                    .32     3       .0000 ! ALLOW ALC 
      ! 2-propanol OH hf/6-31g* torsional surface, adm jr., 3/2/93 
 DIHEDRAL  H    OH1  CT1  CT3    MULTIPLE= 3       1.33     1       .0000 ! ALLOW ALC 
      ! 2-propanol OH hf/6-31g* torsional surface, adm jr., 3/2/93 
                                                    .18     2       .0000 ! ALLOW ALC 
      ! 2-propanol OH hf/6-31g* torsional surface, adm jr., 3/2/93 
                                                    .32     3       .0000 ! ALLOW ALC 
      ! 2-propanol OH hf/6-31g* torsional surface, adm jr., 3/2/93 
 DIHEDRAL  H    OH1  CT2  CT1    MULTIPLE= 3       1.30     1       .0000 ! ALLOW ALC 
      ! ethanol OH hf/6-31g* torsional surface, adm jr., 3/2/93 
                                                    .30     2       .0000 ! ALLOW ALC 
      ! ethanol OH hf/6-31g* torsional surface, adm jr., 3/2/93 
                                                    .42     3       .0000 ! ALLOW ALC 
      ! ethanol OH hf/6-31g* torsional surface, adm jr., 3/2/93 
 DIHEDRAL  H    OH1  CT2  CT3    MULTIPLE= 3       1.30     1       .0000 ! ALLOW ALC 
      ! ethanol OH hf/6-31g* torsional surface, adm jr., 3/2/93 
                                                    .30     2       .0000 ! ALLOW ALC 
      ! ethanol OH hf/6-31g* torsional surface, adm jr., 3/2/93 
                                                    .42     3       .0000 ! ALLOW ALC 
      ! ethanol OH hf/6-31g* torsional surface, adm jr., 3/2/93 
 DIHEDRAL  HA   CA   CA   CA          3.50 {SD=      .291}    2    180.0000 ! ALLOW ARO 
      ! adm jr., 10/02/89 
 DIHEDRAL  HA   CA   CA   CPT         3.50 {SD=      .291}    2    180.0000 ! ALLOW   ARO 
      ! JWK 09/05/89 
 DIHEDRAL  HA   CA   CA   HA          2.50 {SD=      .344}    2    180.0000 ! ALLOW ARO 
      ! ADM JR., 10/02/89 
 DIHEDRAL  HA   CA   CPT  CPT         3.00 {SD=      .314}    2    180.0000 ! ALLOW   ARO 
      ! TRP (JES) 
 DIHEDRAL  HA   CA   CPT  CY          4.00 {SD=      .272}    2    180.0000 ! ALLOW   ARO 
      ! JWK 09/05/89 
 DIHEDRAL  HA   CA   CY   CPT         1.20 {SD=      .497}    2    180.0000 ! ALLOW ARO 
      ! JWK 
 DIHEDRAL  HA   CA   CY   CT2         1.20 {SD=      .497}    2    180.0000 ! ALLOW ARO 
      ! JWK 
 DIHEDRAL  HA   CA   NY   CPT         3.00 {SD=      .314}    2    180.0000 ! ALLOW   ARO 
      ! JWK 09/05/89 
 DIHEDRAL  HA   CA   NY   H           1.00 {SD=      .544}    2    180.0000 ! ALLOW   ARO 
      ! JWK 09/05/89 
 DIHEDRAL  HA   CC   NH2  H           1.40 {SD=      .460}    2    180.0000 !  ALLOW PEP POL 
      ! adm jr. 4/10/91, acetamide update 
 DIHEDRAL  HA   CP3  N    C            .00 {SD=999999.000}    3    180.0000 ! ALLOW PRO PEP 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HA   CP3  N    CP1          .10 {SD=     1.721}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HA   CP3  NP   CP1          .08 {SD=     1.924}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HA   CT1  CT2  CA           .04 {SD=     2.721}    3       .0000 ! ALLOW   ARO 
      ! 2.7 kcal/mole CH3 rot in ethylbenzene, adm jr, 3/7/92 
 DIHEDRAL  HA   CT2  CPH1 CPH1         .00 {SD=999999.000}    3       .0000 ! ALLOW ARO 
      ! 4-methylimidazole 4-21G//6-31G* rot bar. adm jr., 9/4/89 
 DIHEDRAL  HA   CT2  CY   CA           .25 {SD=     1.088}    2    180.0000 ! ALLOW ARO 
      ! JWK 
 DIHEDRAL  HA   CT2  CY   CPT          .25 {SD=     1.088}    2    180.0000 ! ALLOW ARO 
      ! JWK 
 DIHEDRAL  HA   CT2  NH1  C            .00 {SD=999999.000}    3       .0000 ! ALLOW PEP 
      ! LK for autogenerate dihe, sp2-methyl, no dihedral potential 
 DIHEDRAL  HA   CT2  NH1  H            .00 {SD=999999.000}    3       .0000 ! ALLOW PEP 
      ! LK for autogenerate dihe, sp2-methyl, no dihedral potential 
 DIHEDRAL  HA   CT2  S    CT3          .28 {SD=     1.028}    3       .0000 ! ALLOW   ALI SUL ION 
      ! DTN 8/24/90 
 DIHEDRAL  HA   CT3  CPH1 CPH1         .00 {SD=999999.000}    3       .0000 ! ALLOW ARO 
      ! 4-methylimidazole 4-21G//6-31G* rot bar. adm jr., 9/4/89 
 DIHEDRAL  HA   CT3  CS   HA           .16 {SD=     1.360}    3       .0000 ! ALLOW SUL 
      ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92 
 DIHEDRAL  HA   CT3  CT2  CA           .04 {SD=     2.721}    3       .0000 ! ALLOW   ARO 
      ! 2.7 kcal/mole CH3 rot in ethylbenzene, adm jr, 3/7/92 
 DIHEDRAL  HA   CT3  NH1  C            .00 {SD=999999.000}    3       .0000 ! ALLOW PEP 
      ! LK for autogenerate dihe, sp2-methyl, no dihedral potential 
 DIHEDRAL  HA   CT3  NH1  H            .00 {SD=999999.000}    3       .0000 ! ALLOW PEP 
      ! LK for autogenerate dihe, sp2-methyl, no dihedral potential 
 DIHEDRAL  HA   CT3  S    CT2          .28 {SD=     1.028}    3       .0000 ! ALLOW   ALI SUL ION 
      ! DTN 8/24/90 
 DIHEDRAL  HA   CY   CA   CPT         1.20 {SD=      .497}    2    180.0000 ! ALLOW   ARO 
      ! JWK 09/05/89 
 DIHEDRAL  HA   CY   CA   HA          1.20 {SD=      .497}    2    180.0000 ! ALLOW   ARO 
      ! JWK 09/05/89 
 DIHEDRAL  HA   CY   CPT  CA          3.00 {SD=      .314}    2    180.0000 ! ALLOW   ARO 
      ! JWK 09/05/89 
 DIHEDRAL  HA   CY   CPT  CPT         3.00 {SD=      .314}    2    180.0000 ! ALLOW   ARO 
      ! JWK 09/05/89 
 DIHEDRAL  HB   CP1  N    C            .80 {SD=      .608}    3       .0000 ! ALLOW PRO PEP 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HB   CP1  N    CP3          .10 {SD=     1.721}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HB   CP1  NP   CP3          .08 {SD=     1.924}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HB   CT1  NH1  C            .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  HB   CT1  NH1  H            .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  HB   CT2  NH1  C            .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  HB   CT2  NH1  H            .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  HB   CT3  NH1  C            .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  HB   CT3  NH1  H            .00 {SD=999999.000}    1       .0000 ! ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  HC   NP   CP1  C            .08 {SD=     1.924}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HC   NP   CP1  CC           .08 {SD=     1.924}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HC   NP   CP1  CD           .08 {SD=     1.924}    3       .0000 ! ALLOW PRO PEP 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HC   NP   CP1  CP2          .08 {SD=     1.924}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HC   NP   CP1  HB           .08 {SD=     1.924}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HC   NP   CP3  CP2          .08 {SD=     1.924}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HC   NP   CP3  HA           .08 {SD=     1.924}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  HP   CA   CA   CA          4.20 {SD=      .266}    2    180.0000 ! ALLOW   ARO 
      ! JES 8/25/89 benzene 
 DIHEDRAL  HP   CA   CA   CPT         3.00 {SD=      .314}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  HP   CA   CA   CT2         4.20 {SD=      .266}    2    180.0000 ! ALLOW   ARO 
      ! JES 8/25/89 toluene and ethylbenzene 
 DIHEDRAL  HP   CA   CA   CT3         4.20 {SD=      .266}    2    180.0000 ! ALLOW   ARO 
      ! toluene, adm jr., 3/7/92 
 DIHEDRAL  HP   CA   CA   HP          2.40 {SD=      .351}    2    180.0000 ! ALLOW   ARO 
      ! JES 8/25/89 benzene 
 DIHEDRAL  HP   CA   CPT  CPT         3.00 {SD=      .314}    2    180.0000 ! ALLOW   ARO 
      ! JWK indole 05/14/91 
 DIHEDRAL  HP   CA   CPT  CY          3.00 {SD=      .314}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  HP   CA   CY   CPT         2.00 {SD=      .385}    2    180.0000 ! ALLOW ARO 
      !adm jr., 12/30/91, for jwk 
 DIHEDRAL  HP   CA   CY   CT2         1.20 {SD=      .497}    2    180.0000 ! ALLOW ARO 
      ! JWK indole 05/14/91 
 DIHEDRAL  HP   CA   NY   CPT         2.00 {SD=      .385}    2    180.0000 ! ALLOW   ARO 
      !adm jr., 12/30/91, for jwk 
 DIHEDRAL  HP   CA   NY   H            .40 {SD=      .860}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  HP   CY   CA   HP          1.00 {SD=      .544}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  HP   CY   CPT  CA          2.80 {SD=      .325}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  HP   CY   CPT  CPT         2.80 {SD=      .325}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  HR1  CPH1 CPH1 CT2         1.00 {SD=      .544}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  HR1  CPH1 CPH1 CT3         1.00 {SD=      .544}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  HR1  CPH1 CPH1 HR1         1.00 {SD=      .544}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90, his 
 DIHEDRAL  HR1  CPH1 NR3  CPH2        2.50 {SD=      .344}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  HR1  CPH1 NR3  H           3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  HR1  CPH2 NR1  CPH1        3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  HR1  CPH2 NR1  H           1.00 {SD=      .544}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  HR1  CPH2 NR2  CPH1        3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  HR2  CPH2 NR3  CPH1        3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  HR2  CPH2 NR3  H            .00 {SD=999999.000}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90, YES, 0.0 
 DIHEDRAL  HR3  CPH1 CPH1 CT2         2.00 {SD=      .385}    2    180.0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 DIHEDRAL  HR3  CPH1 CPH1 CT3         2.00 {SD=      .385}    2    180.0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 DIHEDRAL  HR3  CPH1 CPH1 HR3         2.00 {SD=      .385}    2    180.0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 DIHEDRAL  HR3  CPH1 NR1  CPH2        3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 DIHEDRAL  HR3  CPH1 NR1  H           1.00 {SD=      .544}    2    180.0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 DIHEDRAL  HR3  CPH1 NR2  CPH2        3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 DIHEDRAL  HS   S    CT2  CT1    MULTIPLE= 3        .24     1       .0000 ! ALLOW   ALI SUL ION 
      ! methanethiol pure solvent, adm jr., 6/22/92 
                                                    .15     2       .0000 ! ALLOW   ALI SUL ION 
      ! methanethiol pure solvent, adm jr., 6/22/92 
                                                    .27     3       .0000 ! ALLOW   ALI SUL ION 
      ! methanethiol pure solvent, adm jr., 6/22/92 
 DIHEDRAL  HS   S    CT2  CT3    MULTIPLE= 3        .24     1       .0000 ! ALLOW   ALI SUL ION 
      ! ethanethiol C-C-S-H surface, adm jr., 4/18/93 
                                                    .15     2       .0000 ! ALLOW   ALI SUL ION 
      ! ethanethiol C-C-S-H surface, adm jr., 4/18/93 
                                                    .27     3       .0000 ! ALLOW   ALI SUL ION 
      ! ethanethiol C-C-S-H surface, adm jr., 4/18/93 
 DIHEDRAL  HS   S    CT2  HA           .20 {SD=     1.217}    3       .0000 ! ALLOW   ALI SUL ION 
      ! methanethiol pure solvent, adm jr., 6/22/92 
 DIHEDRAL  HS   S    CT3  HA           .20 {SD=     1.217}    3       .0000 ! ALLOW   ALI SUL ION 
      ! methanethiol pure solvent, adm jr., 6/22/92 
 DIHEDRAL  N    C    CP1  CP2    MULTIPLE= 2        .40     1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .60     2       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  N    C    CP1  HB     MULTIPLE= 2        .40     1    180.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .60     2       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  N    C    CP1  N      MULTIPLE= 2        .30     1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                   -.30     4       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  N    C    CT1  CT1          .00 {SD=999999.000}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  N    C    CT1  CT2          .00 {SD=999999.000}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  N    C    CT1  CT3          .00 {SD=999999.000}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  N    C    CT1  HB           .00 {SD=999999.000}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  N    C    CT2  HB           .00 {SD=999999.000}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  N    C    CT3  HA           .00 {SD=999999.000}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  N    CT1  CT2  CA           .04 {SD=     2.721}    3       .0000 ! ALLOW   ARO 
      ! 2.7 kcal/mole CH3 rot in ethylbenzene, adm jr, 3/7/92 
 DIHEDRAL  NH1  C    CP1  CP2    MULTIPLE= 2        .40     1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .60     2       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NH1  C    CP1  HB     MULTIPLE= 2        .40     1    180.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .60     2       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NH1  C    CP1  N      MULTIPLE= 2        .30     1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                   -.30     4       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NH1  C    CT1  CT1      MULTIPLE= 2      .00     1       .0000 !   ALLOW PEP 
                                                    .00     2       .0000 !   ALLOW PEP
      ! ala dipeptide corrxn for new C VDW Rmin, 4/10/93 (LK) 
 DIHEDRAL  NH1  C    CT1  CT2      MULTIPLE= 2      .00     1       .0000 !   ALLOW PEP
                                                    .00     2       .0000 !   ALLOW PEP
      ! ala dipeptide corrxn for new C VDW Rmin, 4/10/93 (LK) 
 DIHEDRAL  NH1  C    CT1  CT3      MULTIPLE= 2      .00     1       .0000 !   ALLOW PEP
                                                    .00     2       .0000 !   ALLOW PEP
      ! ala dipeptide corrxn for new C VDW Rmin, 4/10/93 (LK) 
 DIHEDRAL  NH1  C    CT1  HB       MULTIPLE= 2      .00     1       .0000 !   ALLOW PEP
                                                    .00     2       .0000 !   ALLOW PEP
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  NH1  C    CT1  NH1      MULTIPLE= 2      .60     1       .0000 !   ALLOW PEP 
                                                    .00     2       .0000 !   ALLOW PEP
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93 
 DIHEDRAL  NH1  C    CT2  CT2      MULTIPLE= 2      .00     1       .0000 !   ALLOW PEP
                                                    .00     2       .0000 !   ALLOW PEP
      ! from NH1  C    CT1  CT2, for lactams, adm jr. 
 DIHEDRAL  NH1  C    CT2  HA       MULTIPLE= 2      .00     1       .0000 !   ALLOW PEP
                                                    .00     2       .0000 !   ALLOW PEP
      ! LK for autogenerate dihe, sp2-methyl, no dihedral potential 
 DIHEDRAL  NH1  C    CT2  HB       MULTIPLE= 2      .00     1       .0000 !   ALLOW PEP
                                                    .00     2       .0000 !   ALLOW PEP
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  NH1  C    CT2  NH1      MULTIPLE= 2      .60     1       .0000 !   ALLOW PEP 
                                                    .00     2       .0000 !   ALLOW PEP
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93 
 DIHEDRAL  NH1  C    CT3  HA           .00 {SD=999999.000}    3       .0000 ! ALLOW PEP 
      ! LK for autogenerate dihe, sp2-methyl, no dihedral potential 
 DIHEDRAL  NH1  CT1  C    N            .40 {SD=      .860}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NH1  CT2  C    N            .40 {SD=      .860}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NH2  CC   CP1  CP2    MULTIPLE= 2        .40     1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .60     2       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NH2  CC   CP1  HB     MULTIPLE= 2        .40     1    180.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .60     2       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NH2  CC   CP1  N      MULTIPLE= 2        .30     1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                   -.30     4       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NH2  CC   CT1  NH1     MULTIPLE= 2       .05     6    180.0000 ! ALLOW PEP PRO
                                                    .00     2       .0000 !   ALLOW PEP
      ! satisfy multiple dihedral requirement
 DIHEDRAL  NH2  CC   CT1  HB      MULTIPLE= 2       .05     6    180.0000 ! ALLOW PEP PRO
                                                    .00     2       .0000 !   ALLOW PEP
      ! satisfy multiple dihedral requirement
 DIHEDRAL  NH2  CC   CT1  CT1     MULTIPLE= 2       .05     6    180.0000 ! ALLOW PEP PRO
                                                    .00     2       .0000 !   ALLOW PEP
      ! satisfy multiple dihedral requirement
      ! 08-16-96 added the following in oder to describe CB-CA-C-NT
 DIHEDRAL  NH2  CC   CT1  CT2     MULTIPLE= 2       .05     6    180.0000 ! ALLOW PEP PRO
                                                    .00     2       .0000 !   ALLOW PEP
      ! satisfy multiple dihedral requirement
 DIHEDRAL  NH2  CC   CT1  CT3     MULTIPLE= 2       .05     6    180.0000 ! ALLOW PEP PRO
                                                    .00     2       .0000 !   ALLOW PEP
      ! satisfy multiple dihedral requirement
 DIHEDRAL  NH2  CC   CT2  HA           .00 {SD=999999.000}    3    180.0000 ! ALLOW POL 
      ! adm jr., 8/13/90 acetamide geometry and vibrations 
 DIHEDRAL  NH3  CT1  C    N            .40 {SD=      .860}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NH3  CT1  C    NH1          .60 {SD=      .702}    1       .0000 ! ALLOW PEP PRO 
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93 
 DIHEDRAL  NH3  CT1  CC   NH2          .40 {SD=      .860}    1       .0000 ! ALLOW PEP PRO 
      ! Alanine dipeptide; NMA; acetate; etc. backbone param. RLD 3/22/92 
 DIHEDRAL  NH3  CT2  C    N            .40 {SD=      .860}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NH3  CT2  C    NH1          .40 {SD=      .860}    1       .0000 ! ALLOW PEP PRO 
      ! adm jr. 3/24/92, for PRES GLYP 
 DIHEDRAL  NH3  CT2  CC   NH2          .40 {SD=      .860}    1       .0000 ! ALLOW PEP PRO 
      ! Alanine dipeptide; NMA; acetate; etc. backbone param. RLD 3/22/92 
 DIHEDRAL  NP   CP1  C    N            .30 {SD=      .993}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NP   CP1  C    NH1          .30 {SD=      .993}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NP   CP1  CC   NH2          .30 {SD=      .993}    1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  NR1  CPH1 CPH1 CT2         3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FROM NR1 CPH1 CPH1 HA 
 DIHEDRAL  NR1  CPH1 CPH1 CT3         3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FROM NR1 CPH1 CPH1 HA 
 DIHEDRAL  NR1  CPH1 CPH1 HR3         3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 DIHEDRAL  NR1  CPH1 CT2  CT1          .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89 
 DIHEDRAL  NR1  CPH1 CT2  CT2          .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89 
 DIHEDRAL  NR1  CPH1 CT2  CT3          .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89 
 DIHEDRAL  NR1  CPH1 CT2  HA           .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89 
 DIHEDRAL  NR1  CPH1 CT3  HA           .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89 
 DIHEDRAL  NR1  CPH2 NR2  CPH1       14.00 {SD=      .145}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 DIHEDRAL  NR2  CPH1 CPH1 CT2         3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FROM NR2 CPH1 CPH1 HA 
 DIHEDRAL  NR2  CPH1 CPH1 CT3         3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/22/89, FROM NR2 CPH1 CPH1 HA 
 DIHEDRAL  NR2  CPH1 CPH1 HR3         3.00 {SD=      .314}    2    180.0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 DIHEDRAL  NR2  CPH1 CPH1 NR1        14.00 {SD=      .145}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 DIHEDRAL  NR2  CPH1 CT2  CT1          .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! HIS CB-CG TORSION, 
 DIHEDRAL  NR2  CPH1 CT2  CT2          .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! HIS CB-CG TORSION, 
 DIHEDRAL  NR2  CPH1 CT2  CT3          .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! HIS CB-CG TORSION, 
 DIHEDRAL  NR2  CPH1 CT2  HA           .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89 
 DIHEDRAL  NR2  CPH1 CT3  HA           .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89 
 DIHEDRAL  NR2  CPH2 NR1  CPH1       14.00 {SD=      .145}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 DIHEDRAL  NR2  CPH2 NR1  H           1.00 {SD=      .544}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 DIHEDRAL  NR3  CPH1 CPH1 CT2         2.50 {SD=      .344}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  NR3  CPH1 CPH1 CT3         2.50 {SD=      .344}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  NR3  CPH1 CPH1 HR1         2.50 {SD=      .344}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  NR3  CPH1 CPH1 NR3        12.00 {SD=      .157}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  NR3  CPH1 CT2  CT1          .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89 
 DIHEDRAL  NR3  CPH1 CT2  CT2          .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89 
 DIHEDRAL  NR3  CPH1 CT2  CT3          .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89 
 DIHEDRAL  NR3  CPH1 CT2  HA           .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89 
 DIHEDRAL  NR3  CPH1 CT3  HA           .19 {SD=     1.248}    3       .0000 ! ALLOW ARO 
      ! 4-METHYLIMIDAZOLE 4-21G//6-31G* ROT BAR. ADM JR., 9/4/89 
 DIHEDRAL  NR3  CPH2 NR3  CPH1       12.00 {SD=      .157}    2    180.0000 ! ALLOW ARO 
      ! his, ADM JR., 7/20/89 
 DIHEDRAL  NR3  CPH2 NR3  H           1.40 {SD=      .460}    2    180.0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 DIHEDRAL  NY   CA   CY   CPT         4.00 {SD=      .272}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  NY   CA   CY   CT2         3.50 {SD=      .291}    2    180.0000 ! ALLOW ARO 
      ! JWK 
 DIHEDRAL  NY   CA   CY   HA          3.50 {SD=      .291}    2    180.0000 ! ALLOW   ARO 
      ! JWK 09/05/89 
 DIHEDRAL  NY   CA   CY   HP          3.50 {SD=      .291}    2    180.0000 ! ALLOW   ARO 
      ! JWK indole 05/14/91 
 DIHEDRAL  NY   CPT  CA   CA          2.80 {SD=      .325}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  NY   CPT  CA   HA          4.00 {SD=      .272}    2    180.0000 ! ALLOW   ARO 
      ! JWK 09/05/89 
 DIHEDRAL  NY   CPT  CA   HP          3.00 {SD=      .314}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  NY   CPT  CPT  CA         10.00 {SD=      .172}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  NY   CPT  CPT  CY          5.00 {SD=      .243}    2    180.0000 ! ALLOW   ARO 
      ! JWK 05/14/91 fit to indole 
 DIHEDRAL  O    C    CP1  CP2    MULTIPLE= 2        .40     1    180.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .60     2       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  O    C    CP1  HB     MULTIPLE= 2        .40     1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .60     2       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  O    C    CP1  N           -.30 {SD=999999.000}    4       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  O    C    CT1  CT1         1.40 {SD=      .460}    1       .0000 !   ALLOW PEP 
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c 
 DIHEDRAL  O    C    CT1  CT2         1.40 {SD=      .460}    1       .0000 !   ALLOW PEP 
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c 
 DIHEDRAL  O    C    CT1  CT3         1.40 {SD=      .460}    1       .0000 !   ALLOW PEP 
      ! ala dipeptide update for new C VDW Rmin, adm jr., 3/3/93c 
 DIHEDRAL  O    C    CT1  HB           .00 {SD=999999.000}    1       .0000 !   ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  O    C    CT1  NH1          .00 {SD=999999.000}    1       .0000 !   ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  O    C    CT1  NH3          .00 {SD=999999.000}    1       .0000 ! ALLOW PEP PRO 
      ! Backbone parameter set made complete RLD 8/8/90 
 DIHEDRAL  O    C    CT2  CT2         1.40 {SD=      .460}    1       .0000 !   ALLOW PEP 
      ! from O    C    CT1  CT2, for lactams, adm jr. 
 DIHEDRAL  O    C    CT2  HA           .00 {SD=999999.000}    3    180.0000 ! ALLOW POL 
      ! adm jr., 8/13/90 acetamide geometry and vibrations 
 DIHEDRAL  O    C    CT2  HB           .00 {SD=999999.000}    1       .0000 !   ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  O    C    CT2  NH1          .00 {SD=999999.000}    1       .0000 !   ALLOW PEP 
      ! Alanine Dipeptide ab initio calc's (LK) 
 DIHEDRAL  O    C    CT2  NH3          .00 {SD=999999.000}    1       .0000 ! ALLOW PEP PRO 
      ! Backbone parameter set made complete RLD 8/8/90 
 DIHEDRAL  O    C    CT3  HA           .00 {SD=999999.000}    3    180.0000 ! ALLOW POL 
      ! adm jr., 8/13/90 acetamide geometry and vibrations 
 DIHEDRAL  O    C    N    CP1    MULTIPLE= 2       2.75     2    180.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .30     4       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  O    C    N    CP3    MULTIPLE= 2       2.75     2    180.0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .30     4       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  O    C    NH1  CT1         2.50 {SD=      .344}    2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  O    C    NH1  CT2         2.50 {SD=      .344}    2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  O    C    NH1  CT3         2.50 {SD=      .344}    2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  O    C    NH1  H           2.50 {SD=      .344}    2    180.0000 !  ALLOW PEP 
      ! Gives appropriate NMA cis/trans barrier. (LK) 
 DIHEDRAL  O    CC   CP1  CP2    MULTIPLE= 2        .40     1    180.0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .60     2       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  O    CC   CP1  HB     MULTIPLE= 2        .40     1       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
                                                    .60     2       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  O    CC   CP1  N           -.30 {SD=999999.000}    4       .0000 ! ALLOW PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  O    CC   CT2  HA           .00 {SD=999999.000}    3    180.0000 ! ALLOW POL 
      ! adm jr. 4/05/91, for asn,asp,gln,glu and cters 
 DIHEDRAL  O    CC   NH2  H           1.40 {SD=      .460}    2    180.0000 !  ALLOW PEP POL ARO PRO 
      ! adm jr. 4/10/91, acetamide update 
 DIHEDRAL  OB   CD   OS   CT2    MULTIPLE= 2        .96     1    180.0000 ! ALLOW   PEP POL 
      ! adm jr. 3/19/92, from lipid methyl acetate 
                                                   3.85     2    180.0000 ! ALLOW   PEP POL 
      ! adm jr. 3/19/92, from lipid methyl acetate 
 DIHEDRAL  OB   CD   OS   CT3    MULTIPLE= 2        .96     1    180.0000 ! ALLOW   PEP POL 
      ! adm jr. 3/19/92, from lipid methyl acetate 
                                                   3.85     2    180.0000 ! ALLOW   PEP POL 
      ! adm jr. 3/19/92, from lipid methyl acetate 
 DIHEDRAL  OC   CA   CA   CA          3.10 {SD=      .309}    2    180.0000 ! ALLOW   ARO 
      ! adm jr. 8/27/91, phenoxide 
 DIHEDRAL  OC   CA   CA   HP          4.20 {SD=      .266}    2    180.0000 ! ALLOW   ARO 
      ! adm jr. 8/27/91, phenoxide 
 DIHEDRAL  OC   CC   CP1  CP2          .16 {SD=     1.360}    3       .0000 ! ALLOW PEP PRO POL 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  OC   CC   CP1  HB           .16 {SD=     1.360}    3       .0000 ! ALLOW PEP PRO POL 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  OC   CC   CP1  N            .16 {SD=     1.360}    3       .0000 ! ALLOW PEP PRO POL 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  OC   CC   CP1  NP           .16 {SD=     1.360}    3       .0000 ! ALLOW PEP PRO POL 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  OH1  CA   CA   CA          3.10 {SD=      .309}    2    180.0000 ! ALLOW   ARO 
      ! JES 8/25/89 phenol 
 DIHEDRAL  OH1  CA   CA   HP          4.20 {SD=      .266}    2    180.0000 ! ALLOW   ARO 
      ! JES 8/25/89 phenol 
 DIHEDRAL  S    CT2  CT2  HA           .01 {SD=     5.441}    3       .0000 ! ALLOW   ALI SUL ION 
      ! DTN 8/24/90 
 DIHEDRAL  SM   CT2  CT2  HA           .01 {SD=     5.441}    3       .0000 ! ALLOW   ALI SUL ION 
      ! DTN 8/24/90 
 DIHEDRAL  SM   SM   CT2  CT1          .31 {SD=      .977}    3       .0000 ! ALLOW  SUL ALI 
      ! S-S for cys-cys, dummy parameter for now ... DTN  9/04/90 
 DIHEDRAL  SM   SM   CT2  CT2          .31 {SD=      .977}    3       .0000 ! ALLOW  SUL ALI 
      ! S-S for cys-cys, dummy parameter for now ... DTN  9/04/90 
 DIHEDRAL  SM   SM   CT2  HA           .16 {SD=     1.369}    3       .0000 ! ALLOW   ALI SUL ION 
      ! expt. dimethyldisulfide,    3/26/92 (FL) 
 DIHEDRAL  SM   SM   CT3  HA           .16 {SD=     1.369}    3       .0000 ! ALLOW   ALI SUL ION 
      ! expt. dimethyldisulfide,    3/26/92 (FL) 
 DIHEDRAL  SS   CS   CT3  HA           .15 {SD=     1.405}    3       .0000 ! ALLOW SUL 
      ! ethylthiolate 6-31+G* geom/freq, adm jr., 6/1/92 
 DIHEDRAL  X    C    C    X           4.00 {SD=      .272}    2    180.0000 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 DIHEDRAL  X    C    NC2  X           2.25 {SD=      .363}    2    180.0000 ! ALLOW   PEP POL ARO 
      ! 9.0->2.25 GUANIDINIUM (KK) 
 DIHEDRAL  X    CD   OH1  X           2.05 {SD=      .380}    2    180.0000 ! ALLOW   PEP POL ARO ALC 
      ! adm jr, 10/17/90, acetic acid C-Oh rotation barrier 
 DIHEDRAL  X    CD   OS   X           2.05 {SD=      .380}    2    180.0000 ! ALLOW   PEP POL 
      ! adm jr. 3/19/92, from lipid methyl acetate 
 DIHEDRAL  X    CP1  C    X            .00 {SD=999999.000}    6    180.0000 ! ALLOW   POL PEP PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  X    CP1  CC   X            .00 {SD=999999.000}    6    180.0000 ! ALLOW   POL PEP 
      ! changed to 0.0 RLD 5/19/92 
 DIHEDRAL  X    CP1  CD   X            .00 {SD=999999.000}    6    180.0000 ! ALLOW   POL PEP 
      ! Alanine dipeptide; NMA; acetate; etc. backbone param. RLD 3/22/92 
 DIHEDRAL  X    CP1  CP2  X            .14 {SD=     1.454}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  X    CP2  CP2  X            .16 {SD=     1.360}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  X    CP3  CP2  X            .14 {SD=     1.454}    3       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 DIHEDRAL  X    CPA  CPB  X            .00 {SD=999999.000}    2       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): dummy for "auto dihe" (KK 05/13/91) 
 DIHEDRAL  X    CPA  CPM  X            .00 {SD=999999.000}    2       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): dummy for "auto dihe" (KK 05/13/91) 
 DIHEDRAL  X    CPB  C    X           3.00 {SD=      .314}    2    180.0000 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 DIHEDRAL  X    CPB  CPB  X            .00 {SD=999999.000}    2       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): dummy for "auto dihe" (KK 05/13/91) 
 DIHEDRAL  X    CPB  CT2  X            .00 {SD=999999.000}    6       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 DIHEDRAL  X    CPB  CT3  X            .00 {SD=999999.000}    6       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 DIHEDRAL  X    CPT  CPT  X            .00 {SD=999999.000}    2    180.0000 ! ALLOW   ARO 
      ! JWK indole 05/14/91 
 DIHEDRAL  X    CT1  CC   X            .05 {SD=     2.433}    6    180.0000 ! ALLOW   POL PEP 
      ! For side chains of asp,asn,glu,gln, (n=6) from KK(LK) 
 DIHEDRAL  X    CT1  CD   X            .00 {SD=999999.000}    6    180.0000 ! ALLOW   POL PEP 
      ! adm jr. 3/19/92, from lipid methyl acetate 
 DIHEDRAL  X    CT1  CT1  X            .20 {SD=     1.217}    3       .0000 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 DIHEDRAL  X    CT1  CT2  X            .20 {SD=     1.217}    3       .0000 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 DIHEDRAL  X    CT1  CT3  X            .20 {SD=     1.217}    3       .0000 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 DIHEDRAL  X    CT1  NH3  X            .10 {SD=     1.721}    3       .0000 ! ALLOW   ALI POL 
      ! 0.715->0.10 METHYLAMMONIUM (KK) 
 DIHEDRAL  X    CT1  OH1  X            .14 {SD=     1.454}    3       .0000 ! ALLOW   ALI ALC ARO 
      ! EMB  11/21/89 methanol vib fit 
 DIHEDRAL  X    CT1  OS   X           -.10 {SD=999999.000}    3       .0000 ! ALLOW   PEP POL 
      ! adm jr. 3/19/92, from lipid methyl acetate 
 DIHEDRAL  X    CT2  CA   X            .00 {SD=999999.000}    6       .0000 ! ALLOW   ALI ARO 
      ! toluene, adm jr., 3/7/92 
 DIHEDRAL  X    CT2  CC   X            .05 {SD=     2.433}    6    180.0000 ! ALLOW   POL PEP 
      ! For side chains of asp,asn,glu,gln, (n=6) from KK(LK) 
 DIHEDRAL  X    CT2  CD   X            .00 {SD=999999.000}    6    180.0000 ! ALLOW   POL PEP 
      ! adm jr. 3/19/92, from lipid methyl acetate 
 DIHEDRAL  X    CT2  CT2  X            .20 {SD=     1.232}    3       .0000 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 DIHEDRAL  X    CT2  CT3  X            .16 {SD=     1.360}    3       .0000 ! ALLOW   ALI 
      ! rotation barrier in Ethane (SF) 
 DIHEDRAL  X    CT2  NC2  X            .00 {SD=999999.000}    6    180.0000 ! ALLOW   ALI POL 
      ! methylguanidinium, adm jr., 3/26/92 
 DIHEDRAL  X    CT2  NH3  X            .10 {SD=     1.721}    3       .0000 ! ALLOW   ALI POL 
      ! 0.715->0.10 METHYLAMMONIUM (KK) 
 DIHEDRAL  X    CT2  OH1  X            .14 {SD=     1.454}    3       .0000 ! ALLOW   ALI ALC ARO 
      ! EMB  11/21/89 methanol vib fit 
 DIHEDRAL  X    CT2  OS   X           -.10 {SD=999999.000}    3       .0000 ! ALLOW   PEP POL 
      ! adm jr. 3/19/92, from lipid methyl acetate 
 DIHEDRAL  X    CT3  CA   X            .00 {SD=999999.000}    6       .0000 ! ALLOW   ALI ARO 
      ! toluene, adm jr., 3/7/92 
 DIHEDRAL  X    CT3  CC   X            .05 {SD=     2.433}    6    180.0000 ! ALLOW   POL PEP 
      ! For side chains of asp,asn,glu,gln, (n=6) from KK(LK) 
 DIHEDRAL  X    CT3  CD   X            .00 {SD=999999.000}    6    180.0000 ! ALLOW   POL PEP 
      ! adm jr. 3/19/92, from lipid methyl acetate 
 DIHEDRAL  X    CT3  CT3  X            .15 {SD=     1.382}    3       .0000 ! ALLOW   ALI 
      ! alkane update, adm jr., 3/2/92 
 DIHEDRAL  X    CT3  NC2  X            .00 {SD=999999.000}    6    180.0000 ! ALLOW   ALI POL 
      ! methylguanidinium, adm jr., 3/26/92 
 DIHEDRAL  X    CT3  NH2  X            .11 {SD=     1.641}    3       .0000 ! ALLOW   POL 
      ! methylamine geom/freq, adm jr., 6/2/92 
 DIHEDRAL  X    CT3  NH3  X            .09 {SD=     1.814}    3       .0000 ! ALLOW   ALI POL 
      ! fine-tuned to ab initio; METHYLAMMONIUM, KK 03/10/92 
 DIHEDRAL  X    CT3  OH1  X            .14 {SD=     1.454}    3       .0000 ! ALLOW   ALI ALC ARO 
      ! EMB  11/21/89 methanol vib fit 
 DIHEDRAL  X    CT3  OS   X           -.10 {SD=999999.000}    3       .0000 ! ALLOW   PEP POL 
      ! adm jr. 3/19/92, from lipid methyl acetate 
 DIHEDRAL  X    FE   CM   X            .05 {SD=     2.433}    4       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): ligands (KK 05/13/91) 
 DIHEDRAL  X    FE   NPH  X            .00 {SD=999999.000}    2       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): for "ic para" only (KK 05/13/91) 
 DIHEDRAL  X    FE   OM   X            .00 {SD=999999.000}    4       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): ligands (KK 05/13/91) 
 DIHEDRAL  X    NPH  CPA  X            .00 {SD=999999.000}    2       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): dummy for "auto dihe" (KK 05/13/91) 

 IMPROPER  CPB  CPA  NPH  CPA        20.80 {SD=      .119}    0       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 IMPROPER  CPB  X    X    C          90.00 {SD=      .057}    0       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 IMPROPER  CT2  X    X    CPB        90.00 {SD=      .057}    0       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 IMPROPER  CT3  X    X    CPB        90.00 {SD=      .057}    0       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 IMPROPER  HA   C    C    HA         20.00 {SD=      .122}    0       .0000 ! ALLOW   PEP POL ARO 
      ! Heme vinyl substituent (KK, from propene (JCS)) 
 IMPROPER  HA   CPA  CPA  CPM        29.40 {SD=      .100}    0       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 IMPROPER  HA   CPB  C    C          20.00 {SD=      .122}    0       .0000 ! ALLOW HEM ARO 
      ! Heme (6-liganded): substituents (KK 05/13/91) 
 IMPROPER  HA   HA   C    C          20.00 {SD=      .122}    0    180.0000 ! ALLOW   PEP POL ARO 
      ! Heme vinyl substituent (KK, from propene (JCS)) 
 IMPROPER  HR1  NR1  NR2  CPH2         .50 {SD=      .770}    0       .0000 ! ALLOW ARO 
      ! his, adm jr., 7/05/90 
 IMPROPER  HR1  NR2  NR1  CPH2         .50 {SD=      .770}    0       .0000 ! ALLOW ARO 
      ! his, adm jr., 7/05/90 
 IMPROPER  HR3  CPH1 NR1  CPH1         .50 {SD=      .770}    0       .0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 IMPROPER  HR3  CPH1 NR2  CPH1         .50 {SD=      .770}    0       .0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 IMPROPER  HR3  CPH1 NR3  CPH1        1.00 {SD=      .544}    0       .0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 IMPROPER  HR3  NR1  CPH1 CPH1         .50 {SD=      .770}    0       .0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 IMPROPER  HR3  NR2  CPH1 CPH1         .50 {SD=      .770}    0       .0000 ! ALLOW ARO 
      ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 IMPROPER  N    C    CP1  CP3          .00 {SD=999999.000}    0       .0000 ! ALLOW PRO 
      ! 6-31g* AcProNH2 and ProNH2  RLD 5/19/92 
 IMPROPER  NC2  X    X    C          40.00 {SD=      .086}    0       .0000 ! ALLOW   PEP POL ARO 
      ! 5.75->40.0 GUANIDINIUM (KK) 
 IMPROPER  NH1  X    X    H          20.00 {SD=      .122}    0       .0000 ! ALLOW   PEP POL ARO 
      ! NMA Vibrational Modes (LK) 
 IMPROPER  NH2  X    X    H           4.00 {SD=      .272}    0       .0000 ! ALLOW   POL 
      ! adm jr., 8/13/90 acetamide geometry and vibrations 
 IMPROPER  NPH  CPA  CPA  FE        137.40 {SD=      .046}    0       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 IMPROPER  NPH  CPA  CPB  CPB        40.60 {SD=      .085}    0       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 IMPROPER  NPH  CPA  CPM  CPA        18.30 {SD=      .127}    0       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 IMPROPER  NPH  CPM  CPB  CPA        32.70 {SD=      .095}    0       .0000 ! ALLOW HEM 
      ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 IMPROPER  NR1  CPH1 CPH2 H            .45 {SD=      .811}    0       .0000 ! ALLOW ARO 
      ! his, adm jr., 7/05/90 
 IMPROPER  NR1  CPH2 CPH1 H            .45 {SD=      .811}    0       .0000 ! ALLOW ARO 
      ! his, adm jr., 7/05/90 
 IMPROPER  NR3  CPH1 CPH2 H           1.20 {SD=      .497}    0       .0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 IMPROPER  NR3  CPH2 CPH1 H           1.20 {SD=      .497}    0       .0000 ! ALLOW ARO 
      ! his, adm jr., 6/27/90 
 IMPROPER  NY   CA   CY   CPT       100.00 {SD=      .054}    0       .0000 ! ALLOW ARO 
      !adm jr., 5/15/91,  indole 3-21G HE1 out-of-plane surf. 
 IMPROPER  O    CP1  NH2  CC         45.00 {SD=      .081}    0       .0000 ! ALLOW PEP POL PRO 
      ! 6-31g* AcProNH2 and ProNH2  RLD 5/19/92 
 IMPROPER  O    CT1  NH2  CC         45.00 {SD=      .081}    0       .0000 ! ALLOW PEP POL 
      ! adm jr., 8/13/90 acetamide geometry and vibrations 
 IMPROPER  O    CT2  NH2  CC         45.00 {SD=      .081}    0       .0000 ! ALLOW PEP POL 
      ! adm jr., 8/13/90 acetamide geometry and vibrations 
 IMPROPER  O    CT3  NH2  CC         45.00 {SD=      .081}    0       .0000 ! ALLOW PEP POL 
      ! adm jr., 8/13/90 acetamide geometry and vibrations 
 IMPROPER  O    HA   NH2  CC         45.00 {SD=      .081}    0       .0000 ! ALLOW PEP POL PRO 
      ! adm jr., 5/13/91, formamide geometry and vibrations 
 IMPROPER  O    N    CT2  CC        120.00 {SD=      .050}    0       .0000 ! ALLOW PEP POL PRO 
      ! 6-31g* AcProNH2 and ProNH2  RLD 5/19/92 
 IMPROPER  O    NH2  CP1  CC         45.00 {SD=      .081}    0       .0000 ! ALLOW PEP POL PRO 
      ! 6-31g* AcProNH2 and ProNH2  RLD 5/19/92 
 IMPROPER  O    NH2  CT1  CC         45.00 {SD=      .081}    0       .0000 ! ALLOW PEP POL 
      ! adm jr., 8/13/90 acetamide geometry and vibrations 
 IMPROPER  O    NH2  CT2  CC         45.00 {SD=      .081}    0       .0000 ! ALLOW PEP POL 
      ! adm jr., 8/13/90 acetamide geometry and vibrations 
 IMPROPER  O    NH2  CT3  CC         45.00 {SD=      .081}    0       .0000 ! ALLOW PEP POL 
      ! adm jr., 8/13/90 acetamide geometry and vibrations 
 IMPROPER  O    NH2  HA   CC         45.00 {SD=      .081}    0       .0000 ! ALLOW PEP POL 
      ! adm jr., 5/13/91, formamide geometry and vibrations 
 IMPROPER  O    X    X    C         120.00 {SD=      .050}    0       .0000 ! ALLOW   PEP POL ARO 
      ! NMA Vibrational Modes (LK) 
 IMPROPER  OB   X    X    CD        100.00 {SD=      .054}    0       .0000 ! ALLOW   ALC ARO POL 
      ! adm jr., 10/17/90, acetic acid vibrations 
 IMPROPER  OC   X    X    CC         96.00 {SD=      .056}    0       .0000 ! ALLOW   PEP POL ARO ION 
      ! 90.0->96.0 acetate, single impr (KK) 


 NBONDS 

! -----NONBONDED-LIST-OPTIONS-------------------------------
  CUTNB=  13.000 TOLERANCE=    .500 WMIN=   1.500 ATOM   
  INHIBIT=    .250                                       
! -----ELECTROSTATIC OPTIONS--------------------------------
  EPS=   1.000 E14FAC=   1.000 CDIELECTRIC SHIFT         
! -----VAN DER WAALS OPTIONS--------------------------------
  VSWITCH                                                
! -----SWITCHING /SHIFTING PARAMETERS-----------------------
  CTONNB=  10.000 CTOFNB=  12.000                        
! -----EXCLUSION LIST OPTIONS-------------------------------
  NBXMOD=   5                                            
! ----------------------------------------------------------

 END 

 !                    EPS     SIGMA    EPS(1:4)  SIGMA(1:4)

 NONBONDED  C        .1100   3.5636    .1100   3.5636 ! ALLOW   PEP POL ARO  ! NMA pure solvent, adm jr., 3/3/93 
 NONBONDED  CA       .0700   3.5501    .0700   3.5501 ! ALLOW   ARO  ! benzene (JES) 
 NONBONDED  CAL      .1200   3.0469    .1200   3.0469 ! ALLOW  ION   ! optimized before july 89, and need polarization term (BP) 
 NONBONDED  CC       .0700   3.5636    .0700   3.5636 ! ALLOW   PEP POL ARO  ! adm jr. 3/3/92, acetic acid heat of solvation 
 NONBONDED  CD       .0700   3.5636    .0700   3.5636 ! ALLOW  POL   ! adm jr. 3/19/92, acetate a.i. and dH of solvation 
 NONBONDED  CM       .1100   3.7418    .1100   3.7418 ! ALLOW HEM   ! Heme (6-liganded): CO ligand carbon (KK 05/13/91) 
 NONBONDED  CP1      .0200   4.0536    .0100   3.3854 ! ALLOW   ALI ! alkane update, adm jr., 3/2/92 
 NONBONDED  CP2      .0550   3.8754    .0100   3.3854 ! ALLOW   ALI ! alkane update, adm jr., 3/2/92 
 NONBONDED  CP3      .0550   3.8754    .0100   3.3854 ! ALLOW   ALI ! alkane update, adm jr., 3/2/92 
 NONBONDED  CPA      .0900   3.2072    .0900   3.2072 ! ALLOW HEM   ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 NONBONDED  CPB      .0900   3.2072    .0900   3.2072 ! ALLOW HEM   ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 NONBONDED  CPH1     .0500   3.2072    .0500   3.2072 ! ALLOW ARO   ! adm jr., 10/23/91, imidazole solvation and sublimation 
 NONBONDED  CPH2     .0500   3.2072    .0500   3.2072 ! ALLOW ARO   ! adm jr., 10/23/91, imidazole solvation and sublimation 
 NONBONDED  CPM      .0900   3.2072    .0900   3.2072 ! ALLOW HEM   ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 NONBONDED  CPT      .0900   3.2072    .0900   3.3854 ! ALLOW   ARO ! benzene (JES) 
 NONBONDED  CS       .1100   3.9200    .1100   3.9200 ! ALLOW SUL   ! methylthiolate to water and F.E. of solvation, adm jr. 6/1/92 
 NONBONDED  CT1      .0200   4.0536    .0100   3.3854 ! ALLOW   ALI ! isobutane pure solvent properties, adm jr, 2/3/92 
 NONBONDED  CT2      .0550   3.8754    .0100   3.3854 ! ALLOW   ALI ! propane pure solvent properties, adm jr, 2/3/92 
 NONBONDED  CT3      .0800   3.6705    .0100   3.3854 ! ALLOW   ALI ! methane/ethane a.i. and ethane pure solvent, adm jr, 2/3/92 
 NONBONDED  CY       .0700   3.5501    .0700   3.5501 ! ALLOW   ARO ! TRP, JWK 08/29/89 
 NONBONDED  FE       .0001   1.1582    .0001   1.1582 ! ALLOW HEM   ! Heme (6-liganded): Iron atom (KK 05/13/91) 
 NONBONDED  H        .0460    .4000    .0460    .4000 ! ALLOW PEP POL SUL ARO ALC  ! same as TIP3P hydrogen, adm jr., 7/20/89 
 NONBONDED  HA       .0220   2.3520    .0220   2.3520 ! ALLOW PEP ALI POL SUL ARO PRO ALC  ! methane/ethane a.i. and ethane pure solvent, adm jr, 2/3/92 
 NONBONDED  HB       .0220   2.3520    .0220   2.3520 ! ALLOW PEP ALI POL SUL ARO PRO ALC  ! methane/ethane a.i. and ethane pure solvent, adm jr, 2/3/92 
 NONBONDED  HC       .0460    .4000    .0460    .4000 ! ALLOW POL   ! new, small polar Hydrogen, see also adm jr. JG 8/27/89 
 NONBONDED  HP       .0300   2.4200    .0300   2.4200 ! ALLOW ARO   ! JES 8/25/89 values from Jorgensen fit to hydration energy 
 NONBONDED  HR1      .0460   1.6036    .0460   1.6036 ! ALLOW ARO   ! adm jr., 6/27/90, his 
 NONBONDED  HR2      .0460   1.2473    .0460   1.2473 ! ALLOW ARO   ! adm jr., 6/27/90, his 
 NONBONDED  HR3      .0078   2.6157    .0078   2.6157 ! ALLOW ARO   ! adm jr., 3/24/92, maintain old aliphatic H VDW params 
 NONBONDED  HS       .1000    .8018    .1000    .8018 ! ALLOW SUL   ! methanethiol pure solvent, adm jr., 6/22/92 
 NONBONDED  HT       .0460    .4000    .0460    .4000 ! ALLOW WAT   !TIP3P HYDROGEN PARAMETERS, adm jr., NBFIX obsolete 
 NONBONDED  N        .2000   3.2963    .0001   3.2963 ! ALLOW   PRO ! 6-31g* AcProNH2, ProNH2, 6-31g*//3-21g AcProNHCH3 RLD 4/23/93 
 NONBONDED  NC2      .2000   3.2963    .2000   3.2963 ! ALLOW   POL ! JG 8/27/89; note: NH1 in ARG was changed to NC2. 
 NONBONDED  NH1      .2000   3.2963    .2000   2.7618 ! ALLOW   PEP POL ARO  ! This 1,4 vdW allows the C5 dipeptide minimum to exist.(LK) 
 NONBONDED  NH2      .2000   3.2963    .2000   3.2963 ! ALLOW   POL ! adm jr. 
 NONBONDED  NH3      .2000   3.2963    .2000   3.2963 ! ALLOW   POL ! adm jr. 
 NONBONDED  NP       .2000   3.2963    .2000   3.2963 ! ALLOW  PRO  ! N-terminal proline; from 6-31g* +ProNH2  RLD 9/28/90 
 NONBONDED  NPH      .2000   3.2963    .2000   3.2963 ! ALLOW HEM   ! Heme (6-liganded): porphyrin macrocycle (KK 05/13/91) 
 NONBONDED  NR1      .2000   3.2963    .2000   3.2963 ! ALLOW ARO   ! His, adm jr., 9/4/89 
 NONBONDED  NR2      .2000   3.2963    .2000   3.2963 ! ALLOW ARO   ! His, adm jr., 9/4/89 
 NONBONDED  NR3      .2000   3.2963    .2000   3.2963 ! ALLOW ARO   ! His, adm jr., 9/4/89 
 NONBONDED  NY       .2000   3.2963    .2000   3.2963 ! ALLOW   ARO ! trp, JWK 
 NONBONDED  O        .1200   3.0291    .1200   2.4945 ! ALLOW   PEP POL ! This 1,4 vdW allows the C5 dipeptide minimum to exist.(LK) 
 NONBONDED  OB       .1200   3.0291    .1200   2.4945 ! ALLOW   PEP POL ARO  ! adm jr., 10/17/90, acetic acid carbonyl O 
 NONBONDED  OC       .1200   3.0291    .1200   3.0291 ! ALLOW   POL ION ! JG 8/27/89 
 NONBONDED  OH1      .1521   3.1538    .1521   3.1538 ! ALLOW   ALC ARO ! adm jr. 8/14/90, MeOH nonbond and solvent (same as TIP3P) 
 NONBONDED  OM       .1200   3.0291    .1200   3.0291 ! ALLOW HEM     ! Heme (6-liganded): CO ligand oxygen (KK 05/13/91) 
 NONBONDED  OS       .1521   3.1538    .1521   3.1538 ! ALLOW   ALC ARO  ! adm jr. 9/17/90, avoid O* wildcard 
 NONBONDED  OT       .1521   3.1506    .1521   3.1506 ! ALLOW   WAT   !TIP3P OXYGEN PARAMETERS, adm jr., NBFIX obsolete 
 NONBONDED  S        .4500   3.5636    .4500   3.5636 ! ALLOW   SUL ION  ! adm jr., 3/3/92, methanethiol/ethylmethylsulfide pure solvent 
 NONBONDED  SM       .3800   3.5190    .3800   3.5190 ! ALLOW  SUL  ION  ! adm jr., 3/3/92, dimethyldisulphide pure solvent 
 NONBONDED  SS       .4700   3.9200    .4700   3.9200 ! ALLOW  SUL       ! methylthiolate to water and F.E. of solvation, adm jr. 6/1/92 
 NONBONDED  ZN       .2500   1.9422    .2500   1.9422 ! ALLOW  ION       ! RHS March 18, 1990 
 NONBONDED  SOD      .1000   2.7297    .1000   2.7297 ! FROM BENOIT ROUX 
 NONBONDED  CLA      .1500   4.2763    .1500   4.2763 ! chloride 


! HBONDS 

! -----HYDROGEN-BOND-LIST-OPTIONS---------------------------
!  TOLERANCE=    .500                                     
!  DCUTOFF=  7.50 DOFF=  6.50 DON=  5.50         
!  ACUTOFF=100.00 AOFF= 80.00 AON= 60.00           
!  ACCEPTOR= TRUE                                              
   REXP=  6 AEXP=  4 HAEX=  0 AAEX=  0                    
! ----------------------------------------------------------
! END 

 ! HBOND               EMIN     RMIN

 HBOND  N*+* N%       .0000   3.0000  
 HBOND  O*   N*+*     .0000   2.9000  
 HBOND  OH*  N%       .0000   2.8500  
 HBOND  OH*  O*       .0000   2.7500  
 HBOND  S    N%       .0000   3.0000 ! ALLOW   SUL ION 
 HBOND  S    O*       .0000   2.9000 ! ALLOW   SUL ION 

 SET ECHO=TRUE END
