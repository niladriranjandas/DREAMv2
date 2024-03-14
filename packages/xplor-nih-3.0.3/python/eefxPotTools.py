"""
Tools to help create, manipulate and analyze <m eefxPot>.EEFxPot objects
which implements an implicit solvent energy for solution and membrane
environments, based on

 Lazaridis & Karplus Proteins 35, 133 (1999).
 Lazaridis  Proteins 52, 176 (2003).
"""

# Use with eefxPot/eefxPotTools for effective energy function of Xplor-NIH (eefxPot).
# Tian, Schwieters, Opella, Marassi (2013)
# [Tian et al. (2014) J Magn Reson 243:54]
# [Tian et al. (2015) Biophys J  109:574]
#
# Lazaridis-Karplus (LK) parameters for solvation free energy.
# [Lazaridis & Karplus (1999) Proteins 35:133)
# [Privalov & Makhatadze (1990) J Mol Biol 213:385]
# [Privalov & Makhatadze (1992) J Mol Biol 224:715]
# [Privalov & Makhatadze (1993) J Mol Biol 232:660]
# [Makhatadze & Privalov (1993) J Mol Biol 232:639]
#
# (dG-ref(chex) is from Lazaridis IMM parameters (2002):
# [Lazaridis (2003) Proteins 52:176]
#
# Suggested values (in angstrom) of membrane hydrophobic thickness for IMMx:
# (see Marsh, Handbook of lipid bilayers, 2nd Ed. p. 379)
# DMPC -----> eefx.setThickness(25.4)
# POPC -----> eefx.setThickness(27.0)
# DPPC -----> eefx.setThickness(28.6)
# DOPC -----> eefx.setThickness(29.6)
#
# Modified by Ye Tian, 06-13-2013.
# Modified by Ye Tian, 09-01-2016.


# EEFx parameters for charmm36.
# Columns:  1        2         3        4         5      6      7          8       9
#        Volume    dG-ref   dG-free   dH-ref   dCp-ref lamda   radius    dG-ref  dG-free
#                                                              (vdw)     (chex)   (chex)

param_LK36 = {
"C"    : (14.7	,   0.000 ,   0.00 ,   0.000,    0.00 , 3.50 ,  2.000 ,   0.000,  0.00 ),
"CT1"  : (23.7	,  -0.187 ,  -0.25 ,   0.876 ,   0.00 , 3.50 ,  2.275 ,  -0.645, -0.86 ),
"CT2"  : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.175 ,  -0.720, -1.01 ),
"CT3"  : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.060 ,  -0.665, -0.92 ),
"CC"   : (14.7	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  2.000 ,   0.000,  0.00 ),
"NH1"  : ( 4.4	,  -5.950 ,  -8.90 ,  -9.059 ,  -8.80 , 3.50 ,  1.850 ,  -1.145, -1.72 ),
"NR1"  : ( 4.4	,  -5.950 ,  -8.90 ,  -9.059 ,  -8.80 , 3.50 ,  1.850 ,  -1.145, -1.72 ),
"NY"   : ( 4.4	,  -5.950 ,  -8.90 ,  -9.059 ,  -8.80 , 3.50 ,  1.850 ,  -1.145, -1.72 ),
"NR2"  : ( 4.4	,  -3.820 ,  -4.00 ,  -4.654 ,  -8.80 , 3.50 ,  1.850 ,  -1.630, -1.71 ),
"NH2"  : (11.2	,  -5.450 ,  -7.80 ,  -9.028 ,  -7.00 , 3.50 ,  1.850 ,  -1.145, -1.64 ),
"NH3"  : (11.2	, -20.000 , -20.00 , -25.000 , -18.00 , 6.00 ,  1.850 ,  -1.145, -1.15 ),
"NC2"  : (11.2	, -10.000 , -10.00 , -12.000 ,  -7.00 , 6.00 ,  1.850 ,  -0.200, -0.20 ),
"N"    : ( 0.0	,  -1.000 ,  -1.55 ,  -1.250 ,   8.80 , 3.50 ,  1.850 ,  -1.145, -1.72 ),
"CPH1" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.800 ,  -0.410, -0.57 ),
"CPH2" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.800 ,  -0.410, -0.57 ),
"CA"   : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.992 ,  -0.410, -0.57 ),
"CY"   : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.992 ,  -0.410, -0.57 ),
"CPT"  : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.800 ,  -0.410, -0.57 ),
"CP1"  : (23.7	,  -0.187 ,  -0.25 ,   0.876 ,   0.00 , 3.50 ,  2.275 ,  -0.645, -0.86 ),
"CP2"  : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.175 ,  -0.720, -1.01 ),
"CP3"  : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.175 ,  -0.720, -1.01 ),
"OH1"  : (10.8	,  -5.920 ,  -6.70 ,  -9.264 , -11.20 , 3.50 ,  1.770 ,  -0.960, -1.09 ),
"O"    : (10.8	,  -5.330 ,  -5.85 ,  -5.787 ,  -8.80 , 3.50 ,  1.700 ,  -1.270, -1.39 ),
"OC"   : (10.8	, -10.000 , -10.00 , -12.000 ,  -9.40 , 6.00 ,  1.700 ,  -0.900, -0.90 ),
"S"    : (14.7	,  -3.240 ,  -4.10 ,  -4.475 , -39.90 , 3.50 ,  2.000 ,  -1.780, -2.25 ),
"SM"   : (14.7	,  -3.240 ,  -4.10 ,  -4.475 , -39.90 , 3.50 ,  1.975 ,  -1.780, -2.25 ),
"H"    : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.800 ,   0.000,  0.00 ),
"none" : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.800 ,   0.000,  0.00 )
}

# EEFx parameters for protein_eefx.top/par from 2014 [Tian et al. (2014) J Magn Reson 243:54].
# Columns:  1        2         3        4         5      6      7          8       9
#        Volume    dG-ref   dG-free   dH-ref   dCp-ref lamda   radius    dG-ref  dG-free
#                                                              (vdw)     (chex)   (chex)

param_LK = {
"C"    : (14.7	,   0.000 ,   0.00 ,   0.000,    0.00 , 3.50 ,  2.000 ,   0.000,  0.00 ),
"CR"   : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  1.800 ,  -1.350, -2.13 ),
"CH1E" : (23.7	,  -0.187 ,  -0.25 ,   0.876 ,   0.00 , 3.50 ,  2.275 ,  -0.645, -0.86 ),
"CH2E" : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.175 ,  -0.720, -1.01 ),
"CH3E" : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.060 ,  -0.665, -0.92 ),
"CR1E" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.992 ,  -0.410, -0.57 ),
"NH1"  : ( 4.4	,  -5.950 ,  -8.90 ,  -9.059 ,  -8.80 , 3.50 ,  1.850 ,  -1.145, -1.72 ),
"NR"   : ( 4.4	,  -3.820 ,  -4.00 ,  -4.654 ,  -8.80 , 3.50 ,  1.850 ,  -1.630, -1.71 ),
"NP"   : ( 4.4	,  -3.820 ,  -4.00 ,  -4.654 ,  -8.80 , 3.50 ,  1.850 ,   0.000,  0.00 ),
"NH2"  : (11.2	,  -5.450 ,  -7.80 ,  -9.028 ,  -7.00 , 3.50 ,  1.850 ,  -1.145, -1.64 ),
"NH3"  : (11.2	, -20.000 , -20.00 , -25.000 , -18.00 , 6.00 ,  1.850 ,  -1.145, -1.15 ),
"NC2"  : (11.2	, -10.000 , -10.00 , -12.000 ,  -7.00 , 6.00 ,  1.850 ,  -0.200, -0.20 ),
"N"    : ( 0.0	,  -1.000 ,  -1.55 ,  -1.250 ,   8.80 , 3.50 ,  1.850 ,  -1.145, -1.72 ),
"OM"   : (10.8	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  1.600 ,   0.000,  0.00 ),##
"OS"   : ( 0.0	,  -2.900 ,  -3.20 ,  -3.150 ,  -4.80 , 3.50 ,  1.600 ,  -0.960, -1.09 ),##
"OH1"  : (10.8	,  -5.920 ,  -6.70 ,  -9.264 , -11.20 , 3.50 ,  1.770 ,  -0.960, -1.09 ),
"O"    : (10.8	,  -5.330 ,  -5.85 ,  -5.787 ,  -8.80 , 3.50 ,  1.700 ,  -1.270, -1.39 ),
"OC"   : (10.8	, -10.000 , -10.00 , -12.000 ,  -9.40 , 6.00 ,  1.700 ,  -0.900, -0.90 ),
"OT"   : (14.0	,  -6.320 ,  -6.32 ,  -9.974 ,  12.00 , 3.50 ,  1.600 ,  -0.040, -0.04 ),##
"OH2"  : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 , 1.7398 ,   0.000,  0.00 ),##
"LP"   : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 , 0.2245 ,   0.000,  0.00 ),##
"FE"   : (15.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.650 ,   0.000,  0.00 ),##
"S"    : (14.7	,  -3.240 ,  -4.10 ,  -4.475 , -39.90 , 3.50 ,  2.000 ,  -1.780, -2.25 ),
"SM"   : (14.7	,  -3.240 ,  -4.10 ,  -4.475 , -39.90 , 3.50 ,  1.975 ,  -1.780, -2.25 ),
"SH1E" : (21.4	,  -2.050 ,  -2.70 ,  -4.475 , -39.90 , 3.50 ,  2.000 ,  -1.855, -2.44 ),
"H"    : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.800 ,   0.000,  0.00 ),
"none" : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.800 ,   0.000,  0.00 ),
"C5"   : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  1.800 ,  -1.350, -2.13 ),
"CF"   : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  1.800 ,  -1.350, -2.13 ),
"CW"   : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  1.800 ,  -1.350, -2.13 ),
"C5W"  : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  1.800 ,  -1.350, -2.13 ),
"CY"   : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  1.800 ,  -1.350, -2.13 ),
"CY2"  : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  1.800 ,  -1.350, -2.13 ),
"CH2G" : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.175 ,  -0.720, -1.01 ),
"CH2P" : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.175 ,  -0.720, -1.01 ),
"CR1W" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.992 ,  -0.410, -0.57 ),
"CR1H" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.992 ,  -0.410, -0.57 ),
"CRHH" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.992 ,  -0.410, -0.57 ),
"CT"   : (14.7	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  2.000 ,   0.000,  0.00 ),
"CM"   : (14.7	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  2.000 ,   0.000,  0.00 ),
"CS1"  : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  1.800 ,  -1.350, -2.13 ), ##CYSP
"CSN1" : (14.7	,   0.000 ,   0.00 ,   0.000,    0.00 , 3.50 ,  2.000 ,   0.000,  0.00 ), ##CYSP
"CSN2" : (14.7	,   0.000 ,   0.00 ,   0.000,    0.00 , 3.50 ,  2.000 ,   0.000,  0.00 ), ##CYSP
"C3S6" : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.060 ,  -0.665, -0.92 ), ##CYSP
"C3S7" : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.060 ,  -0.665, -0.92 ), ##CYSP
"C3S8" : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.060 ,  -0.665, -0.92 ), ##CYSP
"C3S9" : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.060 ,  -0.665, -0.92 ), ##CYSP
"CS2"  : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.992 ,  -0.410, -0.57 ), ##CYSP
"NO"   : ( 0.0	,  -1.000 ,  -1.55 ,  -1.250 ,   8.80 , 3.50 ,  1.850 ,  -1.145, -1.72 ), ##CYSP
"ON"   : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 , 1.7398 ,   0.000,  0.00 )  ##CYSP
}

# EEFx parameters for CHARMM19.
# Columns:  1        2         3        4         5      6      7          8       9
#        Volume    dG-ref   dG-free   dH-ref   dCp-ref lamda   radius    dG-ref  dG-free
#                                                              (vdw)     (chex)   (chex)
param_LK19 = {
"C"    : (14.7	,   0.000 ,   0.00 ,   0.000,    0.00 , 3.50 ,  2.100 ,   0.000,  0.00 ),
"CR"   : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  2.100 ,  -1.350, -2.13 ),
"CH1E" : (23.7	,  -0.187 ,  -0.25 ,   0.876 ,   0.00 , 3.50 ,  2.365 ,  -0.645, -0.86 ),
"CH2E" : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.235 ,  -0.720, -1.01 ),
"CH3E" : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.165 ,  -0.665, -0.92 ),
"CR1E" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  2.100 ,  -0.410, -0.57 ),
"NH1"  : ( 4.4	,  -5.950 ,  -8.90 ,  -9.059 ,  -8.80 , 3.50 ,  1.600 ,  -1.145, -1.72 ),
"NR"   : ( 4.4	,  -3.820 ,  -4.00 ,  -4.654 ,  -8.80 , 3.50 ,  1.600 ,  -1.630, -1.71 ),
"NP"   : ( 4.4	,  -3.820 ,  -4.00 ,  -4.654 ,  -8.80 , 3.50 ,  1.600 ,   0.000,  0.00 ),
"NH2"  : (11.2	,  -5.450 ,  -7.80 ,  -9.028 ,  -7.00 , 3.50 ,  1.600 ,  -1.145, -1.64 ),
"NH3"  : (11.2	, -20.000 , -20.00 , -25.000 , -18.00 , 6.00 ,  1.600 ,  -1.145, -1.15 ),
"NC2"  : (11.2	, -10.000 , -10.00 , -12.000 ,  -7.00 , 6.00 ,  1.600 ,  -0.200, -0.20 ),
"N"    : ( 0.0	,  -1.000 ,  -1.55 ,  -1.250 ,   8.80 , 3.50 ,  1.600 ,  -1.145, -1.72 ),
"OM"   : (10.8	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  1.600 ,   0.000,  0.00 ),
"OS"   : ( 0.0	,  -2.900 ,  -3.20 ,  -3.150 ,  -4.80 , 3.50 ,  1.600 ,  -0.960, -1.09 ),
"OH1"  : (10.8	,  -5.920 ,  -6.70 ,  -9.264 , -11.20 , 3.50 ,  1.600 ,  -0.960, -1.09 ),
"O"    : (10.8	,  -5.330 ,  -5.85 ,  -5.787 ,  -8.80 , 3.50 ,  1.600 ,  -1.270, -1.39 ),
"OC"   : (10.8	, -10.000 , -10.00 , -12.000 ,  -9.40 , 6.00 ,  1.600 ,  -0.900, -0.90 ),
"OT"   : (14.0	,  -6.320 ,  -6.32 ,  -9.974 ,  12.00 , 3.50 ,  1.600 ,  -0.040, -0.04 ),
"OH2"  : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 , 1.7398 ,   0.000,  0.00 ),
"LP"   : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 , 0.2245 ,   0.000,  0.00 ),
"FE"   : (15.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.650 ,   0.000,  0.00 ),
"S"    : (14.7	,  -3.240 ,  -4.10 ,  -4.475 , -39.90 , 3.50 ,  1.890 ,  -1.780, -2.25 ),
"SM"   : (14.7	,  -3.240 ,  -4.10 ,  -4.475 , -39.90 , 3.50 ,  1.890 ,  -1.780, -2.25 ),
"SH1E" : (21.4	,  -2.050 ,  -2.70 ,  -4.475 , -39.90 , 3.50 ,  1.890 ,  -1.855, -2.44 ),
"H"    : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.800 ,   0.000,  0.00 ),
"none" : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.800 ,   0.000,  0.00 ),
"C5"   : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  2.100 ,  -1.350, -2.13 ),
"CF"   : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  2.100 ,  -1.350, -2.13 ),
"CW"   : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  2.100 ,  -1.350, -2.13 ),
"C5W"  : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  2.100 ,  -1.350, -2.13 ),
"CY"   : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  2.100 ,  -1.350, -2.13 ),
"CY2"  : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  2.100 ,  -1.350, -2.13 ),
"CH2G" : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.235 ,  -0.720, -1.01 ),
"CH2P" : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.235 ,  -0.720, -1.01 ),
"CR1W" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  2.100 ,  -0.410, -0.57 ),
"CR1H" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  2.100 ,  -0.410, -0.57 ),
"CRHH" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  2.100 ,  -0.410, -0.57 ),
"CT"   : (14.7	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  2.100 ,   0.000,  0.00 ),
"CM"   : (14.7	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  2.100 ,   0.000,  0.00 )
}

# EEFx parameters for CHARMM22.
# Columns:  1        2         3        4         5      6      7          8       9
#        Volume    dG-ref   dG-free   dH-ref   dCp-ref lamda   radius    dG-ref  dG-free
#                                                              (vdw)     (chex)   (chex)
param_LK22 = {
"C"    : (14.7	,   0.000 ,   0.00 ,   0.000,    0.00 , 3.50 ,  2.0000 ,   0.000,  0.00 ),
"CT1"  : (23.7	,  -0.187 ,  -0.25 ,   0.876 ,   0.00 , 3.50 ,  2.2750 ,  -0.645, -0.86 ),
"CT2"  : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.1750 ,  -0.720, -1.01 ),
"CT3"  : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.0600 ,  -0.665, -0.92 ),
"CC"   : (14.7	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  2.0000 ,   0.000,  0.00 ),
"NH1"  : ( 4.4	,  -5.950 ,  -8.90 ,  -9.059 ,  -8.80 , 3.50 ,  1.8500 ,  -1.145, -1.72 ),
"NR1"  : ( 4.4	,  -5.950 ,  -8.90 ,  -9.059 ,  -8.80 , 3.50 ,  1.8500 ,  -1.145, -1.72 ),
"NY"   : ( 4.4	,  -5.950 ,  -8.90 ,  -9.059 ,  -8.80 , 3.50 ,  1.8500 ,  -1.145, -1.72 ),
"NR2"  : ( 4.4	,  -3.820 ,  -4.00 ,  -4.654 ,  -8.80 , 3.50 ,  1.8500 ,  -1.630, -1.71 ),
"NH2"  : (11.2	,  -5.450 ,  -7.80 ,  -9.028 ,  -7.00 , 3.50 ,  1.8500 ,  -1.145, -1.64 ),
"NH3"  : (11.2	, -20.000 , -20.00 , -25.000 , -18.00 , 6.00 ,  1.8500 ,  -1.145, -1.15 ),
"NC2"  : (11.2	, -10.000 , -10.00 , -12.000 ,  -7.00 , 6.00 ,  1.8500 ,  -0.200, -0.20 ),
"N"    : ( 0.0	,  -1.000 ,  -1.55 ,  -1.250 ,   8.80 , 3.50 ,  1.8500 ,  -1.145, -1.72 ),
"CPH1" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.8000 ,  -0.410, -0.57 ),
"CPH2" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.8000 ,  -0.410, -0.57 ),
"CA"   : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ),
"CY"   : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ),
"CPT"  : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.8000 ,  -0.410, -0.57 ),
"CP1"  : (23.7	,  -0.187 ,  -0.25 ,   0.876 ,   0.00 , 3.50 ,  2.2750 ,  -0.645, -0.86 ),
"CP2"  : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.1750 ,  -0.720, -1.01 ),
"CP3"  : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.1750 ,  -0.720, -1.01 ),
"OH1"  : (10.8	,  -5.920 ,  -6.70 ,  -9.264 , -11.20 , 3.50 ,  1.7700 ,  -0.960, -1.09 ),
"OT"   : (10.8	,  -5.920 ,  -6.70 ,  -9.264 , -11.20 , 3.50 ,  1.7700 ,  -0.960, -1.09 ),  ##for TIP3 water
"O"    : (10.8	,  -5.330 ,  -5.85 ,  -5.787 ,  -8.80 , 3.50 ,  1.7000 ,  -1.270, -1.39 ),
"OC"   : (10.8	, -10.000 , -10.00 , -12.000 ,  -9.40 , 6.00 ,  1.7000 ,  -0.900, -0.90 ),
"S"    : (14.7	,  -3.240 ,  -4.10 ,  -4.475 , -39.90 , 3.50 ,  2.0000 ,  -1.780, -2.25 ),
"SM"   : (14.7	,  -3.240 ,  -4.10 ,  -4.475 , -39.90 , 3.50 ,  1.9750 ,  -1.780, -2.25 ),
"H"    : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.8000 ,   0.000,  0.00 ),
"none" : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.8000 ,   0.000,  0.00 ),
## Residue LYR retinal schiff base.
"NCH1" : (11.2	, -10.000 , -10.00 , -12.000 ,  -7.00 , 6.00 ,  1.8500 ,  -0.200, -0.20 ), ##LYR
"CR15" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ), ##LYR
"CR14" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ), ##LYR
"CR13" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ), ##LYR
"CR12" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ), ##LYR
"CR11" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ), ##LYR
"CR10" : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ), ##LYR
"CR9"  : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ), ##LYR
"CR8"  : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ), ##LYR
"CR7"  : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ), ##LYR
"CR6"  : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ), ##LYR
"CR5"  : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.9924 ,  -0.410, -0.57 ), ##LYR
## Residue CYSP  [Sammalkorpi & Lazaridis (2007) Biophys J 92:10]
"CS1"  : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  1.800 ,  -1.350, -2.13 ), ##CYSP
"CSN1" : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  1.800 ,  -1.350, -2.13 ), ##CYSP
"CSN2" : ( 8.3	,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  1.800 ,  -1.350, -2.13 ), ##CYSP
"C3S6" : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.060 ,  -0.665, -0.92 ), ##CYSP
"C3S7" : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.060 ,  -0.665, -0.92 ), ##CYSP
"C3S8" : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.060 ,  -0.665, -0.92 ), ##CYSP
"C3S9" : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.060 ,  -0.665, -0.92 ), ##CYSP
"CS2"  : (18.4	,   0.057 ,   0.08 ,  -0.973 ,   6.90 , 3.50 ,  1.992 ,  -0.410, -0.57 ), ##CYSP
"NO"   : ( 4.4	,  -1.000 ,  -1.55 ,  -1.250 ,   8.80 , 3.50 ,  1.850 ,  -1.145, -1.72 ), ##CYSP
"ON"   : ( 10.8	,  -4.337 ,  -4.74 ,  -4.687 ,  -7.13 , 3.50 ,  1.739 ,  -1.270, -1.39 )  ##CYSP
}

#
# EEFx-Lipid parameters for CHARMM22. 
# Columns:  1        2         3        4         5      6      7          8       9
#        Volume    dG-ref   dG-free   dH-ref   dCp-ref lamda   radius    dG-ref  dG-free
#                                                              (vdw)     (chex)   (chex)
param_LK22_LIPID = {
"CL"   : (14.7	,   0.000 ,   0.00 ,   0.000,    0.00 , 3.50 ,  2.000 ,   0.000,  0.00 ),
"CTL1" : (23.7	,  -0.187 ,  -0.25 ,   0.876 ,   0.00 , 3.50 ,  2.275 ,  -0.645, -0.86 ),
"CTL2" : (22.4	,   0.372 ,   0.52 ,  -0.610 ,  18.60 , 3.50 ,  2.175 ,  -0.720, -1.01 ),
"CTL3" : (30.0	,   1.089 ,   1.50 ,  -1.779 ,  35.60 , 3.50 ,  2.060 ,  -0.665, -0.92 ),
"OBL"  : (10.8	,  -5.330 ,  -5.85 ,  -5.787 ,  -8.80 , 3.50 ,  1.700 ,  -1.270, -1.39 ),
"O2L"  : (10.8	, -10.000 ,   0.00 , -12.000 ,  -9.40 , 6.00 ,  1.700 ,  -0.900, -0.90 ),
"OHL"  : (10.8	,  -5.920 ,  -6.70 ,  -9.264 , -11.20 , 3.50 ,  1.770 ,  -0.960, -1.09 ),
"OSL"  : (10.8	,  -5.920 ,  -6.70 ,  -9.264 , -11.20 , 3.50 ,  1.770 ,  -0.960, -1.09 ),
"OSLP" : (10.8	,  -5.920 ,   0.00 ,  -9.264 , -11.20 , 3.50 ,  1.770 ,  -0.960, -1.09 ),
"NH3L" : (11.2	, -20.000 , -20.00 , -25.000 , -18.00 , 6.00 ,  1.850 ,  -1.145, -1.15 ),
"NTL"  : (11.2	, -13.000 , -26.00 , -25.000 , -18.00 , 6.00 ,  1.850 ,  -1.145, -1.15 ),
"SL"   : (14.7	,  -3.240 ,  -4.10 ,  -4.475 , -39.90 , 3.50 ,  1.975 ,  -1.780, -2.25 ),
"PL"   : (17.6	, -13.000 , -19.00 , -25.000 , -18.00 , 6.00 ,  2.150 ,  -1.145, -1.15 ),  #calculated by YT
"none" : ( 0.0	,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.800 ,   0.000,  0.00 )
}

# MOPED (Method for Optimizing Physical Energy Parameters Using Decoys) parameters.
# [Seok et. al (2002) J Comp Chem 24:89]
#
# Modified by Ye Tian, 06-13-2013.
#
# Columns:  1        2         3        4         5      6      7          8       9
#        Volume    dG-ref   dG-free   dH-ref   dCp-ref lamda   radius    dG-ref  dG-free
#                                                              (vdw)     (chex)   (chex)
param_MOPED = {
"C"    : (15.4 ,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  2.100 ,  0.000 , 0.000 ),
"CR"   : ( 8.3 ,  -0.890 ,  -1.40 ,   2.220 ,   6.90 , 3.50 ,  2.150 , -1.350 , -2.12 ),
"CH1E" : (25.8 ,  -0.187 ,  -0.09 ,   0.876 ,   0.00 , 4.62 ,  2.265 , -0.645 , -0.31 ),
"CH2E" : (17.8 ,   0.372 ,   4.52 ,  -0.610 ,  18.60 , 4.86 ,  2.635 , -0.720 , -8.75 ),
"CH3E" : (30.0 ,   1.089 ,   3.05 ,  -1.779 ,  35.60 , 3.31 ,  2.574 , -0.665 , -1.86 ),
"CR1E" : (25.8 ,   0.057 ,  -1.44 ,  -0.973 ,   6.90 , 3.16 ,  2.000 , -0.410 , 10.37 ),
"NH1"  : ( 6.5 ,  -5.950 ,  -8.34 ,  -9.059 ,  -8.80 , 2.97 ,  1.500 , -1.145 , -1.61 ),
"NR"   : ( 4.4 ,  -3.820 ,  -6.05 ,  -4.654 ,  -8.80 , 3.60 ,  1.950 , -1.630 , -2.58 ),
"NH2"  : (11.4 ,  -5.450 ,  -7.83 ,  -9.028 ,  -7.00 , 3.20 ,  1.551 , -1.145 , -1.65 ),
"NH3"  : (12.0 , -20.000 , -16.89 , -25.000 , -18.00 , 4.17 ,  1.250 , -1.145 , -0.97 ),
"NC2"  : (10.4 , -10.000 ,  -8.99 , -12.000 ,  -7.00 , 3.59 ,  1.100 , -0.200 , -0.18 ),
"N"    : ( 0.0 ,  -1.000 ,  -1.45 ,  -1.250 ,   8.80 , 2.62 ,  1.550 , -1.145 , -1.66 ),
"OH1"  : (12.5 ,  -5.920 , -12.41 ,  -9.264 , -11.20 , 3.94 ,  2.300 , -0.960 , -2.01 ),
"O"    : (10.2 ,  -5.330 ,  -8.69 ,  -5.787 ,  -8.80 , 2.92 ,  1.685 , -1.270 , -2.07 ),
"OC"   : ( 7.4 , -10.000 , -10.72 , -12.000 ,  -9.40 , 4.45 ,  1.800 , -0.900 , -0.97 ),
"S"    : (14.7 ,  -3.240 ,  -5.01 ,  -4.475 , -39.90 , 2.62 ,  1.440 , -1.780 , -2.76 ),
"SH1E" : (21.5 ,  -2.050 ,  -4.43 ,  -4.475 , -39.90 , 3.00 ,  2.090 , -1.855 , -4.01 ),
"H"    : ( 0.0 ,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.800 ,  0.000 ,  0.00 ),
"none" : ( 0.0 ,   0.000 ,   0.00 ,   0.000 ,   0.00 , 3.50 ,  0.800 ,  0.000 ,  0.00 )
}


useTorsionDB=True
def initEEFx(topparSet=None,
             ):
    """
    Initialize topology and parameters for use with the EEFx implicit solvent
    force field. If present, topparSet specifies a two membered tuple containing
    the names of topology and parameter files.


    Finally, if topparSet is not specified and the EEFX_TOPPAR environment
    variable is present, it is interpreted as either a special name, or the
    topology and parameter names separated by the colon character.

    Alternately, topparSet can be one of the special names: "default" or
    "charmm22".

    If the parameter file sets the XPLOR variable protein_par_vers to a string
    containing "charmm22", the module-local variable useTorsionDB is set to
    False. Otherwise, this value takes the value True.
    
    """
    if not topparSet:
        from os import environ as env
        if "EEFX_TOPPAR" in list(env.keys()):
            val = env["EEFX_TOPPAR"]
            if ":" in val:
                topparSet = val.split(":")
            else:
                topparSet = val
                pass
            pass
        else:
            topparSet="default"
            pass
        pass

    if topparSet=="default":
        topparSet = ("protein_eef.top",
                     "protein_eef.par")
    elif topparSet=="charmm22":
        topparSet = ("eefx/charmm22/topallh22x-eefx.pro",
                     "eefx/charmm22/parallh22x-eefx.pro")
    elif type(topparSet)==type("string"):
        raise Exception("Invalid value for topparSet.")

    
    
    import protocol
    protocol.topology['protein']  = topparSet[0]
    protocol.parameters['protein']= topparSet[1]
    protocol.initParams('protein')

    if "charmm22" in protocol.parVersion['protein']:
        useTorsionDB=False
        pass

    return


def create_EEFxPot(name,
                   selection="not PSEUDO",
                   immx=False,
                   thickness=20,
                   profileN=10,
                   a=0.85,
                   paramSet=None,
                   moveTol=0.5,
                   temperature=298.15,
                   potList=None,
                   selPairs=[],
                   asList=False,
                   verbose=False):
    """
    Create an instance of <m eefxPot>.EEFxPot, with moveTol and temperature
    settings as specified in the argument list.

    paramSet specifies a solvation parameter set to use. If it is not
    specified, a value is deduced from the XPLOR variable
    protein_par_vers specified in the parameter file. If set
    explicitly, valid values in this module are:

       param_LK36
       param_LK
       param_LK19
       param_LK22
       param_LK22_LIPID
       param_MOPED

    This value depends on the XPLOR protein_par_vers variable, which should be 
    set in the .par parameter file. If protein_par_vers has a value of eef-1.1,
    eef-1.2, eef-2.0, or contains the string "charmm22" paramsSet is set to
    param_LK22. Otherwise, it is set to param_LK.

    If the potList argument is given, separate terms will be added to
    the specified potList for electrostatic, solvation, and VdW
    terms. In this case, the returned potential term should not be
    appended to this potList.

    The selPairs argument is an optional sequence of pairs of <m
    atomSel>.AtomSel objects or strings selection strings which
    specify particular pairs of atoms to compute interactions between.

    If asList is specified, a PotList object will be returned
    containing entries for the electrostatic, solvation, and VdW
    terms.

    Setting immx to True enables IMMx membrane mode, where set following
    parameters are also set and have meaning:

         thickness - membrane thickness
         profileN  - membrane profile parameter
         a         - membrane dielectric screening parameter
         

  """

    from selectTools import convertToAtomSel
    sel = convertToAtomSel(selection)

    sim = sel.simulation()
    if len(selPairs):
        from atomSel import union, AtomSel
        pairs=[]
        for sel1,sel2 in selPairs:
            sel1 = convertToAtomSel(sel1,sim)
            sel2 = convertToAtomSel(sel2,sim)
            pairs.append( ( sel1 , sel2 ) ) 
            pass
        pass


    import protocol
    if not paramSet:
        if  protocol.topVersion['protein']=="eef-1.1" or \
           protocol.topVersion['protein']=="eef-1.2" or \
           protocol.topVersion['protein']=="eef-2.0" or \
           "charmm22" in protocol.parVersion['protein']:
            paramSet=param_LK22
        else:
            paramSet=param_LK
            pass
        pass
    


    groupDefs=[]
    paramkeys=list(paramSet.keys())

    from atomSel import AtomSel
    for atom in AtomSel("all",sel.simulation()):
        if sel.containsAtom(atom):
            try:
                key = atom.chemType()
                if key not in paramkeys:
                    if not key.startswith("H"):
                        print("Warning: create_EEFxPot: ", end=' ')
                        print("chem type not present in parameter set:", key)
                        pass
                    groupDefs.append("none")
                else:
                    groupDefs.append( key )
                    pass
                pass
            except IndexError:
                raise Exception("Atom %s contains unknown residue or atom type" %
                                atom.string())
            pass
        else:
            groupDefs.append("none")
            pass
        pass
            
    from nonBondTools import getNBParamsFromXplor
    params = getNBParamsFromXplor(sim)


##=======
    #now do something with nbfixParams: they are  A, B, A14, B14
    #indexed by the string T1*T2, where T1 and T2 are chemical types
    #nonbonded parameters are generated using combination rules from
    # entries given by NONBon specifications, except for those pairs
    # explicitly specified by NBFIx    

##====================
                               
    if verbose:
        print("Using %d atoms in %d groups" % (sel.size(),
                                               len(list(paramSet.keys()))))
    from eefxPot import EEFxPot
    pot = EEFxPot(name,groupDefs,paramSet,
                  params.startIndices,
                  params.nbParams,
                  params.resExclude,
                  params.nbfixParams,
                  sel)

    if immx:
        pot.setIMMx(immx)
        pot.setThickness(20)
        pot.setProfileN(10)
        pot.setA(0.85)
        pass
    pot.setVerbose(verbose)
    pot.setMoveTol(moveTol)
    pot.setTemperature(temperature)

    if ("eef-1.2" in protocol.parVersion['protein'] or 
        "eef-2.0" in protocol.parVersion['protein'] or 
        "charmm22" in protocol.parVersion['protein'] or 
        "charmm36" in protocol.parVersion['protein']  ):
        pot.setE14FullScale(True)
        pass

    if potList!=None:
        potList.append( pot.elecPot() )
        potList.append( pot.solvPot() )
        potList.append( pot.vdwPot()  )
        pass

    if asList:
        base = pot
        pot = PotProxy(EEFx_PotList(name))
        base.resetInstanceName(name+"-base")
        for term in (base.elecPot(),base.solvPot(), base.vdwPot()):
            pot.append(term)
        pass

    if selPairs:
        pot.deleteSelectionPair(0)
        for pair in pairs:
            pot.addSelectionPair( *pair )
            #pot.savedSelPairs.append( pair )
            pass
        pass

    return pot

import potList
class EEFx_PotList(potList.realPotList):
    """Object type returned if the asList argument is True in create_EEFxPot above.
    """
    def __init__(self,*args):
        potList.realPotList.__init__(self,*args)
        return
    def __getattr__(self, name):
        print("arg: " + name)
        if name=="scale" or name=="setScale":
            return getattr(self[0].baseRef, name)
            
        try:
            return potList.realPotList.__getattribute__(self, name)
        except:
            #print "base class method: " + name
            try:
                return getattr(self[0].baseRef, name)
            except:
                raise Exception("EEFx_PotList: could not find attribute: " + name)
        pass
    pass
    

def EEFx_Shadow__getattr__(self, name):
    try:
        return getattr(self, name)
    except:
        return getattr(self[0].baseRef, name)
    pass


def eefxAtomInfo(pot="",sel=""):
    """Needs to be documented appropriate.
    """
    from atomSel import AtomSel
    try:
        atoms=AtomSel(sel)
        if len(atoms)!=0:
            for i in atoms:
                pot.printAtomInfo(i.index())
        else:
            print("No atoms selected.")
    except:
        print("Error reporting atom infomation for EEFx.")
    
def eefxResEnergy(pot="",sel="",output=""):
    """Needs to be documented appropriate.
    """
    from atomSel import AtomSel
    try:
        res=[]
        resnames=[]
        for i in AtomSel(sel):
            if i.residueNum() not in res:
                res.append(i.residueNum())
                resnames.append(i.residueName())
        allE=[]
        allAtomE=pot.atomEnergies()
        for r in res:
            resE={'Evdw':0.,'Eelec':0.,'Eslv':0.,'Egc':0.,'Evolt':0.}
            for j in AtomSel("resid "+str(r)):
                resE['Evdw'] += allAtomE[j.index()].Evdw
                resE['Eelec']+= allAtomE[j.index()].Eelec
                resE['Eslv'] += allAtomE[j.index()].Eslv
                resE['Egc']  += allAtomE[j.index()].Egc
                resE['Evolt']+= allAtomE[j.index()].Evolt
            allE.append(resE)
        if len(output)!=0:
            try:
                fo=open(output,"w")
                fo.write('Res\tEvdw\tEelec\tEslv\tEgc\tEvolt\n')
                for i in range(len(res)):
                    fo.write(resnames[i]+str(res[i])+'\t')
                    fo.write("%.2f"%allE[i]['Evdw']+'\t'+"%.2f"%allE[i]['Eelec']+'\t'+"%.2f"%allE[i]['Eslv']+'\t'+"%.2f"%allE[i]['Egc']+'\t'+"%.2f"%allE[i]['Evolt']+'\n')
                fo.close()
            except:
                print("Can't export EEFx energies to file.")
        else:
            print('Res\tEvdw\tEelec\tEslv\tEgc\tEvolt')
            for i in range(len(res)):
                print(resnames[i]+str(res[i])+'\t',"%.2f"%allE[i]['Evdw'],'\t',"%.2f"%allE[i]['Eelec'],'\t',"%.2f"%allE[i]['Eslv'],'\t',"%.2f"%allE[i]['Egc'],'\t',"%.2f"%allE[i]['Evolt'])

    except:
        print("Error reporting EEFx energy by residue.")



def eefxAtomEnergy(pot="",output=""):
    """Needs to be documented appropriate.
    """
    from atomSel import AtomSel
    try:
        allE=[]
        allAtomE=pot.atomEnergies()
        print('Res\tAtom\tType\tEvdw\tEelec\tEslv\tEgc\tEvolt')
        for i in AtomSel("ALL"):
            AtomE={'Evdw':0.,'Eelec':0.,'Eslv':0.,'Egc':0.,'Evolt':0.}
            AtomE['Evdw'] += allAtomE[i.index()].Evdw
            AtomE['Eelec']+= allAtomE[i.index()].Eelec
            AtomE['Eslv'] += allAtomE[i.index()].Eslv
            AtomE['Egc']  += allAtomE[i.index()].Egc
            AtomE['Evolt']+= allAtomE[i.index()].Evolt

            print(i.residueName()+str(i.residueNum())+'\t'+i.atomName()+'\t'+i.chemType()+'\t',"%.3f"%AtomE['Evdw'],'\t',"%.3f"%AtomE['Eelec'],'\t',"%.3f"%AtomE['Eslv'],'\t',"%.3f"%AtomE['Egc'],'\t',"%.3f"%AtomE['Evolt'])
            
            allE.append(AtomE)
    except:
        print("Error reporting EEFx energy by residue.")


def create_EEFxBoundaryPot(name='EEFxWall',
                           Xsize=None,Ysize=None,Zsize=None,
                           radius=0):
    """Create a boundary potential term to restrain atoms from flying off.

    If radius>0 a spherical boundary will be created. Otherwise,
    <m planeDistPot> plane-distance potentials will be created at +/-
    values given in the Xsize, Ysize, and Zsize arguments. Any of
    these set to None will cause no potential in that direction.
    """
    ret=None

    from atomSel import AtomSel
    sel = AtomSel("all and (not pseudo) and (not name H*)")
    restraints=""
    if radius>0:
        from noePotTools import create_NOEPot
        for a in sel:
            restraints += "\n" + \
                          "assign (all and not pseudo and not name H*) (resid " +\
                          +str(a.residueNum())+" and name "+ \
                          a.atomName()+") "+str(radius)+" "+str(radius)+" 0"
            pass
        ret = create_NOEPot(name,restraints=restraints)
        ret.setAveType('center')
        pass
    else:
        from planeDistTools import create_PlaneDistPot
        from potList import PotList
        ret = PotList(name) 
        if Xsize!=None:
            restraints=""
            for a in sel: 
                restraints+="assign (resid "+str(a.residueNum())+" and name "+\
                             a.atomName()+") "+str(Xsize)+" " +\
                             str(Xsize)+" 0"+ "\n"
                pass
            xWall1=create_PlaneDistPot("xMinus",A=1.,B=0,C=0,D=-Xsize/2.,
                                       restraints=restraints)
            xWall2=create_PlaneDistPot("xPlus",A=1.,B=0,C=0,D= Xsize/2.,
                                       restraints=restraints)
            xWall1.setFreedom("fix")
            xWall2.setFreedom("fix")
            ret.append(xWall1)
            ret.append(xWall2)
            pass
        if Ysize!=None:
            restraints=""
            for a in sel: 
                restraints+= "\n" +"assign (resid "+str(a.residueNum())+\
                             " and name "+a.atomName()+") "+str(Ysize)+" "+\
                             str(Ysize)+" 0"
                pass
            yWall1=create_PlaneDistPot("yMinus",A=0,B=1.,C=0,D=-Ysize/2.,
                                       restraints=restraints)
            yWall2=create_PlaneDistPot("yPlus",A=0,B=1.,C=0,D= Ysize/2.,
                                       restraints=restraints)
            yWall1.setFreedom("fix")
            yWall2.setFreedom("fix")
            ret.append(yWall1)
            ret.append(yWall2)
            pass
        if Zsize!=None:
            restraints=""
            for a in sel: 
                restraints+= "\n" +"assign (resid "+str(a.residueNum())+\
                             " and name "+a.atomName()+") "+str(Zsize)+" "+\
                             str(Zsize)+" 0"
                pass
            zWall1=create_PlaneDistPot("zMinus",A=0,B=0,C=1.,D=-Zsize/2.,
                                       restraints=restraints)
            zWall2=create_PlaneDistPot("zPlus" ,A=0,B=0,C=1.,D= Zsize/2.,
                                       restraints=restraints)
            zWall1.setFreedom("fix")
            zWall2.setFreedom("fix")
            ret.append(zWall1)
            ret.append(zWall2)
            pass
        pass
    return ret




#for backward compatibility
from membraneTools import setCenter, setCenterXY, flipAll
