Example of generation of a potential file.  
------------------------------------------

Specifically, this directory contains tools to recreate ../rna09_v0.dat.


File: suites_chi_res3.0_b60_allringatoms_epsilonfilter.csv
The RNA09 structrue database (excluding alternate conformations and clashing atoms) was
used to extract parameters (e.g., torsion angles, B-factors) of suites with:
* resolution of 3.0 Angstroms of better, 
* all atoms with B-factors less or equal to 60,
* no suspicious sugar puckers,
* 155 <= epsilon < -50 deg.
The results are in this file.


Each directory, except 'overall', contains the tools to generate the different
terms of the potential.  'overall' has a script that collects all the individual
tems into one file: rna09_v0.dat.






