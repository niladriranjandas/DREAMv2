# DREAMv2
Distance restraint and energy-assisted modeling from NMR data


Packages needed:

	GROMACS (version:2023) link:https://manual.gromacs.org/

	MODELLER (version 10.4) link:https://salilab.org/modeller/10.4/release.html
	(free for academic usage)	

	DSSP (version 2.0.4) https://ssbio.readthedocs.io/en/latest/instructions/dssp.html


#############################################################################################

Set the following environment variables.

export DSSP="dssp"
export GMX="/usr/local/gromacs/bin/gmx"
export MOD="mod10.4"


############################################################################################

Required files:
	.seq : amino acid sequence
	.upl : experimental upper bounds 
	.aco : diherdal angle bounds (if present)


DREAMv2 can either be run in the "expert" mode or the "normal" mode. 
	1. Expert mode: Parameters required by DREAMv2 are supplied the user. These include,
		aug: 		Augment upper bounds with triangle inequality if sum is less than "aug"
		etalo:		Lower bound for the density of fragments in divide stage
		inclneigh:	Expand residue groups to include any other residue if it shares at least "inclneigh" number of upper bounds with the group in total
		grpexp:		Expand atom groups to include atoms having marginal upper bounds with the group

	2. Normal mode: The predictive module of DREAMv2 estimates the required parameters.



1. Preliminary preparation before running the structure computation module
	1.1 For expert mode, run the command:

	./setPreliminary_withparams.sh --name <protein_name> --seq <seq_file> --upl <upl_file> --ang <dihedral_file> --aug <augment_flag> --etalo <subgraph_density> --inclneigh <include_neighbour> --grpexp <group_expand>

	e.g.:
	./setPreliminary_withparams.sh --name 2m4k --seq ../protein/2m4k_param_2m4k/2m4k.seq --upl ../protein/2m4k_param_2m4k/2m4k_concat_dist.upl --ang ../protein/2m4k_param_2m4k/2m4k_concat_dihed.aco --aug 0 --etalo 0.7 --inclneigh 8 --grpexp 30


	1.2 For normal mode, run the command:
	./setPreliminary_withparams.sh --name <protein_name> --seq <seq_file> --upl <upl_file> --ang <dihedral_file>

	e.g.:
	./setPreliminary_withparams.sh -n 2m4k -s ../protein/2m4k_param_2m4k/2m4k.seq -u ../protein/2m4k_param_2m4k/2m4k_concat_dist.upl -a ../protein/2m4k_param_2m4k/2m4k_concat_dihed.aco


	The above steps,
		1. Generates a random number ("random_No.") and appends to name of the protein ("name") given by the user
		2. Creates a folder directory DREAMv2/protein/<name>_<random_No.>
		3. Creates an xml file (DREAMv2/protein/<name>_<random_No.>/<name>_<random_No.>.xml) containing all the inputs. The user is requested to check this file in all the input names are in order. The user may refer to any ".xml" in the existing folders in "DREAMv2/protein".

2. Run the code:
	./run_test.sh <name>_<random_No.>


3. Viewing the results:
	The final output files are in the folder DREAMv2/protein/finalpdbfiles. They consist of:
	1. Protein coordinate files (.pdb files)
	2. Experimental distance violations and Ramachandran plots corresponding to each of the protein structure coordinate files (.png files).


Running DREAMv2 for already tested protein instances
	The required input files are already present in the folder. DREAMv2 can be run for the already tested examples by running the command:

	./run_test.sh <protein_name>
	Here protein_name can be 2m4k, 2dk6, 2cqo, 1t17, 1v9w, 1pbu, 5x1x, 2d9b, 2cpr, or 1xxe



