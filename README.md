# DREAMv2
Distance restraint and energy-assisted modeling from NMR data

# Packages needed:
-	GROMACS (version:2023)  [https://manual.gromacs.org/]
-	MODELLER (version 10.4) [https://salilab.org/modeller/10.4/release.html] (free for academic usage)	
-	DSSP (version 2.0.4) [https://ssbio.readthedocs.io/en/latest/instructions/dssp.html]
-	MATLAB (version R2019a)
-	Pymol  (version 2.5.0) [https://storage.googleapis.com/pymol-storage/installers/index.html]

The following Python (version 3.8.5) packages are needed:
-	Numpy
-	Pandas
-	Plotly
-	matplotlib

The following packages are already included:
-	XPLOR (version 3.0.3): for thin layer water-refinement in the post-processing stages
-	PdbStat: for format conversion of the .upl file
-	Reducer: for adding hydrogen in the gap-filling stage
-	SDPT3: for computing solution of the underlying semidefinite programming convex optimization formulation
-	CVX: for calling SPDT3

Set the following environment variables.

    export DSSP=<your path to the DSSP command>
    export GMX=<path to GROMACS gmx command>
    export MOD=<path to the MODELLER command>


# Required files and modes:

-   .seq file:  amino acid sequence
-   .upl file:  experimental upper bounds 
-   .aco file:  diherdal angle bounds (if present)


For user-specified input files, DREAMv2 can run either in the "expert" mode or the "normal" mode. 
-   Expert mode: Parameters required by DREAMv2 are supplied by the user. These include,
	-   aug: 	Augment upper bounds with the triangle inequality if the sum is less than "aug".
	-   etalo:	Lower bound of the density of fragments in the divide stage.
	-   inclneigh:	Expand residue groups to include any other residue if it shares at least "inclneigh" number of upper bounds with the group in total.
	-   grpexp:	Expand atom groups to include atoms having marginal upper bounds with the group.

-   Normal mode: The predictive module of DREAMv2 estimates the required parameters.



### Preliminary preparation before running the structure computation module
Go to the folder DREAMv2/code
-	For expert mode, run the command:

./setPreliminary_withparams.sh --name <protein_name> --seq <seq_file> --upl <upl_file> --ang <dihedral_file> --aug <augment_flag> --etalo <subgraph_density> --inclneigh <include_neighbour> --grpexp <group_expand>

E.g.,

	cd DREAMv2/code
 
	./setPreliminary_withparams.sh --name 2m4k --seq ../protein/2m4k_param_2m4k/2m4k.seq --upl ../protein/2m4k_param_2m4k/2m4k_concat_dist.upl --ang ../protein/2m4k_param_2m4k/2m4k_concat_dihed.aco --aug 0 --etalo 0.7 --inclneigh 8 --grpexp 30


-   For normal mode, run the command:

./setPreliminary_withparams.sh --name <protein_name> --seq <seq_file> --upl <upl_file> --ang <dihedral_file>

E.g.,

	cd DREAMv2/code
 
	./setPreliminary_withparams.sh --name 2m4k --seq ../protein/2m4k_param_2m4k/2m4k.seq --upl ../protein/2m4k_param_2m4k/2m4k_concat_dist.upl --ang ../protein/2m4k_param_2m4k/2m4k_concat_dihed.aco


The above steps,
-   Generates a random number ("random No.") and appends it to the name of the protein ("name") given by the user
-   Creates a folder directory DREAMv2/protein/<"name"\_"random No.">
-   Creates an xml file (DREAMv2/protein/<"name"\_"random No.">/<"name"\_"random No.">.xml) containing all the inputs. The user is requested to check this file to see that all the input names are in order. The user may refer to any ".xml" in the existing folders in "DREAMv2/protein".

### Run the code:
Go to the folder DREAMv2/code

	./run_test.sh <"name"_"random No.">

Replace name with the protein name supplied and random No. with the one generated in the preliminary step

### Viewing the results:
The final output files are in the folder DREAMv2/protein/finalpdbfiles. They consist of:
-   Protein coordinate files (.pdb files)
-   Experimental distance violations and Ramachandran plots corresponding to each of the protein structure coordinate files (.png files).

### Running DREAMv2 for already tested protein instances
The required input files are already present in the folder. DREAMv2 can be run for the already tested examples by running the command:

	./run_test.sh <protein_name>

Here "protein_name" can be 2m4k, 2dk6, 2cqo, 1t17, 1v9w, 1pbu, 5x1x, 2d9b, 2cpr, or 1xxe.




