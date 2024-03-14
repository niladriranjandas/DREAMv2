Information about usage of the "average3d.py" module:

- launch PyMOL as usual

- load your structure as usual (e.g., type "load 1T3V.pdb, 1t3v" at the pymol prompt)

- import the average3d.py module by typing "run /some/path/to/average3d.py"

- use the avgStates() function to average the conformers, calculate RMSDs,
  write data to files, etc... See the file "1T3V_average3d_example.pml.py" 
  for examples of usage of avgStates() with the 1T3V structure, and see the
  file "1C89_average3d_example.pml.py" for examples with the 1C89 structure.

- after importing "average3d.py" into pymol, you may also type help(avgStates) 
  for info about the syntax of the avgStates() function... (or just see the file).
