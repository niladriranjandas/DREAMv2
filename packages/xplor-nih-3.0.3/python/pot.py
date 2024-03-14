
#present here for documentation purposes only - no code

"""
Pot, a generic energy term.

All Xplor-NIH energy terms derive from this object. It contains the following
member functions:

  instanceName() - name given to the energy term when it is created.
  potName()      - the type of energy term

  calcEnergy()                 - calculate energy, returns the energy value.
  calcEnergyAndDerivs(derivs)  - calculate energy, derivatives, 
  			         returns the energy value.

  scale()         - energy scale factor
  setScale(val)   - set this value.

  numRestraints()  - return the number of restraints defined for this term.
  		     Returns -1 if no reasonable value can be assigend.
  violations()     - return number of violated restraints. Returns -1 if no
  		     reasonable value can be assigend.
  rms()            - return root mean square deviation for this term. Returns -1
  		     if no reasonable value can be assigend.

  threshold()       - threshold deviation above which a restraint is considered
  		      violated.
  setThreshold(val) - set this value.
"""

