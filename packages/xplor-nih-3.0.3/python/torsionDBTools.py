"""Module for extracting torsion angle information from a structure database.

Additional functionality is provided to estimate probability density using the
histogram approach.


When using this module, please cite:

Bermejo, G.A., Clore, G.M., and Schwieters, C.D. (2012). Smooth statistical
torsion angle potential derived from a large conformational database via
adaptive kernel density estimation improves the quality of NMR protein
structures. Protein Sci. 21, 1824-1836.

"""
# Written by Guillermo A. Bermejo, September 2009.
# Based on preliminary code by John Kuszewski.

import math
import itertools

import psfGen
import simulation
import atomSel
import dihedral
import selectTools
import protocol

#______________________________________________________________________________
# This section contains the torsion angle definitions for each residue.
# Eventually this could be merged with torsionTools module.

# torsion_info is the data structure that contains all the corresponding torsion
# angle information.

# The order in the tuple containing the atom selections is important (hence the
# use of tuple and not list) since they define the torsion around the bond
# between the atoms selected by the two central selections.

torsion_info = {'ALA': {},
                'GLY': {},
                'SER': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name og)'),
                                 (-180.0, 180.0))},
                'THR': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name og1)'),
                                 (-180.0, 180.0))},
                'VAL': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg1)'),
                                 (-180.0, 180.0))},
                'ASN': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name od1)'),
                                 (-180.0, 180.0))},
                'ASP': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name od1)'),
                                 (-90.0, 90.0))}, # compatible with selectTools
                'HIS': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name nd1)'),
                                 (-180.0, 180.0))},
                'ILE': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg1)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg1)',
                                  '(segid "%s" and resid %d and name cd*)'),
                                 (-180.0, 180.0))},
                'LEU': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd1)'),
                                 (-180.0, 180.0))},
                'PHE': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd1)'),
                                 (0.0, 180.0))}, # compatible with selectTools
                'TRP': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd1)'),
                                 (-180.0, 180.0))},
                'TYR': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd1)'),
                                 (0.0, 180.0))}, # compatible with selectTools
                'GLU': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)'),
                                 (-180.0, 180.0)),
                        'chi3': (('(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)',
                                  '(segid "%s" and resid %d and name oe1)'),
                                 (-90.0, 90.0))}, # compatible with selectTools
                'GLN': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)'),
                                 (-180.0, 180.0)),
                        'chi3': (('(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)',
                                  '(segid "%s" and resid %d and name oe1)'),
                                 (-180.0, 180.0))},
                'MET': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name sd)'),
                                 (-180.0, 180.0)),
                        'chi3': (('(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name sd)',
                                  '(segid "%s" and resid %d and name ce)'),
                                 (-180.0, 180.0))},
                'ARG': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)'),
                                 (-180.0, 180.0)),
                        'chi3': (('(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)',
                                  '(segid "%s" and resid %d and name ne)'),
                                 (-180.0, 180.0)),
                        'chi4': (('(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)',
                                  '(segid "%s" and resid %d and name ne)',
                                  '(segid "%s" and resid %d and name cz)'),
                                 (-180.0, 180.0))},
                'LYS': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)'),
                                 (-180.0, 180.0)),
                        'chi3': (('(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)',
                                  '(segid "%s" and resid %d and name ce)'),
                                 (-180.0, 180.0)),
                        'chi4': (('(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)',
                                  '(segid "%s" and resid %d and name ce)',
                                  '(segid "%s" and resid %d and name nz)'),
                                 (-180.0, 180.0))},
                'PRO': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)'),
                                 (-180.0, 180.0)),
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)'),
                                 (-180.0, 180.0)),
                        'chi3': (('(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)',
                                  '(segid "%s" and resid %d and name n)'),
                                 (-180.0, 180.0)),
                        'chi4': (('(segid "%s" and resid %d and name cg)',
                                  '(segid "%s" and resid %d and name cd)',
                                  '(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)'),
                                 (-180.0, 180.0)),  # Added by GAB - Check!
                        'chi5': (('(segid "%s" and resid %d and name cd)',
                                  '(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)'),
                                 (-180.0, 180.0))}, # Added by GAB - Check!
                'CYS': {'chi1': (('(segid "%s" and resid %d and name n)',
                                  '(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name sg)'),
                                 (-180.0, 180.0)),
                        # The following angles apply only for disulfide bridges.
                        'chi2': (('(segid "%s" and resid %d and name ca)',
                                  '(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name sg)',
                                  '(name sg and bondedto (segid "%s" and resid %d and name sg))'),
                                 (-180.0, 180.0)),
                        # I'm not sure whose CYS the angles below belong to.
                        # Also, chi5 (the chi1 of the other CYS) is missing.
                        # (GAB)
                        'chi3': (('(segid "%s" and resid %d and name cb)',
                                  '(segid "%s" and resid %d and name sg)',
                                  '(name sg and bondedto (segid "%s" and resid %d and name sg))',
                                  '(name cb and bondedto (name sg and bondedto (segid "%s" and resid %d and name sg)))'),
                                 (-180.0, 180.0)),
                        'chi4': (('(segid "%s" and resid %d and name sg)',
                                  '(name sg and bondedto (segid "%s" and resid %d and name sg))',
                                  '(name cb and bondedto (name sg and bondedto (segid "%s" and resid %d and name sg)))',
                                  '(name ca and bondedto (name cb and bondedto (name sg and bondedto (segid "%s" and resid %d and name sg))))'),
                                 (-180.0, 180.0))},
                
                # Nucleic acids. 
                'GUA': {'chi': (('(segid "%s" and resid %d and name O4\')',
                                  '(segid "%s" and resid %d and name C1\')',
                                  '(segid "%s" and resid %d and name N9)',
                                  '(segid "%s" and resid %d and name C4)'),
                                 (-180.0, 180.0))},
                'ADE': {'chi': (('(segid "%s" and resid %d and name O4\')',
                                  '(segid "%s" and resid %d and name C1\')',
                                  '(segid "%s" and resid %d and name N9)',
                                  '(segid "%s" and resid %d and name C4)'),
                                 (-180.0, 180.0))},
                'CYT': {'chi': (('(segid "%s" and resid %d and name O4\')',
                                  '(segid "%s" and resid %d and name C1\')',
                                  '(segid "%s" and resid %d and name N1)',
                                  '(segid "%s" and resid %d and name C2)'),
                                 (-180.0, 180.0))},
                'THY': {'chi': (('(segid "%s" and resid %d and name O4\')',
                                  '(segid "%s" and resid %d and name C1\')',
                                  '(segid "%s" and resid %d and name N1)',
                                  '(segid "%s" and resid %d and name C2)'),
                                 (-180.0, 180.0))},
                'URI': {'chi': (('(segid "%s" and resid %d and name O4\')',
                                  '(segid "%s" and resid %d and name C1\')',
                                  '(segid "%s" and resid %d and name N1)',
                                  '(segid "%s" and resid %d and name C2)'),
                                 (-180.0, 180.0))}}

# The following lists have the keys of torsion_info split into protein and
# nucleic acids. This is needed to properly complete torsion_info below.
protein_resnames = ['ALA', 'GLY', 'SER', 'THR', 'VAL', 'ASN', 'ASP', 'HIS',
                    'ILE', 'LEU', 'PHE', 'TRP', 'TYR', 'GLU', 'GLN', 'MET',
                    'ARG', 'LYS', 'PRO', 'CYS']

nucleic_resnames = ['GUA', 'ADE', 'CYT', 'THY', 'URI']

# Complete torsion_info with backbone torsions. (The following works assuming
# all residues in torsion_info belong to protein.)

# Proteins.
for resname in protein_resnames:
    torsion_info[resname]['phi'] = (('(name c and bondedto (segid "%s" and resid %d and name n))',
                                     '(segid "%s" and resid %d and name n)',
                                     '(segid "%s" and resid %d and name ca)',
                                     '(segid "%s" and resid %d and name c)'),
                                    (-180.0, 180.0))
    torsion_info[resname]['psi'] = (('(segid "%s" and resid %d and name n)',
                                     '(segid "%s" and resid %d and name ca)',
                                     '(segid "%s" and resid %d and name c)',
                                     '(name n and bondedto (segid "%s" and resid %d and name c))'),
                                    (-180.0, 180.0))
    torsion_info[resname]['omega'] = (('(segid "%s" and resid %d and name ca)',
                                       '(segid "%s" and resid %d and name c)',
                                       '(name n and bondedto (segid "%s" and resid %d and name c))',
                                       '(name ca and bondedto (name n and bondedto (segid "%s" and resid %d and name c)))'),
                                      (-180.0, 180.0))

# Nucleic acids.
for resname in nucleic_resnames:
    # Standard torsions.
    # ["O3'i-1", "P", "O5'", "C5'"]
    torsion_info[resname]['alpha'] = (('(name O3\' and bondedto (segid "%s" and resid %d and name P))',
                                     '(segid "%s" and resid %d and name P)',
                                     '(segid "%s" and resid %d and name O5\')',
                                     '(segid "%s" and resid %d and name C5\')'),
                                    (-180.0, 180.0))
    # ["P", "O5'", "C5'", "C4'"]
    torsion_info[resname]['beta'] = (('(segid "%s" and resid %d and name P)',
                                     '(segid "%s" and resid %d and name O5\')',
                                     '(segid "%s" and resid %d and name C5\')',
                                     '(segid "%s" and resid %d and name C4\')'),
                                    (-180.0, 180.0))
    # ["O5'", "C5'", "C4'", "C3'"]
    torsion_info[resname]['gamma'] = (('(segid "%s" and resid %d and name O5\')',
                                     '(segid "%s" and resid %d and name C5\')',
                                     '(segid "%s" and resid %d and name C4\')',
                                     '(segid "%s" and resid %d and name C3\')'),
                                    (-180.0, 180.0))
   # ["C5'", "C4'", "C3'", "O3'"]
    torsion_info[resname]['delta'] = (('(segid "%s" and resid %d and name C5\')',
                                     '(segid "%s" and resid %d and name C4\')',
                                     '(segid "%s" and resid %d and name C3\')',
                                     '(segid "%s" and resid %d and name O3\')'),
                                    (-180.0, 180.0))
    # ["C4'", "C3'", "O3'", "Pi+1"]
    torsion_info[resname]['epsilon'] = (('(segid "%s" and resid %d and name C4\')',
                                     '(segid "%s" and resid %d and name C3\')',
                                     '(segid "%s" and resid %d and name O3\')',
              '(name P and bondedto (segid "%s" and resid %d and name O3\'))'),
                                    (-180.0, 180.0))
    # ["C3'", "O3'", "Pi+1", "O5'i+1"]
    torsion_info[resname]['zeta'] = (('(segid "%s" and resid %d and name C3\')',
                                     '(segid "%s" and resid %d and name O3\')',
              '(name P and bondedto (segid "%s" and resid %d and name O3\'))',
                                     '(name O5\' and bondedto (name P and bondedto (segid "%s" and resid %d and name O3\')))'),
                                    (-180.0, 180.0))
    # ["C4'", "O4'", "C1'", "C2'"]
    torsion_info[resname]['nu0'] = (('(segid "%s" and resid %d and name C4\')',
                                     '(segid "%s" and resid %d and name O4\')',
                                     '(segid "%s" and resid %d and name C1\')',
                                     '(segid "%s" and resid %d and name C2\')'),
                                    (-180.0, 180.0))
    # ["O4'", "C1'", "C2'", "C3'"]
    torsion_info[resname]['nu1'] = (('(segid "%s" and resid %d and name O4\')',
                                     '(segid "%s" and resid %d and name C1\')',
                                     '(segid "%s" and resid %d and name C2\')',
                                     '(segid "%s" and resid %d and name C3\')'),
                                    (-180.0, 180.0))
    # ["C1'", "C2'", "C3'", "C4'"]
    torsion_info[resname]['nu2'] = (('(segid "%s" and resid %d and name C1\')',
                                     '(segid "%s" and resid %d and name C2\')',
                                     '(segid "%s" and resid %d and name C3\')',
                                     '(segid "%s" and resid %d and name C4\')'),
                                    (-180.0, 180.0))
    # ["C2'", "C3'", "C4'", "O4'"]
    torsion_info[resname]['nu3'] = (('(segid "%s" and resid %d and name C2\')',
                                     '(segid "%s" and resid %d and name C3\')',
                                     '(segid "%s" and resid %d and name C4\')',
                                     '(segid "%s" and resid %d and name O4\')'),
                                    (-180.0, 180.0))
    # ["C3'", "C4'", "O4'", "C1'"]
    torsion_info[resname]['nu4'] = (('(segid "%s" and resid %d and name C3\')',
                                     '(segid "%s" and resid %d and name C4\')',
                                     '(segid "%s" and resid %d and name O4\')',
                                     '(segid "%s" and resid %d and name C1\')'),
                                    (-180.0, 180.0))

    torsion_info[resname]['chi'] = (('(segid "%s" and resid %d and name O4\')',
                                     '(segid "%s" and resid %d and name C1\')',
                                     '(segid "%s" and resid %d and (((resn ADE or resn GUA) and name N9) or ((resn CYT or resn URI) and name N1)))',
                                     '(segid "%s" and resid %d and (((resn ADE or resn GUA) and name C4) or ((resn CYT or resn URI) and name C2)))'),
                                    (-180.0, 180.0))    
    
    # Pseudotorsions.
    # ["C4'i-1", "P", "C4'", "Pi+1"]
    torsion_info[resname]['eta'] = (('(name C4\' and bondedto (name C3\' and bondedto (name O3\' and bondedto (segid "%s" and resid %d and name P))))',
                                     '(segid "%s" and resid %d and name P)',
                                     '(segid "%s" and resid %d and name C4\')',
                                     '(name P and bondedto (name O3\' and bondedto (name C3\' and bondedto (segid "%s" and resid %d and name C4\'))))'),
                                    (-180.0, 180.0))
    # ["P", "C4'", "Pi+1", "C4'i+1"]
    torsion_info[resname]['theta'] = (('(segid "%s" and resid %d and name P)',
                                     '(segid "%s" and resid %d and name C4\')',
                                     '(name P and bondedto (name O3\' and bondedto (name C3\' and bondedto (segid "%s" and resid %d and name C4\'))))',
                                     '(name C4\' and bondedto (name C5\' and bondedto (name O5\' and bondedto (name P and bondedto (name O3\' and bondedto (name C3\' and bondedto (segid "%s" and resid %d and name C4\')))))))'),
                                    (-180.0, 180.0))
    

# The following functions represent the interface to handle the above defined
# torsion_info data structure.

# Whether torsion_info is kept or modified so that it looks something like its
# more succint equivalent in torsionTools (scTorsionAtoms dictionary), the
# functions below should provide the same interface.

def get_resnames():
    """Return a list of residue names whose torsions are defined in this module."""
    return list(torsion_info.keys())


def get_torsion_names(resname):
    """Return a list of names of torsion angles in residue name resname.

    The argument resname, and the contents of the returned list are strings.
    The order of torsion names in the returned list are subjected to the
    following hierarchy:

    'phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4', 'chi5'
    (proteins)

    'alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi', 'nu0', 'nu1',
    'nu2', 'nu3', 'nu4', 'eta', 'theta'
    (nucleic acids)

    """
    # Note the order list is hardwired: it has to be updated if more torsions
    # are defined in the future.  
    if resname in protein_resnames:
        order = ['phi', 'psi', 'omega', 'chi1', 'chi2', 'chi3', 'chi4', 'chi5']
    elif resname in nucleic_resnames:
        order = ['alpha', 'beta', 'gamma', 'delta', 'epsilon', 'zeta', 'chi',
                 'nu0', 'nu1', 'nu2', 'nu3', 'nu4', 'eta', 'theta']
    resname = resname.strip().upper()
    return sorted(list(torsion_info[resname].keys()), key=order.index)


def get_torsion_selections(resname, torsion, resid=None, segid=''):
    """Return a tuple of 4 atom selections defining torsion in resname.

    The argument resname (string) is a residue name and torsion (string) is the
    name of a torsion angle within the residue (for available names see function
    get_torsion_names).  Alternatively, torsion can be a sequence of torsion
    angle names, in which case a tuple of 4-atom-selection tuples is returned
    (for all torsion angles in the residue simply set torsion to the string
    "all").  For example, torsion = (torsion_i, torsion_j) would return

        ((sel_i1, sel_i2, sel_i3, sel_i4), (sel_j1, sel_j2, sel_j3, sel_j4))
        
    where the elements of torsion tuple are appropriate torsion names, and each
    sel_xy is an atom selection string that selects a single atom.  As hinted
    in the above example, the order in the returned tuple obeys to that of the
    input torsion sequence; when torsion="all" the order is that given by
    get_torsion_names function.
    Arguments resid and segid are optional residue number (integer) and segment
    (string) to incorporate into each XPLOR-style atom selection (string); if
    not provided the corresponding fields in each selection are empty (e.g.,
    'segid "%s" and resid %d and name n').

    """
    resname = resname.strip().upper()
    if type(torsion) == str:
        torsion = torsion.strip().lower()
        if torsion == 'all': torsion = get_torsion_names(resname) # a list

    if type(torsion) == str:  # single torsion requested
        if resid == None:
            return torsion_info[resname][torsion][0]
        else:
            return tuple([x % (segid, resid) for x in  
                          torsion_info[resname][torsion][0]])
    else:  # several torsions requested
        result = []
        for name in torsion:
            name = name.strip().lower()
            if resid == None:
                result.append(torsion_info[resname][name][0])
            else:
                result.append(tuple([x % (segid, resid) for x in  
                                     torsion_info[resname][name][0]]))
        return tuple(result)


def get_torsion_range(resname, torsion):
    """Return a tuple with the min and max possible values (float) of torsion.

    The argument resname (string) is a residue name and torsion (string) is the
    name of a torsion angle within the residue (for available names see function
    get_torsion_names).  Alternatively, torsion can be a sequence of torsion
    angle names, in which case a tuple of (min, max) tuples is returned (for all
    torsion angles in the residue simply set torsion to the string "all").  For
    example, torsion = (torsion_i, torsion_j) would return

        ((min_i, max_i), (min_j, max_j))
        
    where the elements of torsion tuple are appropriate torsion names.  As
    hinted in the above example, the order in the returned tuple obeys to that
    of the input torsion sequence; when torsion="all" the order is that given by
    get_torsion_names function.
    
    """
    resname = resname.strip().upper()
    if type(torsion) == str:
        torsion = torsion.strip().lower()
        if torsion == 'all': torsion = get_torsion_names(resname)  # a list

    if type(torsion) == str:  # single torsion requested
        torsion = torsion.strip().lower()
        return torsion_info[resname][torsion][1]
    else:  # several torsions requested
        result = []
        for name in torsion:
            name = name.strip().lower()
            result.append(torsion_info[resname][name][1])
        return tuple(result)
    
# Function to compare side chains of input model(s) relative to a reference.
# This is also general and could be moved to torsionTools module.

def count_correct_chis(resname, model, reference, accuracy=40.0,
                       excluderesid=[], radialexclude=[], bcutoff=40.0,
                       occupancycutoff=1.0, applyboth=False, wantall=False,
                       verbose=False):
    """
    Get number of correct sidechain chi angles in model relative to a reference.

    'resname' (string) is the residue name (3-letter code) whose side chain
    torsion (chi) angles will be checked against the reference.  Side chains are
    extracted from 'model' (string), a PDB filename, and compared to those of
    'reference' (another PDB file name).  Alternatively, more than one model can
    be specified by making 'model' a sequence of PDB filenames.  For example,
    the models may be NMR structures compared to an X-ray reference.  Two side
    chains are compared only if they share the same residue number and segment
    id (segid) in the model and reference, respectively (WARNING: care must be
    taken when one structure uses the chain id field and the other segid to
    label chains).

    A torsion angle in a model is considered correct if it lies within
    'accuracy' degrees from that in the reference.  A specific residue can be
    explicitly omitted from consideration if its number is included in the
    sequence 'excluderesid' (e.g., we may know a side chain is involved in
    clashes in the reference, and therefore cannot be trusted).  Furthermore, a
    side chain is disregarded if any of its atoms has a B-factor > 'bcutoff' or
    an occupancy < 'occupancycutoff' in the reference structure.
    'radialexclude' is another way of excluding side chains from consideration,
    in this case if they are within angular exclusion zones in the reference
    structure (e.g., badly fit leucines; see docstring of get_residue_torsions
    function for details).  If the B-factor, occupancy, and radial exclusion
    filters need to be additionally applied to the models, 'applyboth' should be
    True.   

    This function returns a tuple.  Assuming a side chain with four chi angles
    (e.g., resname="lys") the returned tuple can be expressed as (assuming
    wantall=True; see below):

    (nchi1, nchi2, nchi3, nchi4, nchi12, nchi123, nchi1234, n)

    where nchix is the number of side chains with accurate chix angles, nchixy
    the number of side chains with both chix and chiy correct, etc., and n is
    the total number of side chains compared.

    If argument wantall is False the function returns:

    (nchi1, nchi12, nchi123, nchi1234, n)

    If no acceptable side chain is found in the reference structure (e.g., all
    have B-factors > bcutoff), an empty tuple is returned.

    If argument verbose is True, detailed information on the torsion-angle
    accuracy of individual side chains is printed.

    The returned numbers can be used to implement a scoring scheme such as that
    described by [Shapovalov et. al (2011) Structure, 19:844].

    """
    # Warning: reference and trial structures may have different
    # chain ids (e.g., reference is xray model and trials are Xplor
    # models where chain id has been removed).

    resname = resname.strip().upper()
    
    if type(model) is str:
        models = [model]  # single input model
    else:
        models = model  # several input models
                
    # List chi names.
    chi_names = get_torsion_names(resname)[3:] # slice out backbone torsions

    # Special case: CYS - don't consider torsions accross disulfide bridge.
    if resname == 'CYS': chi_names = ['chi1']

    nchi = len(chi_names) # number of torsion angles in side chain
    if nchi == 0: raise Exception('no side chain torsion angles in %s' %resname)

    
    protocol.initStruct()  # clear previous structure info, if any
    
    # Load reference structure: structure info is obtained from 'reference', as
    # well as the initial coordinates.
    b_occ = protocol.loadPDB(reference, processSSBonds=False,
                             correctSymmetricSidechains=True,
                             deleteUnknownAtoms=True)

    bfactors = b_occ.bfactors
    occupancies = b_occ.occupancies
    ressel = 'resname ' + resname
    if excluderesid:
        exclude = [str(x) for x in excluderesid]
        ressel = ressel + ' and not (resid ' + ' or resid '.join(exclude) + ')'

    # Get reference side chain torsions.
    refsidechains = get_residue_torsions(torsion_names=chi_names,
                                         radialexclude=radialexclude,
                                         ressel=ressel,
                                         bfactors=bfactors,
                                         occupancies=occupancies,
                                         bfactorcutoff=bcutoff,
                                         occupancycutoff=occupancycutoff)
    if len(refsidechains) == 0: return ()  # no acceptable side chain found

    ncmp = 0  # total number of side chains compared thus far

    # List of type [nchi1, nchi2, nchi3, nchi4, nchi12, nchi123, nchi1234].
    # nchi slots are for the individual angles, and nchi-1 for the "cumulative"
    # counts (i.e., nchi12, etc.).
    counts = (2*nchi-1) * [0]
    
    # Loop through input models.
    for model in models:

        # Get coords of model; it relies on structure info obtained from ref.
        b_occ = protocol.initCoords(files=[model],
                                    correctSymmetricSidechains=True)
        if applyboth:
            bfactors = b_occ.bfactors       # all 0 if model is not X-ray-based
            occupancies = b_occ.occupancies # all 1 if model is not X-ray-based
            # radialexclude is fine as it is.
        else:       # the following scope could be outside loop
            bfactors = None
            occupancies = None
            radialexclude = []  

        # For example, in case of an NMR structure, we don't want to supply
        # 'radialexclude' here (set applyboth=False): if a side chain is wrong
        # we should detect it and reflect it in the counts.
        # However, if the model is another an X-ray structure, we may want to
        # consider only properly fit side chains and supply 'radialexclude' here
        # (set applyboth=False).
        # Whatever the case 'radialexclude' is always applied to the reference.
        
        sidechains = get_residue_torsions(torsion_names=chi_names,
                                          radialexclude=radialexclude,
                                          ressel=ressel, bfactors=bfactors,
                                          occupancies=occupancies,
                                          bfactorcutoff=bcutoff,
                                          occupancycutoff=occupancycutoff)
        # Header for verbose printout.
        msg='structure , x-resname-x resid (segid) , [nchi1,nchi2,...,nchi1234]'
        if verbose: print(msg)
        
        # Match side chains in model and reference.
        for ref_sc in refsidechains:
            for sc in sidechains:
                
                if sc == ref_sc: # effectively, True if segids and resids match

                    ncmp += 1
    
                    # Determine accuracy of individual torsions.
                    chi_counts = nchi * [0] # 0 <=> chi outside 'accuracy' range
                    info = []  # for verbose
                    for chi in chi_names:
                        x = sc.torsions[chi].value
                        ref = ref_sc.torsions[chi].value
                        upper = max(sc.torsions[chi].range)
                        lower = min(sc.torsions[chi].range)
                        domain = upper - lower
                        diff = get_smallangle(x, ref, domain)
                        info.append((chi, x, ref, diff))  # for verbose
                        if diff < accuracy:
                            chi_counts[chi_names.index(chi)] = 1 # flag as
                                                                 # accurate                                   
                    # Determine "cumulative" counts.
                    cum = (nchi-1) * [0]
                    for i in range(nchi-1):
                        if chi_counts[i] + chi_counts[i+1] == 2:
                            cum[i] = 1
                        else:
                            break  # e.g., if either chi1 or chi2 is outside
                                   # 'accuracy', then not only chi12 is, but
                                   # also chi123 and chi1234

                    # Both 'chi_counts' and 'cum' have binary items (0 or 1).

                    # Heart of verbose printout.
                    if verbose:
                        print(model, ',', sc, ',', chi_counts+cum)
                        for x in info:
                            print(x[0], '(model, ref, diff):', x[1], x[2], x[3])

                    # Update counts with current side chain.
                    counts = [x+y for (x,y) in zip(counts, chi_counts+cum)]

    # Done counting angles. 

    # If required, keep only nchi1 and cumulative counts (if they exist).
    if wantall is False:
        if nchi > 1: counts = [counts[0]] + counts[-(nchi-1):]
        # else 'counts' already has all the info we want: nchi1.

    counts.append(ncmp) # total number of comparisons done
    
    return tuple(counts)    # order is important, hence the tuple

 
#______________________________________________________________________________              



class ResidueTorsionConfig:
    """A residue's point in torsion-angle configuration (sub)space.

    The torsion-angle space associated with a particular residue is given by
    all the torsion-angle dimensions that define the residue configuration,
    namely omega, phi, psi, chi1, etc.  A residue can be specified by a point
    in such space or any of the possible subspaces (e.g., Ramachandran
    subspace).

    """
    def __init__(self, pdbid=None, segid=None, resid=None, resname=None,
                 prevresname=None, postresname=None, torsions={}):
        self.pdbid = pdbid
        self.segid = segid
        self.resid = resid
        self.resname = resname
        self.prevresname = prevresname
        self.postresname = postresname
        self.torsions = torsions  # dict with Torsion.name: Torsion
    def __eq__(self, other):
        return (self.pdbid == other.pdbid and self.segid == other.segid and
                self.resid == other.resid)
    def __repr__(self):
        return '%s-%s-%s %i (%s%s)' % (self.prevresname, self.resname,
                                       self.postresname, self.resid,
                                       self.pdbid, self.segid)
        

class Torsion:
    """Represents a torsion angle."""
    def __init__(self, name=None, value=None, anglerange=[], maxB=None,
                 minoccupancy=None):
        self.name = name  # phi, psi, ...
        self.value = value
        self.range = anglerange # sequence with min and max angle values
        self.maxB = maxB  # largest B-factor of the 4 atoms defining the angle
        self.minoccupancy = minoccupancy  # min occupancy of the 4-atom set



def accepted_torsion(torsion_name, segid, resid, resname, sim, bfactors,
                     occupancies, bfactorcutoff, occupancycutoff):
    """Return a Torsion instance if B-factor and occupancy criteria are met.

    Raises an exception if any atom that defines the torsion_name torsion angle
    in input residue (specified via its segid, resid, and resname) is missing,
    has a B-factor larger than bfactorcutoff, or an occupancy less than
    occupancycutoff.  bfactors and occupancies are the corresponding attributes
    of an InitCoordsResult instance.  sim is a simulation.Simulation instance. 

    """
    # The following will raise a KeyError exeption when asking for an
    # inexistent torsion angle in resname (desired behavior).
    seltemplates = get_torsion_selections(resname, torsion_name)

    sel0 = atomSel.AtomSel(seltemplates[0] % (segid, resid), sim)
    sel1 = atomSel.AtomSel(seltemplates[1] % (segid, resid), sim)
    sel2 = atomSel.AtomSel(seltemplates[2] % (segid, resid), sim)
    sel3 = atomSel.AtomSel(seltemplates[3] % (segid, resid), sim)

    # Check for missing atoms.
    if len(sel0) == 0 or len(sel1) == 0 or len(sel2) == 0 or len(sel3) == 0:
#        print 'Missing atoms'  # test
        raise Exception('Missing atoms')

    atom0 = selectTools.convertToAtom(sel0)
    atom1 = selectTools.convertToAtom(sel1)
    atom2 = selectTools.convertToAtom(sel2)
    atom3 = selectTools.convertToAtom(sel3)

    # Check all atoms defining the torsion-angle satisfy B-factor criterion.
    if (bfactors[atom0.index()] > bfactorcutoff or
        bfactors[atom1.index()] > bfactorcutoff or
        bfactors[atom2.index()] > bfactorcutoff or
        bfactors[atom3.index()] > bfactorcutoff):
#        print 'B-factor too high'  # test
        raise Exception('B-factor too high')

    # The largest B-factor (lower that bfactorcutoff) in the atom set.
    maxB = max(bfactors[atom0.index()], bfactors[atom1.index()],
               bfactors[atom2.index()], bfactors[atom3.index()])

    # Check all atoms defining the torsion-angle satisfy occupancy criterion.
    if (occupancies[atom0.index()] < occupancycutoff or
        occupancies[atom1.index()] < occupancycutoff or
        occupancies[atom2.index()] < occupancycutoff or
        occupancies[atom3.index()] < occupancycutoff):
#        print 'Occupancy too low'  # test
        raise Exception('Occupancy too low')

    # The smallest occupancy (larger than occupancycutoff) in the atoms set.
    minoccupancy = min(occupancies[atom0.index()], occupancies[atom1.index()],
                       occupancies[atom2.index()], occupancies[atom3.index()])

    torsion_value_rad = dihedral.Dihedral(sel0, sel1, sel2, sel3) # angle in rad
    torsion_value = torsion_value_rad.value() * 180.0/math.pi # angle in degrees

    torsion_range = get_torsion_range(resname, torsion_name)
    
    torsion = Torsion(name=torsion_name, value=torsion_value,
                      anglerange=torsion_range, maxB=maxB,
                      minoccupancy=minoccupancy)
    return torsion


def neighborresname(segid, resid, targetresnames_selstring, shift, sim):
    """Return name of residue neighboring (segid, resid) if in targetresnames_selstring

    The residue number of the neighboring residue is resid+shift.
    targetresnames_selstring is an XPLOR resname selection string of the type
    '(resname ALA or resname VAL)'.  The name of the neighbor residue is
    returned if it is in targetresnames_selstring or targetresnames_selstring=
    '(resname *)'.  An exception of the latter case is when the neighboring
    residue number doesn't have any atoms (e.g, falls outside the sequence),
    where '   ' is returned.  Otherwise, returns None.  (Note that if
    targetresnames_selstring is different from '(resname *)' and resid+shift
    falls outside the sequence, it returns None.)

    """
    temp = 'segid "%s" and resid %d and ' % (segid, resid+shift)
    selstring = temp + targetresnames_selstring    
    neighborressel = atomSel.AtomSel(selstring, sim)
    if neighborressel:
        neighborresname = neighborressel[0].residueName()
    else:
        if targetresnames_selstring == '(resname *)': # all neighbor types requ-
            neighborresname = '   ' # ested but no one found (e.g., Nterm-resid)
        else:
            neighborresname = None
    return neighborresname


def resname_selstring(resnameseq):
    """Return an XPLOR resname selection string.

    Returns '(resname a or resname b or resname c)', given the input sequence
    of strings 'a', 'b', and 'c'.

    """
    seltemplate = 'resname %s or '
    selstring = ''
    for resname in resnameseq:
        selstring += seltemplate % resname
    selstring = '(' + selstring[:-4] + ')'
    return selstring


def complement_resnames(resnames):
    """Return a list of all residue names (strings), except for those in input.

    resnames is either a residue name or a sequence of residue names, each
    being a 3-letter code residue name string.

    """
    if type(resnames) == str: resnames = [resnames]
    resnames = set([x.strip().upper() for x in resnames])
    all_resnames = set(get_resnames())
    return list(all_resnames.difference(resnames))
    
    

def get_smallangle(angle1, angle2, domain):
    """Return the shortest distance between two vectors in angle space.

    The vectors are specified by their multidimensional angular coordinates
    given by sequences angle1 and angle2.  The domain in each dimension is
    given by the domain sequence.  Alternatively, in a one-dimensional case,
    each argument may be a single number.
    
    """
    def smallangle(angle1, angle2, domain):
        """Return small angle in one dimension."""
        diff = abs(angle1 - angle2)
        if diff < domain/2.0:
            return diff
        else:
            return domain - diff
        
    if type(angle1) in (float, int):  # assume other arguments in same format
        angle1 = [angle1]
        angle2 = [angle2]
        domain = [domain]

    data = list(zip(angle1, angle2, domain))
    return sum([smallangle(x[0],x[1],x[2])**2.0 for x in data])**0.5


# Make sure argument default values match those in function
# get_residue_torsions_from_files() below.
def get_residue_torsions(pdbid='', torsion_names=[], radialexclude=[],
                         radialinclude=[], ressel='tag', prevresnames=['*'],
                         postresnames=['*'], bfactors=None, occupancies=None,
                         bfactorcutoff=30.0, occupancycutoff=1.0,
                         sim=simulation.currentSimulation()):
    """Return a list of ResidueTorsionConfig instances from current coordinates.

    torsion_names is a sequence of one or more torsion angle names (for
    available names see function get_torsion_names; Exception is raised if
    torsion_names is not provided) queried from residue selection ressel (in
    XPLOR format) in the protein optionally labeled by its pdbid.  Optional
    lists of previous and posterior residue names (prevresnames and
    postresnames, respectively) can be specified; for example, ressel=
    'segid A and resname ALA' with postresnames=['PRO'] selects all pre-Proline
    Alanines in segment A.  If not specified, all previous and posterior residue
    names are considered.  bfactors and occupancies are the corresponding
    attributes of the <m protocol>.InitCoordsResult instance generated upon
    coordinate initialization with <m protocol>.initCoords; if not provided all
    atoms are given a bfactor of 0.0 and an occupancy of 1.0.  Torsion angles
    involving atoms with B-factors larger than bfactorcutoff and/or occupancies
    smaller than occupancycutoff are not output.  Similarly, torsion angles are
    not output if within a user-specified radial exclusion zone.  One or more
    exclusion zones are specified via the optional radialexclude argument, a
    sequence containing one or more sequences with the center coordinates and
    the radius of each zone.  For example, two exclusion zones in a
    two-dimensional case would be specified as:
    
        radialexclude = [((x1i, x2i), radiusi), ((x1j, x2j), radiusj)]

    Similarly, an optional radial inclusion zone can be specified via the
    argument radialinclude (with format identical to radialexclude); if not
    specified the entire space is assumed wanted.
    
    IMPORTANT: before running this function remove atoms with arbitrary
    coordinates (i.e., "unknown atoms" in XPLOR jargon) from the input
    structure, e.g., by setting argument deleteUnknownAtoms=True in
    <m protocol>.initCoords() function.

    """
    # If we don't do the following, an emtpy list will pass the else clause
    # in the loop over torsion_names because it won't hit the break!!!
    if torsion_names == []: raise Exception('no torsion_names provided')
    
    # Set default B-factors and occupancies.
    if bfactors == None: bfactors = [0.0] * sim.numAtoms()
    if occupancies == None: occupancies = [1.0] * sim.numAtoms()

    # Dictionary with segid: [(resid, resname), ...] elements.
    segments = selectTools.getSegsResidues(ressel, sim)

    # List to contain ResidueTorsionConfig instances returned by this function.
    result = []

    prevresnames_selstring = resname_selstring(prevresnames)
    postresnames_selstring = resname_selstring(postresnames)    

    for segid in list(segments.keys()):
        for (resid, resname) in segments[segid]:
                                                          
            prevresname = neighborresname(segid, resid, prevresnames_selstring,
                                          -1, sim)
            postresname = neighborresname(segid, resid, postresnames_selstring,
                                          +1, sim)
            
            # Get torsion angles only if the previous and posterior
            # residues are among the queried types.
            if prevresname and postresname:
                torsions = {}  # to hold queried torsion angles as Torsion inst.
                torsion_values = []  # to hold queried torsion-angle values
                torsion_domains = [] # to hold the domains of queried torsions
                for torsion_name in torsion_names:
                    try:
                        torsion = accepted_torsion(torsion_name, segid, resid,
                                                   resname, sim, bfactors,
                                                   occupancies, bfactorcutoff,
                                                   occupancycutoff)
                        torsions[torsion.name] = torsion
                        torsion_values.append(torsion.value)
                        
                        bounds = get_torsion_range(resname, torsion_name)
                        torsion_domains.append(max(bounds)-min(bounds))
                    except Exception:
                        break  # discard if not all queried angles can be got
                else:  # only if torsions list is complete:

                    accept = True

                    for (center, radius) in radialexclude:
                        if get_smallangle(center, torsion_values,
                                          torsion_domains) < radius:
                            accept = False

                    for (center, radius) in radialinclude:
                        if get_smallangle(center, torsion_values,
                                          torsion_domains) > radius:
                            accept = False
                        
                    # Accept only if outside exclusion and inside inclusion zone
                    if accept: 
                        
                        # By default, segid is assigned the one-character
                        # chainID field in the pdb file, overriding the (longer)
                        # segID pdb field.  However, when there's only one chain
                        # and chainID='', segid is assigned the segid pdb field,
                        # which may hold the pdbID; we handle such case so that
                        # the ResidueTorsionConfig instance has an appropriately
                        # labeled segment.  (Note that simply clearing the segID
                        # field in the pdb file itself before its input works
                        # too.)
                        if len(segid) > 1:
                            segidlabel = ''     # do not change segid variable
                        else:
                            segidlabel = segid  # do not change segid variable

                        residueconfig = ResidueTorsionConfig(pdbid=pdbid,
                                                             segid=segidlabel,
                                                             resid=resid,
                                                             resname=resname,
                                                        prevresname=prevresname,
                                                        postresname=postresname,
                                                             torsions=torsions)
    #                    print residueconfig  # test
                        result.append(residueconfig)
    return result


# Make sure argument default values match those in get_residue_torsions.
def get_residue_torsions_from_files(pdbfiles='', pdbfiledir='.',
                                    torsion_names=None, radialexclude=[],
                                    radialinclude=[], ressel='tag',
                                    prevresnames=['*'], postresnames=['*'],
                                    bfactorcutoff=30.0, occupancycutoff=1.0):
    """Return a list of ResidueTorsionConfig instances from input pdb files.

    pdbfiles is a string with pdb file names, one per line, from which to get
    the returned ResidueTorsionConfig instances.  pdbfiledir is a string with
    the directory path where the files are found.  The remaining arguments are
    those needed to run function get_residue_torsions, at the heart of this
    function (see get_residue_torsions for more information).

    """
    result = []
    for pdbfile in pdbfiles.splitlines():
        if pdbfile:
            pdbid = pdbfile.split('.pdb')[0]  # all before extention '.pdb'
            filename = pdbfiledir + '/' + pdbfile

            protocol.initStruct()  # clear structure info
            
            B_occ = protocol.loadPDB(file=filename, processSSBonds=False,
                                     correctSymmetricSidechains=True,
                                     deleteUnknownAtoms=True)

            residueconfigs = get_residue_torsions(pdbid=pdbid,
                                                  torsion_names=torsion_names,
                                                  radialexclude=radialexclude,
                                                  radialinclude=radialinclude,
                                                  ressel=ressel,
                                                  prevresnames=prevresnames,
                                                  postresnames=postresnames,
                                                  bfactors=B_occ.bfactors,
                                                  occupancies=B_occ.occupancies,
                                                  bfactorcutoff=bfactorcutoff,
                                                occupancycutoff=occupancycutoff)
            result.extend(residueconfigs)
    return result


def separate_torsion_outliers(residues, radius, percent):
    """Return a tuple with (list of non-outlier residues, list of outliers).

    residues is a sequence of ResidueTorsionConfig instances.  For each
    d-dimensional torsion angle in the sequence, the number of torsion angles
    within the specified radius is computed; if such number represents less than
    the specified percent of the data, the torsion angle is considered an
    "outlier".  The returned tuple contains as first and second elements a list
    of non-outlier and outlier ResidueTorsionConfig instances, respectively.

    """
    n = len(residues)
    outliers = []
    inliers = []
    for i in range(n):
        count = 0
        items = list(residues[i].torsions.items())
        items.sort()  # by torsion name
        torsion1 = [x[1].value for x in items]
        domain = [max(x[1].range)-min(x[1].range) for x in items]
        for j in range(n):
            if i != j:
                items = list(residues[j].torsions.items())
                items.sort()  # by torsion name
                torsion2 = [x[1].value for x in items]                
                distance = get_smallangle(torsion1, torsion2, domain)
                if distance < radius: count += 1
        if 100.0*count/n < percent:
            outliers.append(residues[i])
        else:
            inliers.append(residues[i])
    return (inliers, outliers)
        
    
                               




def get_binvalue(binnumber, binwidth, bounds):
    """Return bin value given bin number, bin width and overall value bounds."""
    binvalue = min(bounds) + (binnumber + 0.5) * binwidth
    return binvalue


def set_table(ranges, binwidths, usebinnumbers=False):
    """Return a dictionary with coordinates-0 key-value pairs.

    The returned dictionary represents an n-dimensional contingency table,
    where the key is a tuple of numbers that specifies the bin coordinates and
    the value is an integer with the bin count (count=0 in this case, i.e.,
    empty contingency table).  binwidths is a sequence of the required bin width
    in each dimension.  ranges is a sequence with the min and max values that
    define the domain of each dimension.  The dimensionality is given by
    len(binwidths) (=len(ranges)).  Bin widths and ranges are matched by
    positional offset: the order is important!  Output bin coordinates can be
    expressed in terms of bin numbers (if usebinnumbers=True) or the value at
    center of bin (default).

    """
    (d, b) = (len(binwidths), len(ranges))
    if d != b:
        raise Exception('binwidths has dimension %s, and ranges %s' % (d, b))
    
    bin_counts = ()  # to hold the total number of bins in each dimension
    for (binwidth, bounds) in zip(binwidths, ranges):
        bin_count = (max(bounds) - min(bounds))/float(binwidth)
        if bin_count != int(bin_count):
            raise Exception('Non-integral number of bins')
        bin_counts += (int(bin_count),)
    
##    test = {}             # test 1 (uncomment "test 1" lines for testing)
    table = {}
    table_points = tuple([list(range(x)) for x in bin_counts])
    for binnumbers in itertools.product(*table_points):
        if usebinnumbers:  
            table[binnumbers] = 0  # zero count (returned table is empty)
        else:               
            binvalues = ()
            for (binnumber, binwidth, bounds) in zip(binnumbers, binwidths,
                                                     ranges):
                binvalue = get_binvalue(binnumber, binwidth, bounds)
                binvalues += (binvalue,)
            table[binvalues] = 0  # zero count (returned table is empty)
##            test[binnumbers] = binvalues  # test 1
##    return test                           # test 1
    return table


def bin_data(data, ranges, binwidths, pdf=False, logpdf=False,
             usebinnumbers=False):
    """Return a dictionary with coordinates-bin_value key-value pairs.

    If pdf=False and logpdf=False, the returned dictionary represents an
    n-dimensional contingency table, where the key is a tuple of numbers that
    specifies the bin coordinates and the value is an integer with the bin
    count.  If pdf=True, the returned dictionary is a probability density
    function (pdf) estimate, with bin counts divided by the total counts and the
    bin volume, giving the densities.  If logpdf=True, the negative of the
    natural log of the densities is calculated; an arbitrary small value is
    assigned for zero densities.  An exception is raised if both pdf=True and
    logpdf=True.
    data is a sequence of n-dimensional data points to be binned.  Each point is
    a sequence of values representing the coordinates of the point; the length
    of the sequence determines the dimensionality of the contingency table.
    binwidths is a sequence of the required bin width in each dimension.
    ranges is a sequence with the min and max values that define the domain of
    each dimension.  Data points, bin widths, and ranges are matched by
    positional offset: the order is important!  Output bin coordinates can be
    expressed in terms of bin numbers (if usebinnumbers=True) or the value at
    center of bin (default).

    """
    # Checks.
    (d, b) = (len(data[0]), len(binwidths))
    if d != b:
        raise Exception('data has dimension %s, and binwidths %s' % (d, b))

    # (Dimensional compatibility among data points not checked.)

    if pdf and logpdf:
        raise Exception('Choose either pdf or logpdf, not both')
    
    # Binning.
    table = set_table(ranges, binwidths, usebinnumbers)  # empty table
    for slot in list(table.keys()):  # slot = bin
        for datum in data:
            dimensions = list(zip(slot, binwidths, ranges, datum))
            for (coordinate, binwidth, bounds, value) in dimensions:
                if usebinnumbers:
                    binvalue = get_binvalue(coordinate, binwidth, bounds)
                else:
                    binvalue = coordinate
                if not ((binvalue - binwidth/2.0) < value and
                        value <= (binvalue + binwidth/2.0)):
                    break  # if not in this dimension of bin, don't check others
            else:
                table[slot] += 1  # increment bin count

    # At this point we have a contingency table. If either a pdf or logpdf is
    # required, modify it accordingly before returning it.
    if pdf or logpdf:
        total_count = sum(table.values())
        bin_volume = 1.0
        for binwidth in binwidths:
            bin_volume *= binwidth
        for slot in list(table.keys()):
            count = table[slot]
            #FIX: check integer division - CDS 2020/02/12
            if pdf: table[slot] = count/(total_count * bin_volume)
            if logpdf:
                if count == 0: count = 0.5  # arbitrary; to avoid log error
                #FIX: check integer division - CDS 2020/02/12
                table[slot] = -math.log(count/(total_count * bin_volume))

    return table  # contingency table, pdf, or -log(pdf)
    
    
def bin_torsions(configs, resolutions={}, binwidth=5.0, pdf=False, logpdf=False,
                 usebinnumbers=False):
    """Return a dictionary with coordinates-bin_value key-value pairs.

    If pdf=False and logpdf=False, the returned dictionary represents an
    n-dimensional contingency table, where the key is a tuple of numbers that
    specifies the bin coordinates and the value is an integer with the bin
    count.  If pdf=True, the returned dictionary is a probability density
    function (pdf) estimate, with bin counts divided by the total counts and the
    bin volume, giving the densities.  If logpdf=True, the negative of the
    natural log of the densities is calculated; an arbitrary small value is
    assigned for zero densities.  An exception is raised if both pdf=True and
    logpdf=True.
    configs is a sequence of ResidueTorsionConfig instances whose associated
    torsion angle values are to be binned.  The size of the torsions attribute
    determines the dimensionality of the contigency table.  The coordinates of
    each bin are ordered by the alphabetical order of the associated torsion
    angle names: 'chi1', 'chi2', 'chi3', 'chi4', 'phi', 'psi', 'omega'.
    resolutions is a dictionary with the bin widths for each torsion dimension,
    e.g., {'phi': 5, 'psi': 10}.  For any torsion angle not specified in
    resolutions, the bin width is set to binwidth (degrees).  Output bin
    coordinates can be expressed in terms of bin numbers (if usebinnumbers=True)
    or the value at center of bin (default).

    """
    # (Check for size of data?)

    # Set order of torsion-angle dimensions according to alphabetical
    # order of the torsion names.
    order = list(configs[0].torsions.keys())
    order.sort()
    order = tuple(order)

    # Format parameters for input to bin_data().
    ranges = ()
    binwidths = ()
    for torsion_name in order:
        ranges += (configs[0].torsions[torsion_name].range,)
        binwidths += (resolutions.get(torsion_name, binwidth),)

    data = []
    for config in configs:
        datum = ()
        for torsion_name in order:
            datum += (config.torsions[torsion_name].value,)
        data.append(datum)

    # Binning and return.
    return bin_data(data, ranges, binwidths, pdf, logpdf, usebinnumbers)



    
 





    

    
            
    

    




    
