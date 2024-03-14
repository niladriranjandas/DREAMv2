
"""tools to aid in structure analysis. As of this writing, this consists of
torsion angle analysis.

This is a preliminary version and not well tested.
"""

class torsionAngles:
    def nucleicAcid(sel):
        """ returns a dictionary of torsion angle values keyed by the torsion
        angle name for each torsion angle in
        the first residue referred to in the <m atomSel>.AtomSel argument.
        """

        torsionList = [("alpha",  ("Prev-O3'", "P", "O5'", "C5'")),
                       ("beta",   ("P", "O5'", "C5'", "C4'")),
                       ("gamma",  ("O5'", "C5'", "C4'", "C3'")),
                       ("delta",  ("C5'", "C4'", "C3'", "O3'")),
                       ("epsilon",("C4'", "C3'", "O3'", "Next-P")),
                       ("zeta",   ("C3'", "O3'", "Next-P", "Next-O5'"))]
        resname = sel[0].residueName()

        if (resname=='URA' or resname=='URI' or
            resname=='CYT' or resname=='THY'  ):
            torsionList.append(('chi',("O4'", "C1'", "N1", "C2")))
            pass
        if (resname=='ADE' or resname=='GUA'):
            torsionList.append(('chi',("O4'", "C1'", "N9", "C4")))
            pass
        
        return getAngles(sel,torsionList)
    def protein(sel):
        """ returns a dictionary of torsion angle values keyed by the torsion
        angle name for each torsion angle in
        the first residue referred to in the <m atomSel>.AtomSel argument.
        """
        torsionList = [("phi", ("Prev-C", "N", "CA", "C")),
                       ("psi", ("N", "CA", "C", "Next-N")),
                       ("omega",("CA", "C", "Next-N", "Next-CA"))]

        resname = sel[0].residueName()

        try:
            torsionList += scTorsionAtoms[resname]
        except KeyError:
            pass

        return getAngles(sel,torsionList)

    nucleicAcid = staticmethod(nucleicAcid)
    protein     = staticmethod(protein)
    pass

scTorsionAtoms = {}
scTorsionAtoms['GLY'] = []
scTorsionAtoms['ALA'] = []
scTorsionAtoms['SER'] = [('CHI1',('N', 'CA', 'CB', 'OG')),
                         ('CHI2',('CA','CB', 'OG', 'HG'))]
scTorsionAtoms['THR'] = [('CHI1' , ('N' , 'CA', 'CB', 'OG1')),
                         ('CHI21', ('CA', 'CB', 'OG1', 'HG1'))]
scTorsionAtoms['LYS'] = [('CHI1', ('N',  'CA', 'CB', 'CG' )),
                         ('CHI2', ('CA', 'CB', 'CG', 'CD' )),
                         ('CHI3', ('CB', 'CG', 'CD', 'CE' )),
                         ('CHI4', ('CG', 'CD', 'CE', 'NZ' ))]
scTorsionAtoms['CYS'] =[('CHI1', ('N',  'CA', 'CB', 'SG')),
                        ('CHI2', ('CA', 'CB', 'SG', 'HG'))]
scTorsionAtoms['MET'] = [('CHI1', ('N',  'CA', 'CB', 'CG' )),
                         ('CHI2', ('CA', 'CB', 'CG', 'SD' )),
                         ('CHI3', ('CB', 'CG', 'SD', 'CE' ))]
scTorsionAtoms['VAL'] =[('CHI1',  ('N',  'CA', 'CB',  'CG1' ))]
scTorsionAtoms['ILE'] = [('CHI1',  ('N',  'CA',  'CB',  'CG1')),
                         ('CHI21', ('CA', 'CB',  'CG1', 'CD1'))]
scTorsionAtoms['LEU'] = [('CHI1',  ('N',  'CA',  'CB',  'CG'   )),
                         ('CHI2',  ('CA', 'CB',  'CG',  'CD1'  ))]
scTorsionAtoms['ASP'] = [('CHI1',  ('N',  'CA', 'CB',  'CG'  )),
                         ('CHI2',  ('CA', 'CB', 'CG',  'OD1' ))]
scTorsionAtoms['ASN'] = [('CHI1',  ('N',  'CA', 'CB',  'CG'  )),
                         ('CHI2',  ('CA', 'CB', 'CG',  'OD1' ))]
scTorsionAtoms['GLU'] = [('CHI1',  ('N',  'CA', 'CB',  'CG'  )),
                         ('CHI2',  ('CA', 'CB', 'CG',  'CD'  )),
                         ('CHI3',  ('CB', 'CG', 'CD',  'OE1' ))]
scTorsionAtoms['GLN'] = [('CHI1',  ('N',  'CA', 'CB',  'CG'  )),
                         ('CHI2',  ('CA', 'CB', 'CG',  'CD'  )),
                         ('CHI3',  ('CB', 'CG', 'CD',  'OE1' ))]
scTorsionAtoms['ARG'] = [('CHI1',  ('N',  'CA', 'CB',  'CG'  )),
                         ('CHI2',  ('CA', 'CB', 'CG',  'CD'  )),
                         ('CHI3',  ('CB', 'CG', 'CD',  'NE'  )),
                         ('CHI4',  ('CG', 'CD', 'NE',  'CZ'  ))]
scTorsionAtoms['PRO'] = []
scTorsionAtoms['HIS'] = [('CHI1',  ('N',  'CA', 'CB', 'CG'  )),
                         ('CHI2',  ('CA', 'CB', 'CG', 'ND1' ))]
scTorsionAtoms['PHE'] = [('CHI1',  ('N',  'CA', 'CB', 'CG' )),
                         ('CHI2',  ('CA', 'CB', 'CG', 'CD1'))]
scTorsionAtoms['TYR'] = [('CHI1',   ('N',   'CA', 'CB', 'CG' )),
                         ('CHI2',   ('CA',  'CB', 'CG', 'CD1')),
                         ('CHI6',   ('CE1', 'CZ', 'OH', 'HH' ))]
scTorsionAtoms['TRP'] = [('CHI1',  ('N',  'CA', 'CB', 'CG' )),
                         ('CHI2',  ('CA', 'CB', 'CG', 'CD1'))]

            
    
def getAngles(sel,torsionList):
    """ returns a dictionary of torsion angle values (in Radians) keyed by
    the torsion angle name for each torsion angle in
    the first residue referred to in the <m atomSel>.AtomSel argument.
    """

    from math import pi
    from atomSel import AtomSel
    from dihedral import Dihedral

    resid = sel[0].residueNum()
    sim = sel.simulation()
    segid = sel[0].segmentName()

    ret = {}
    for (name,atoms) in torsionList:
        atoms=[]
        for atom in atoms:
            resNum=resid
            if atom.startswith("Prev-"):
                resNum -= 1
                atom = atom[5:]
                pass
            if atom.startswith("Next-"):
                resNum += 1
                atom = atom[5:]
                pass
            selStr="name %s and resid %d" % (atom,resNum)
            if segid:
                selStr += " and segid %4s" % segid
                pass
            atoms.append( AtomSel(selStr,sim) )
            pass
                
        try:
            d= Dihedral(*atoms)
            # convert to radians
            ret[name] = d.value() * pi/180
        except:
            pass
        pass
    return ret
