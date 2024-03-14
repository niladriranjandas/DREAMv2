"""Tools for surface area calculation
"""
from surfaceArea import SurfaceArea

def create_SurfaceArea(rSolvent=1.2,
                       sel="known",
                       radiusType="allatom"
                       ):
    """create a <m surfaceArea>.SurfaceArea object with the specified
    solvent radius. Surface area calculations will include atoms specified
    by the <m atomSel>.AtomSel sel.

    if radiusType is allatom, radii are set to be rSolvent + 0.5 * sigma of
    the NONBond parameter - of all XPLOR type-based parameters read in thus
    far.

    if radiusType is heavyatom, the radii from McCammon, Wolynes and Karplus,
    Biochem. 18, 927-942 (1979) are used, and protons positions are not
    considered.
    
    """
    radii = []
    from atomSel import AtomSel
    if type(sel) == type('string'):
        sel = AtomSel(sel)
        pass

    if radiusType=='allatom':
        global chemTypeRadiusLookup
        for atom in AtomSel("all"):
            try:
                radius = chemTypeRadiusLookup[ atom.chemType() ]
            except KeyError:
                initializeChemTypeRadii()
                radius = chemTypeRadiusLookup[ atom.chemType() ]
                pass
            radii.append(radius + rSolvent)
            pass
        pass
    elif radiusType=='heavyatom':
        selection=sel
        heavyAtoms = AtomSel("(%s) and (not name H*)" % selection.string())

        # vdw radii of heavy atoms with protons included
        # taken from McCammon, Wolynes and Karplus, Biochem. 18, 927-942 (1979)
        # for each element, the entries are for 0,1,... attached protons
        chemRadii={}
        chemRadii['O'] = [1.60, 1.70]
        chemRadii['N'] = [1.65,1.65,1.70,1.75]
        chemRadii['C'] = [1.80,1.85,1.90,1.95]
        chemRadii['S'] = [1.90,1.90]
        chemRadii['P'] = [1.95,1.95] #Guess CDS 2016/07/06

        radii = [0]*selection.simulation().numAtoms()
        for atom in heavyAtoms:
            numH = len( AtomSel("bondedto id %d and name H*" %
                                (atom.index()+1)) )
            elem = atom.atomName()[0]
            index = atom.index()

            if elem=='C' and numH==1:
                # check for aromatic CH group
                if len( AtomSel("bondedto id %d" % (atom.index()+1)) ) == 3:
                    numH+=1
                    pass
                pass
            if numH>= len(chemRadii[elem]):
                print('create_DensPot: WARNING: unsupported number of protons', end=' ')
                print('(%d) for %s' % (numH, elem))
                print('\tatom:', atom.string())
                numH -= 1
                pass
            
            radii[index] = chemRadii[elem][numH]+rSolvent
        
            pass
        sel = heavyAtoms
        pass
    else:
        raise Exception("invalid radiusType: " + radiusType)
            
    sa = SurfaceArea(radii,sel)
#    sa.setRSolvent(rSolvent)
    return sa

def initializeChemTypeRadii():
    from tempfile import mktemp #FIX: possible security problem
    import os
    from xplorSimulation import getXplorSimulation
    xSim = getXplorSimulation()
    tmpFilename=mktemp(suffix=".xplor")
    outputState=xSim.disableOutput()
    xSim.command("write params output=%s end" % tmpFilename)
    xSim.enableOutput(outputState)
    readXplorRadii(tmpFilename)
    os.unlink(tmpFilename)
    return

chemTypeRadiusLookup = {}
def readXplorRadii(filename):
    import re
    global chemTypeRadiusLookup
    for line in open(filename).readlines():
        #print line
        match = re.search(r"^\s*nonb.*\s+([a-z0-9]+)\s+([0-9.]+)\s+([0-9.]+)",
                          line,re.IGNORECASE)
        if match:
            name = match.group(1)
            sigma = float( match.group(3) )
            #print "found %s %f" % (name, sigma)
            chemTypeRadiusLookup[name] = 0.5 * sigma
            pass
        pass
    return

def calcSA(rSolvent=1.2,
           sel="known"):
    """calculate the total effective surface area of the given selection.
    Radii are calculated as documented for create_SurfaceArea.
    """
    from atomSel import AtomSel
    sa = create_SurfaceArea(rSolvent,sel)
    if type(sel)==type("string"): sel = AtomSel(sel)
    ret=0.
    for atom in sel:
        ret += sa.singleAtom(atom)
    return ret
    
def dumpVMDRadii(sa):
    """dump radius info in VMD format to file named 'radius.vmd'
    """
    radii={}
    from atomSel import AtomSel
    for atom in AtomSel('known'):
        radii[ atom.atomName() ] = sa.radius(atom)
        pass
    ofile = open("radius.vmd","w")
    for name in list(radii.keys()):
        ofile.write('[atomselect top "name %s"] set radius %f\n' %
                    (name, radii[name]))
        pass
    return
