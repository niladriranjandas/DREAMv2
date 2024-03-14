"""
Tools to help create, manipulate and analyze distance restraints measured
between an atom and a plane.
"""

import protocol
protocol.addPseudoResName("ANI")


def create_PlaneDistPot(name,restraintsFile="",
                        oAtom=None,xAtom=None,yAtom=None,
                        replicatePlaneAtoms=True,
                        A=None, B=None, C=None, D=None,
                        restraints=""
                        ):
    """ create a <m planeDistPot>.PlaneDistPot using a plane defined by the
    equation
        A x + B y + C z + D = 0
    or
        by three atoms.

    If one of A, B, C, or D is not specified, it is taken to be zero.

    If the atom-type definition is used, replicatePlaneAtoms is consulted
    as to whether the plane remains fixed to the plane-defining atoms or
    the plane is decoupled from those atoms.

    restraintsFile is a file name of a restraints assignment table, while
    restraints is a string containing these restraints.
    
    """

    import simulationTools
    simulationTools.gcRegisteredTopoTerms()

    if replicatePlaneAtoms or not oAtom:
        (oAtom,xAtom,yAtom) = addPlaneAtoms(oAtom,xAtom,yAtom,
                                            A, B, C, D)
    else:
        from selectTools import convertToAtom
        oAtom = convertToAtom(oAtom)
        xAtom = convertToAtom(xAtom)
        yAtom = convertToAtom(yAtom)

        pass

        

    if restraintsFile:
        restraints += open(restraintsFile).read()
        pass
    
    from planeDistPot import PlaneDistPot
    print(name,oAtom.string(),xAtom.string(),yAtom.string(),len(restraints))
    pot =  PlaneDistPot(name,
                        oAtom,xAtom,yAtom,
                        restraints)
    if not replicatePlaneAtoms:
        pot.setFreedom('ignore')
        pass

    terms=simulationTools.getRegisteredTopoTerms(pot.potName())
    if name in [term.instanceName() for term in terms]:
        raise Exception("name %s is already in use" % name)

    from simulationTools import registerTopoTerm
    registerTopoTerm(pot)
    protocol.updatePseudoAtoms()

    return pot

default_segid="PLAN"
default_resid=1100
def addPlaneAtoms(oAtom=None,xAtom=None,yAtom=None,
                  A=None, B=None, C=None, D=None,
                  resid=None,segid=None):
    """
    create the atoms defining the restraint plane.

    Either (oAtom, xAtom, yAtom ) or (A,B,C,D) should be specified.
    """

    global default_segid, default_resid
    if not segid:
        segid = default_segid
        

    from xplorSimulation import getXplorSimulation
    xSim = getXplorSimulation()
    import protocol
    protocol.initParams("axis")

    from atomSel import AtomSel
    if not resid:
        while len( AtomSel("segid %s and resid %d" %
                           (segid,default_resid),xSim) )>0:
            default_resid +=1 
            pass
        resid = default_resid
        pass

    if len( AtomSel("segid %s and resid %s" % (segid,resid),xSim) ) > 0:
        raise Exception("segid %s and resid %d already exist." % (segid,resid))

    cmd = psfTemplate.replace('_n__','%-4d'%resid)
    cmd = cmd.replace('SGMT','%-4s'%segid)
    xSim.fastCommand(cmd)
    xSim.syncFrom()
    protocol.updatePseudoAtoms(xSim)

    xSim.sync()
    from simulation import currentSimulation
    currentSimulation().sync()


    from vec3 import unitVec, dot, cross, Vec3, norm
    from selectTools import convertToAtom
    if oAtom:
        oAtom = convertToAtom(oAtom,xSim)
        xAtom = convertToAtom(xAtom,xSim)
        yAtom = convertToAtom(yAtom,xSim)

        p1=oAtom.pos()
        p2=xAtom.pos()
        p3=yAtom.pos()
        n =  unitVec(cross(p2-p1,p3-p1))
        (A,B,C) =n
        D = -dot(n,p1)
        pass
        
    # set initial coords for the new atoms
    # closest point on plane to center of mass
    n = unitVec(Vec3(A if A!=None else 0,
                     B if B!=None else 0,
                     C if C!=None else 0,))
    D = D if D!=None else 0
    smallVal=1e-5
    if norm(n)<smallVal:
        raise Exception("plane normal magnitude too small")

    from atomAction import centerOfMass
    pCM = centerOfMass("known and not PSEUDO")
    
    pO = pCM - (dot(n,pCM)+D) * n
    x = Vec3(1,0,0); y = Vec3(0,1,0) ; z = Vec3(0,0,1)

    xDir = cross(n,x)
    yDir = cross(n,y)
    zDir = cross(n,z)
    if norm(xDir)>smallVal:
        yDir = cross(n,xDir)
        pX = pO + unitVec(xDir)
        pY = pO + unitVec(yDir)
    elif norm(yDir)>smallVal:
        yDir = cross(n,zDir)
        pX = pO + unitVec(yDir)
        pY = pO + unitVec(zDir)
    else:
        raise Exception("error determining plane atom positions")
                
    resSel = "segid %s and resid %d" % (segid,resid)

    oAtom = AtomSel(resSel + " and name OO" ,xSim)[0]
    xAtom = AtomSel(resSel + " and name X " ,xSim)[0]
    yAtom = AtomSel(resSel + " and name Y " ,xSim)[0]
    oAtom.setPos( pO )
    xAtom.setPos( pX )
    yAtom.setPos( pY )

    currentSimulation().sync()

    return (oAtom,xAtom,yAtom)

# the string _n__ is replaced by the residue number
# the string SGMT is replaced by the segment name
psfTemplate = """
structure
PSF

       1 !NTITLE
 REMARKS   planeDistTools.py: auto-generated structure parameters

       3 !NATOM
       1 SGMT _n__ ANI  OO   OOO    0.000000E+00   10.0000           0
       2 SGMT _n__ ANI  X    XXX    0.000000E+00   10.0000           0
       3 SGMT _n__ ANI  Y    YYY    0.000000E+00   10.0000           0

       2 !NBOND: bonds
       1       2       1       3

       1 !NTHETA: angles
       2       1       3

       0 !NPHI: dihedrals


       0 !NIMPHI: impropers


       0 !NDON: donors


       0 !NACC: acceptors


       0 !NNB

       0       0       0       0

       1       0 !NGRP
       0       0       0
end       
"""

def massSetup(list=[],axisMass=300):
    """
    appropriately setup plane-atom masses to axisMass.

    if list is not specified, then pseudoatoms associated with
    all registered PlaneDistPot objects are configured. 
    """
    try:
        len(list)
    except TypeError:
        list = [list]
        pass

    from simulationTools import getRegisteredTopoTerms
    if not list: list = getRegisteredTopoTerms("PlaneDistPot")

    for t in list:
        for a in [t.oAtom(),t.xAtom(),t.yAtom()]:
            a.setMass(axisMass)
            pass
        pass
    return

def topologySetup(ivm,list=[]):
    """
    configure the given <m ivm>.IVM object's topology setup using the
    freedom string for each PlaneDistPot in list.
    This function should be called prior to
    <m ivm>.IVM.autoTorsion() or <m protocol>.torsionTopology()

    The freedom specifier should be one of the following keywords:
       fix
       group
       translate
       rotate
       free
       ignore

     if list is not specified, then pseudoatoms associated with
     all registered PlaneDistPot objects are configured. 
     """

    try:
        len(list)
    except TypeError:
        list = [list]
        pass

    from simulationTools import getRegisteredTopoTerms
    if len(list)==0: list = getRegisteredTopoTerms("PlaneDistPot",
                                                   ivm.simulation)

    for p in list:
        planeSel = "id %d or id %d or id %d" %( p.oAtom().index()+1,
                                                p.xAtom().index()+1,
                                                p.yAtom().index()+1)
        keyword = p.freedom().split()[0]
        if keyword=='fix':
            ivm.fix(planeSel)
        elif keyword=='group':
            ivm.group(planeSel)
        elif keyword=='translate':
            ivm.hinge("translate",planeSel)
        elif keyword=='rotate':
            ivm.hinge("rotate",planeSel)
        elif keyword=='free':
            ivm.breakAllBondsIn(planeSel)
            ivm.setBaseAtoms(planeSel)
        elif keyword=='ignore':
            pass
        else:
            raise Exception("term: %s: bad freedom keyword: %s" %
                            (p.instanceName, p.freedom()))
        pass
    return
    
        
    

#FIX: add
# torsionTopology, cartesianTopology will call this routine -
#   need to add hook to make sure this happens.
# 

def analyze(potList):
    """perform analysis of PlaneDistPot terms and return nicely formatted
    summary"""

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'PlaneDistPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret += "%-9s  %6s  %6s \n" % \
           ( "", "RMS", "Viols")

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        print(term.showViolations())
        print(term.info())

        ret += "%-9s  %6.3f  %6d\n" % \
               (name , term.rms(), term.violations() )
        pass
    
    return ret


import simulationTools
simulationTools.registerTerm(analyze,"PlaneDist terms","PlaneDist",
r"""
For each term the following are reported:
  RMS   - root mean square deviation between calculated and target values.
  Viols - number of violations.
""")
