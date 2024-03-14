
"""
add docs
"""

import protocol
protocol.addPseudoResName("ANI")

def create_EnsWeights(name,
                      members=None,
                      baseSize=None,
                      derivedSpecs=[],
                      sim=None):
    """
    create an EnsWeights object with the given name to allow optimization of
    ensemble weights for specified potential terms.

    The optional members argument should be a list of ensemble indices to
    optimize.

    The baseSize and derivedSpecs arguments are used to specify an
    interrelationship between ensemble weights. If baseSize is omitted,
    derivedSpecs specifies a set weights dependent on the
    previous len(members)-len(derivedSpecs) weights. If baseSize is specified,
    the number of independent weights is baseSize, and derivedSpecs should be of
    len(members).

    derivedSpecs is a list of tuples with three arguments each.The
    populations of the dependent (derived) w_k is

       w_k = A_k + B_k * w_b + w_b^T * C_k * w_b

    where w_b are base weights, A_k is a constant, B_k a vector and C_k a
    matrix. The parameters A,B and C are specified by the derivedSpecs entry as a
    sequence of three values (A,B,C). The specifications for B and C are in sparse
    notation by index of nonzer entries, so the second gives the nonzero members
    of the B vector in the equation above in the form [index, value]. The third
    entry gives the nonzero values of the C matrix in the form
    [index1, index2, value].

    """

    from simulation import currentSimulation
    if not sim: sim = currentSimulation()

    import ensembleSimulation
    if not members:
        try:
            size = ensembleSimulation.fromSimulation(sim).size()
        except:
            size=1
            pass
        members = list(range(size))
        pass
    
    esim = ensembleSimulation.fromSimulation(sim)

    lderivedSpecs=[]
    from cdsMatrix import RMat
    from cdsVector import CDSVector_double as RVec
    if not baseSize: #auto-fill first baseSize elements if baseSize arg omitted
        baseSize = len(members) - len(derivedSpecs)
        for i in range(baseSize):
            lderivedSpecs.append( (0., ((i,1.0),), ()) )
            pass
        pass
    derivedSpecs = lderivedSpecs + derivedSpecs 
    lderivedSpecs=[]
    for (i,a,b) in derivedSpecs:
        la = RVec(baseSize,0.)
        for index,val in a:
            la[index] = val
            pass
        lb = RMat(baseSize,baseSize,0.)
        for m,n,val in b:
            lb[m,n] = val
            pass
        lderivedSpecs.append( (i,la,lb) )
        pass

    phiList=[]
    atoms = addPseudoAtoms(esim,baseSize)
    from bondAngle import BondAngle
    for i in range(baseSize-1):
        phi = BondAngle(atoms[0],atoms[2*i+1],atoms[2*i+2])
        phiList.append(phi)
        pass

    from ensWeights import EnsWeights
    eweights = EnsWeights(name,members,phiList,baseSize,lderivedSpecs,sim)

#    # set initial atom positions - weights will be normalized
#    eweights.setWeights([1]*len(members))

    curWeights=[]
    
    for i in sorted(members):
        curWeights.append(esim.weight(i))
        pass        
        
    eweights.setWeights(curWeights)
    registerTerm(eweights)
    
    return eweights

segid="ENSW"
resid=1100
lastp=0
initialized=False
def addPseudoAtoms(esim,size):
    global resid
    import ensWeightsTools

    if size==1:
        return

    #these pseudoatoms will live in ensemble member 1 (second member)
    emem = esim.members(1)


    from xplorSimulation import getXplorSimulation
    xSim = getXplorSimulation()
    import protocol
    protocol.initParams("axis")

    if not initialized:
        xSim.fastCommand("topology %s end" % patch)  #only do this once
        xSim.syncFrom()
        ensWeightsTools.initialized=True
        pass
        

    from atomSel import AtomSel
    retSel = AtomSel("segid %s and resid %s" % (segid,resid),emem)
    if len( retSel ) > 0:
        raise Exception("segid %s and resid %d already exist." % (segid,resid))

    if esim.member().memberIndex()==1:
        cmd = psfTemplate.replace('_n__','%-4d'%resid)
        cmd = cmd.replace('SGMT','%-4s'%segid)
#        print 'cmd:', size, cmd
        xSim.fastCommand(cmd)

        for i in range(1,size-1):
        
            xSim.fastCommand(
                'patch ewea reference=1=(segid "%s" and resid %d) end' %
                         (segid,resid))
            xSim.fastCommand(
                'vector do ( name = P%d ) ( ATOM "%s" "%d" PY )' %
                (i,segid,resid))
            xSim.fastCommand(
                'vector do ( name = OO%d ) ( ATOM "%s" "%d" OOY )' %
                (i,segid,resid))
            pass
        pass

    from simulation import currentSimulation
    xSim.syncFrom()
    currentSimulation().sync()
    esim.barrier()
    retSel.reevaluate()

#    print emem.name()
#    print len(AtomSel("all"))
#    print AtomSel("all")[-1].string()
#    print len(AtomSel("segid %s and resid %s" % (segid,resid),emem))
#    print len(AtomSel("segid %s and resid %s" % (segid,resid),xSim))
    
    if esim.member().memberIndex()==1:

        from vec3 import Vec3
        retSel[0].setPos( Vec3(1,0,0) )
        for i in range(size-1):
            retSel[2*i+1].setPos( Vec3(0,0,0) )
            retSel[2*i+2].setPos( Vec3(0,1,0) )
            pass
        pass

    esim.barrier()

    protocol.updatePseudoAtoms()


    resid += 1

    return retSel

patch=r"""
mass YYY 100
mass OOO 100

presidue ewea
  group
    add       atom 1OOY type=OOO end
    add       atom 1PY  type=YYY end

  add bond 1OOY 1PY

  add angle 1X 1OOY 1PY
end
"""

# the string _n__ is replaced by the residue number
# the string SGMT is replaced by the segment name
psfTemplate = """
structure
PSF

       1 !NTITLE
 REMARKS   ensWeightsTools.py: auto-generated structure parameters

       3 !NATOM
       1 SGMT _n__ ANI  X    XXX    0.000000E+00   10.0000           0
       2 SGMT _n__ ANI  OO   OOO    0.000000E+00   10.0000           0
       3 SGMT _n__ ANI  P0   YYY    0.000000E+00   10.0000           0

       2 !NBOND: bonds
       1       2       2       3

       1 !NTHETA: angles
       1       2       3

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


def massSetup(list=[],axisMass=1000):
    """
    appropriately setup pseudo-atom masses to axisMass.

    if list is not specified, then pseudoatoms associated with
    all registered EnsWeight objects are configured. 
    """
    try:
        len(list)
    except TypeError:
        list = [list]
        pass

    from simulation import currentSimulation
    global registeredTerms
    if not list:
        try:
            list = list(registeredTerms[currentSimulation().lookupID()].values())
        except:
            pass
        pass

    for t in list:
        for ba in t.angles():
            for a in [ba.atom0(),ba.atom1(),ba.atom2()]:
                a.setMass(axisMass)
                pass
            pass
        pass
    return

def topologySetup(ivm,list=[]):
    """
    configure the given <m ivm>.IVM object's topology setup using the
    freedom string for each EnsWeight in list.
    This function should be called prior to
    <m ivm>.IVM.autoTorsion() or <m protocol>.torsionTopology()

    The freedom specifier should be one of the following keywords:
       fix
       bend
       ignore

     if list is not specified, then pseudoatoms associated with
     all registered EnsWeights objects are configured. 
     """

    try:
        len(list)
    except TypeError:
        list = [list]
        pass

    from simulation import currentSimulation
    global registeredTerms
    if len(list)==0:
        try:
            list = registeredTerms[currentSimulation().lookupID()].values()
        except:
            pass
        pass

    
    import ensembleSimulation
    if ensembleSimulation.singleThread(1):
        for p in list:

            keyword = p.freedom()
            if keyword=="ignore":
                continue

            angles = p.angles()

            if len(angles)<1: continue

            xAtom = angles[0].atom0()
#        print 'len', len(p.angles())
#        print 'angle0', p.angles()[0].value()
#        print 'angle', angles[0].value()

            ivm.fix(xAtom)

#        print 'topologySetup: processing:',p.instanceName()
#        print oAtom, type(oAtom)
#        print oAtom.string()
#        print angles[0].atom0()

            for i in range(len(angles)):
                oAtom = angles[i].atom1()
                yAtom = angles[i].atom2()
#            print angles[i].atom2()
                #print yAtom.string()
                if keyword=='fix':
                    ivm.fix( (oAtom,yAtom ))
                elif keyword=='bend' or keyword=='vary':
                    ivm.group( (oAtom,yAtom) )
                    ivm.hinge("bend", oAtom,xAtom,yAtom)
                    pass
                else:
                    raise Exception("term: %s: bad freedom keyword: %s" %
                                    (p.instanceName, p.freedom()))
                pass
        pass
    ensembleSimulation.multiThread()
    return
    
def optimizeWeights(ensWeights,
                    potList   ,
                    context=None):
    """optimize ensemble weights given by the ensWeights argument using the
    supplied potList, and reset target values used during
    simulated annealing"""
    from ivm import IVM
    minimizeWeights=IVM()
    minimizeWeights.fix("not segid ENSW")
    import protocol
    protocol.torsionTopology(minimizeWeights)

    from simulationTools import convertToPotList
    potList = convertToPotList(potList)

    print('initial energy, weights:', potList.calcEnergy(),ensWeights.weights())
    if context: context()
    for i in range(5):
        protocol.initMinimize(minimizeWeights,
                              numSteps=100,
                              potList=potList)
        minimizeWeights.run()
        weights=ensWeights.weights()
        if context: context()
        print('optimized energy,weights:', potList.calcEnergy(),weights)
        pass
    ensWeights.setTargetWeights( [ weights[k] for k in
                                   sorted(weights.keys()) ] )
    return

registeredTerms={}
def registerTerm(term):
    """
    add the given EnsWeights object to a list associated with its
    Simulation. These objects will be automatically processed by topologySetup
    and massSetup.
    """
    global registeredTerms

    sim = term.simulation()
    id = sim.lookupID()
    if not id in registeredTerms: registeredTerms[id] = {}

    # use a instanceName-keyed dictionary to avoid referring to
    # dead terms. This is a hack needed until objects are deleted properly
    # (then they can unregister themselves).
    registeredTerms[id][term.instanceName()]=term
    return

def getEnsemblePots(potList):
    ret=[]
    from simulationTools import flattenPotList
    for term in flattenPotList(potList):
        if hasattr(term,'ensWeight') and hasattr(term,'simulation'):
            if term.simulation() and term.simulation().size()>1:
                ret.append(term)
            pass
        pass
    return ret
    

def analyze(potList):
    "perform analysis of EnsWeights terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    ensWeightsList = getPotTerms(potList,'EnsWeights')

    #FIX should replace the next line with all EnsemblePots
    supportedTerms = getEnsemblePots(potList)

    for term in supportedTerms:
        ensWeightsList += [ew[0] for ew in term.getEnsWeights()]
        pass

    if ensWeightsList:
    
        instanceNames = [x.instanceName() for x in ensWeightsList]
        instanceNames.sort()
    
        ret+= "%9s  Weights\n" % \
              (" " ,)
    
        prevName=None
        for name in instanceNames:
            if name==prevName: continue
            prevName=name
            ensWeights = [x for
                          x in ensWeightsList if x.instanceName()==name][0]
    
            print("EnsWeights term: ", ensWeights.instanceName())
            print(ensWeights.info())
    
            weights=sorted(ensWeights.weights().items())
            pretty=''
            for weight in weights:
                pretty += " %d: %6.4f," % weight
                pass
    
            ret += "%-12s  {%s }\n" % (name ,pretty[:-1])
            pass
        pass
    
    if len(supportedTerms):
        ret += "\nper-term Ensemble Weights for this structure:\n"
        for term in supportedTerms:
            if term.potName()=="EnsWeights":  # these are already reported
                continue
            eIndex = term.simulation().member().memberIndex()
            ret += "%20s: Ensemble Weight: %8.6f\n" % (term.instanceName(),
                                                       term.ensWeight(eIndex))
            pass
        pass
    
    return ret

import simulationTools
simulationTools.registerTerm(analyze,"Ensemble Weights","EnsWeights",
r"""Ensemble weights (populations) reported for each energy term.

The <m ensWeights>.EnsWeights term allows for varying/optimizing ensemble
weights. For each ensemble member the associated weight is printed. Also, for
each energy term, the associated ensemble weight is reported. This is to
emphasize that the appropriate ensemble weight can be different for different
energy terms, perhaps due to different experimental conditions. Terms averaged
using an <m avePot>.AvePot should all have equal weight.
""")


def getWeights(term):
    return list(term.weights().values())

simulationTools.registerExtraStats("EnsWeights","Weights",getWeights)
