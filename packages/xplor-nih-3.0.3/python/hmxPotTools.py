
def create_HMXPot(name,
                  fileName=None,
                  data="",
                  rnaSegid=None,
                  selection="not PSEUDO",
                  l=1.4,r=5,
                  radiusScale=2.3,
                  ):
    """
    data is input via lines of comma-separated values having the columns
    resid, resname, contrib1, contrib2, hmx
    Lines starting with a # are skipped.
    """

    if fileName: data += open(fileName).read()
    
    stringData=[line.split(',') for line in data.split('\n')
            if not (line.startswith('#') or not ',' in line)]
    
    data={}
    resids=[]
    from atomSel import AtomSel
    from vec3 import unitVec
    dummyAtoms=[]
    for resid,resname,s1,s2,hmx,so in stringData:
        resid = int(resid)
        resids.append(resid)
        data[resid] = 100*float(hmx)

        sel = "name O2' and resid %d" % resid
        if rnaSegid: sel += 'and segid "%s"' % rnaSegid
        try:
            o2pAtom = AtomSel(sel)[0]
        except IndexError:
            raise Exception("Could not find atom for "+sel)
        o2pPos = o2pAtom.pos()
        segid = o2pAtom.segmentName()


        c2pPos = AtomSel("name C2' and resid %d and segid \"%s\""%
                         (resid,segid)
                         )[0].pos()
        dummyAtoms.append( addDummyAtom(resid,
                                        segid=segid,
                                        pos=o2pPos +
                                             l * unitVec(o2pPos-c2pPos)) )
        pass

    import xplor
    xplor.command("param nonbonded DUMY   0.145   %f  0.145  %f end" %
              tuple([r]*2))

    from atomSel import union
    selection = union(selection,"resname HMXA")
    ret = HMXPot(name,selection,dummyAtoms,resids, data, radiusScale)
    registerTerm(ret)
    return ret



    

from pyPot import PyPot
class HMXPot(PyPot):
    def __init__(s,name,selection,dummyAtoms,resids, data, repel=2.1):
        """
        """

        PyPot.__init__(s,name)
        s.resetPotName('HMXPot')

        s.resids = resids
        s.data = data
        s.corr = -1

        s.obs = []
        for resid in resids:
            s.obs.append( data[resid] )
            pass

        selPairs=[]
        s.restraints_ = []
        simulation = selection.simulation()
        from atomSel import AtomSel, intersection
        for atom in dummyAtoms:
            name = atom.atomName()
            resid = atom.residueNum()
            segid = atom.segmentName()
            o2pAtom=AtomSel('name O2\' and segid "%s" and resid %d' %(segid,
                                                                      resid),
                            simulation)[0]
            
            interacting = intersection(selection,
                                       'not (resid %d:%d and segid "%s")' %
                                       (resid-1,resid+1,segid))
            interacting = intersection(interacting,"not resname HMXA")
            selPairs.append( ('ATOM "%s" %d %s' % (segid,resid,name),
                              interacting) )
            s.restraints_.append( Restraint(s,segid,o2pAtom.residueName(),
                                            resid,data[resid]) )
            pass                     
        
        from repelPotTools import create_RepelPot
        s.repel=create_RepelPot('repel',selection=selection,
                                selPairs=selPairs,repel=repel)
        s.repel.setCalcSelPairDerivs(True)

        s.selPairMap={} # map resid to selection pair index
        for i,sel in enumerate([pair.a for pair in s.repel.selectionPairs()]):
            for resid in resids:
                if sel[0].residueNum()==resid:
                    s.selPairMap[resid] = i
                    pass
                pass
            pass


        import types
        if type(selection) == bytes:
            from atomSel import AtomSel
            selection = AtomSel(selection)
            pass
        return

    def __del__(self,destroy=0):
        PyPot.__del__(self)
        return

    def simulation(s): return s.repel.defaultSelection().simulation()
    def selection(s): return s.repel.defaultSelection()
    def numRestraints(s): return len(s.restraints_)
    def restraints(s): return s.restraints_
    def scaleFactor(s): return s.scaleFactor_
    def offset(s): return s.offset_
    
    def calcEnergy(s):
        return s.calcEnergyAndDerivList()

    def correlation(s):
        return s.corr

    def calcEnergyAndDerivList(s,derivs=None):

        calcd=[]
        for resid in s.resids:
            calcd.append( s.repel.selectionPairEnergy( s.selPairMap[resid] ) )
            pass

        sim = s.selection().simulation()

        from cdsVector import CDSVector_double, sum
        calcd = CDSVector_double(calcd)
        obs  = CDSVector_double(s.obs)
        
        ave_obs   = sum(obs) / len(obs)
        ave_calcd = sum(calcd) / len(calcd)
        num = sum( (calcd-ave_calcd) * (obs - ave_obs) )
        denom1  = sum( (calcd-ave_calcd)**2 )
        sumObs2 = sum( (obs-ave_obs)**2 )
        
        from math import sqrt
        denom = sqrt(denom1 * sumObs2)
        s.corr = num / denom

        if derivs:
            for i,resid in enumerate(s.resids):
                dC = (obs[i]-ave_obs - num/denom1*(calcd[i]-ave_calcd)) /denom
                contrib = dC * s.repel.selectionPairDerivs( s.selPairMap[resid] )
                derivs[sim] -= s.scale() * contrib
                pass
            pass
        else:
            pass
            
        from cdsMatrix import CDSMatrix_double, inverse
        M=CDSMatrix_double([ [sum(calcd**2), sum(calcd)],
                             [sum(calcd)   , len(calcd)]] )
        b=CDSVector_double( [sum(obs*calcd), sum(obs)] )
        slope,intercept=tuple(inverse(M)*b)
        s.scaleFactor_ = slope
        s.offset_ = intercept
        calcd *= s.scaleFactor()
        calcd += s.offset()
        for val,restraint in zip(calcd,s.restraints()):
            restraint.calcd = val
            pass
        # FIX: check that the correlation hasn't changed

        return s.scale() * (1.0 - s.corr)

    def showRestraints(s):
        ret  = "\n     HMX Restraints: term: %s\n\n" % s.instanceName()
        ret += "Index     Segid  Resname Resid       "
        ret += "Calcd    Obs     Diff\n"
        ret += "-"*79
        ret += '\n'
        cnt=0
        for r in s.restraints():
            ret += "%3d "%cnt
            cnt += 1
            ret += "     %4s     %4s  %5d " % (r.segid,
                                                 r.resname,
                                                 r.resid)
            ret +=  "     %7.2f %7.2f %7.2f\n" %(r.calcd,r.obs,
                                              r.diff())
            pass
        return ret
    pass


class Restraint:
    def __init__(s,pot,segid,resname,resid,obs):
        s.pot=pot
        s.segid = segid
        s.resname = resname
        s.resid = resid
        s.obs=obs
        s.calcd=0
        return
    def name(s):
        ret = 'segid: %s  resid: %d' % (s.segid,s.resid)
        return ret
    def diff(s):
        return s.calcd - s.obs
    def violated(s):
        if abs(s.diff()) > s.pot.threshold():
            return True
        else:
            return False
        pass
    pass


default_segid="ADUC"
def addDummyAtom(resid,pos,segid=None):
    """
    """

    global default_segid
    if not segid:
        segid = default_segid
        pass

    from xplorSimulation import getXplorSimulation
    xSim = getXplorSimulation()
    import protocol
    protocol.initParams("axis")


    from atomSel import AtomSel

    if len( AtomSel("name OO and segid %s and resid %s" % (segid,resid))) > 0:
        raise Exception("segid %s and resid %d already exist." % (segid,resid))

    atom = xSim.addAtom(atomName='OO',
                        residueName='HMXA',
                        residueNum=resid,
                        chemicalType='DUMY',
                        segmentName=segid,
                        pos=pos,
                        isPseudo=True)

    from simulation import currentSimulation
    currentSimulation().sync()

    return atom

def topologySetup(ivm,list=[]):
    """
    Configure <m ivm>.IVM object s.t. the dummy atom remains in its position
    relative to the O2' and C2' atoms.
    """

    try:
        len(list)
    except TypeError:
        list = [list]
        pass

    from simulation import currentSimulation
    if not list:
        try:
            list = getRegisteredTerms( currentSimulation() )
        except:
            pass
        pass

    from atomSel import AtomSel
    for p in list:
        for r in p.restraints():
            simulation = p.simulation()
            atom    = AtomSel('ATOM "%s" %d OO'   % (r.segid,r.resid),
                              simulation)[0]
            o2pAtom = AtomSel('ATOM "%s" %d O2\'' % (r.segid,r.resid),
                              simulation)[0]
            c2pAtom = AtomSel('ATOM "%s" %d C2\'' % (r.segid,r.resid),
                              simulation)[0]
            ivm.group((o2pAtom,c2pAtom,atom))
            pass
        pass
    return
import protocol
protocol.addTopologySetup(topologySetup,"hmxPotTools")

registeredTerms={}
def registerTerm(term):
    """
    add the given HMXPot object to a list associated with its esim.
    These objects will be automatically processed by topologySetup and
    massSetup.
    """
    global registeredTerms

    esim = term.simulation()
    id = esim.lookupID()
    if not id in registeredTerms: registeredTerms[id] = []

    registeredTerms[id].append(term)
    return
def unRegisterTerm(term):
    "remove the given HMXPot object from a list associated with its esim"
    
    global registeredTerms

    esim = varTensor.simulation()
    id = esim.lookupID()
    
    registeredTerms[id].remove(term)
    return

def getRegisteredTerms(sim):
    """
    return a list of all HMXPot objects registered to
    the specified simulation.
    """
    global registeredTerms
    ret=[]
    try:
        ret = registeredTerms[sim.lookupID()]
    except:
        pass
    try:
        if sim.type()=="EnsembleSimulation":
            from ensembleSimulation import EnsembleSimulation_currentSimulation
            ret += registeredTerms[
                EnsembleSimulation_currentSimulation().member().lookupID()]
            pass
        pass
    except:
        pass
    return ret




def analyze(potList):
    "perform analysis of HMXPot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'HMXPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret+= "%9s  %7s\n" % \
     (" " , "Correlation")
    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        ret += "%9s  %7.3f\n" % (name,term.correlation())

        print(term.showRestraints())

        print("scaleFactor: %8.3e"% term.scaleFactor())
        print("offset:      %7.4f"% term.offset())
        pass

    
    return ret
#
# print print analysis
import simulationTools
simulationTools.registerTerm(analyze,
                             "HMX Analysis","HMXPot",
r"""
For each term, the correlation between back-calculated and observed values is
printed.
""")
        
def correlation(potTerm): return potTerm.correlation()
simulationTools.registerExtraStats("HMXPot","correlation",correlation)
