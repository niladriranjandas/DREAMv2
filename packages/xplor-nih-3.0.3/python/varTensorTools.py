"""Tools for manipulation of tensors represented with 7 pseudo-atoms.

This module presents functions to simplify the creation, manipulation and
analysis of <m varTensor>.VarTensor objects, which are used to represent
orientation tensors.
"""


from atomSel import AtomSel
from varTensor import VarTensor
from cdsVector import CDSVector_double
import ensembleSimulation

import protocol
protocol.addPseudoResName("ANI")


def create_VarTensor(name,axis=0,esim=0):
    """Create a <m varTensor>.VarTensor object.

    name (a string) is its assigned name.  axis is a sequence with the segment
    id (segid; a string) and residue number (resid; a number), in that order,
    associated with the axis atoms.  esim is an optional ensemble simulation
    argument.

    """
    from simulation import currentSimulation
    if not esim: esim = currentSimulation()

    if not axis: axis = addAxisAtoms()
    
    selStr = "resid %d" % axis[1]
    if axis[0]!="":
        selStr += " and segid %s" % axis[0]
        pass

    import simulationTools
    simulationTools.gcRegisteredTopoTerms()

    ret =  VarTensor(name,
                     AtomSel(selStr + " and name OO ",esim),
                     AtomSel(selStr + " and name X  ",esim),
                     AtomSel(selStr + " and name Y  ",esim),
                     AtomSel(selStr + " and name Z  ",esim),
                     AtomSel(selStr + " and name OO2",esim),
                     AtomSel(selStr + " and name PA1",esim),
                     AtomSel(selStr + " and name PA2",esim),
                     esim)

#    terms=simulationTools.getRegisteredTopoTerms(ret.potName(),esim)
#    if name in [term.instanceName() for term in terms]:
#        raise Exception("name %s is already in use" % name)

    from simulationTools import registerTopoTerm
    registerTopoTerm(ret)
    protocol.updatePseudoAtoms()

    return ret

def deletePseudoAtoms(obj):
    """
    Delete the pseudoatoms associated with the specified VarTensor obj.
    """
    esim = obj.simulation()
    from xplorSimulation import fromSimulation
    xsim = fromSimulation(esim.member().subSim())
    outputState=xsim.disableOutput()
    for sel in (obj.oAtomSel() ,
                obj.xAtomSel() ,
                obj.yAtomSel() ,
                obj.zAtomSel() ,
                obj.o2AtomSel(),
                obj.p1AtomSel(),
                obj.p2AtomSel()):
        #XplorSimulation specific - we need this
        #  because simulation.deleteAtoms will fail due to the presence of
        #  dependents.
        xsim.wrap().deleteAtoms( AtomSel(sel,xsim).indices() )
        pass
    xsim.enableOutput(outputState)
    protocol.updatePseudoAtoms()
    
    return

def saupeToVarTensor(mat,Dmax=1):
    """given a Saupe matrix as a <m mat3>.SymMat3, generate a 
    <m varTensor>.VarTensor object with whose orientation axis, Da, and Rh
    correspond to the input matrix. The Saupe matrix can be premultiplied by
    Dmax or, optionally, it can be supplied as a separate argument.

    If an <m ensembleSimulation>.EnsembleSimulation is active, the input
    tensors can be different for each ensemble member.

    Note that the generated VarTensor object is accompanied by pseudoatoms,
    and the term is registered within varTensorTools, so that full
    cleanup of the VarTensor object requires a call to
    unregisterTerm and deletePseudoAtoms.
    """

    if not "SymMat3" in str(type(mat)):
        raise Exception("mat must be a mat3.SymMat3 object")
    
    ten=create_VarTensor('ten')
    ten.setEnsembleAxis(True)
    ten.setEnsembleDa(True)
    ten.setEnsembleRh(True)

    from mat3 import eigen
    eigenPairs = eigen(mat)
    
    xInd = 0; yInd=1; zInd=2;
            
    if abs(eigenPairs[zInd].value()) < abs(eigenPairs[yInd].value()):
        zInd,yInd = yInd,zInd
        pass
    if abs(eigenPairs[zInd].value()) < abs(eigenPairs[xInd].value()):
        zInd,xInd = xInd,zInd
        pass
    if abs(eigenPairs[yInd].value()) < abs(eigenPairs[xInd].value()):
        yInd,xInd = xInd,yInd
        pass
    
    Azz = 2.0/3.0 * eigenPairs[zInd].value()
    Ayy = 2.0/3.0 * eigenPairs[yInd].value()
    Axx = 2.0/3.0 * eigenPairs[xInd].value()
    
    Aa = 3 * Azz / 2.
    Ar = (Axx-Ayy)
            
    Da = 0.5 * Aa * Dmax
    
    Rh = Ar/Aa 
    from vec3 import Vec3
    xPos = Vec3(eigenPairs[xInd].vector())
    yPos = Vec3(eigenPairs[yInd].vector())
    zPos = Vec3(eigenPairs[zInd].vector())

    from vec3 import dot, cross
    zPos *= dot(cross(xPos,yPos) , zPos) # make right-handed

    oPos = ten.oAtom().pos()
    ten.xAtom().setPos( xPos+oPos )
    ten.yAtom().setPos( yPos+oPos )
    ten.zAtom().setPos( zPos+oPos )
    ten.setDa(Da)
    ten.setRh(Rh)

    return ten

def registerExptToVarTensor(varTensor,expt):
    """
    add the term expt to the <m varTensor>.VarTensor's expts list.
    """
    varTensor.expts().append(expt)
    varTensor.expts()[-1].thisown=False
    return

from selectTools import convertToAtom
def configIVM(p,ivm):
    """first argument is an potential term containing a varTensor and the
    second is an IVM. This routine does the required initial topology setup."""

    # remove the o-x and o-z bonds
    ivm.breakBond( [convertToAtom(p.xAtom(),ivm.simulation),
                    convertToAtom(p.oAtom(),ivm.simulation)] )

    ivm.breakBond( [convertToAtom(p.zAtom(),ivm.simulation),
                    convertToAtom(p.oAtom(),ivm.simulation)] )

    memberIndex= p.simulation().member().memberIndex()

    if memberIndex==0 or p.ensembleAxis():
        ivm.setBaseAtoms( ivm.baseAtoms() +
                          [convertToAtom(p.zAtom(),ivm.simulation)] )
        ivm.group( (p.xAtom(), p.yAtom(), p.zAtom()) )
        ivm.hinge( "rotate", (p.xAtom(), p.yAtom(),  p.zAtom()) )
    else:
        ivm.fix( [convertToAtom(p.xAtom(),ivm.simulation),
                  convertToAtom(p.yAtom(),ivm.simulation),
                  convertToAtom(p.zAtom(),ivm.simulation)] )
        pass

    if memberIndex>0 and not p.ensembleDa():
        ivm.fix( [convertToAtom(p.oAtom() ,ivm.simulation),
                  convertToAtom(p.p1Atom(),ivm.simulation)] )
        pass

    if memberIndex>0 and not p.ensembleRh():
        ivm.fix( [convertToAtom(p.o2Atom() ,ivm.simulation),
                  convertToAtom(p.p2Atom(),ivm.simulation)] )
        pass

    return
        

def configIVM_fix(p,ivm):
    """first argument is an RDCPot1 and the second is an IVM.
    this routine fixes all atoms such that the tensor does not vary."""
    ivm.fix( [convertToAtom(atom ,ivm.simulation) for atom in [p.oAtom(), p.xAtom(), p.yAtom(), p.zAtom() ,
                  p.o2Atom(), p.p1Atom(), p.p2Atom()]] )

def configIVM_fixAxis(p,ivm):
    """first argument is an RDCPot1 and the second is an IVM.
    this routine fixes all atoms such that the tensor orientation does
    not vary."""
    ivm.fix( [convertToAtom(atom ,ivm.simulation) for atom in [p.xAtom(), p.yAtom(), p.zAtom() ]] )

def configIVM_fixDa(p,ivm):
    """first argument is an RDCPot1 and the second is an IVM.
    this routine does topology setup for fixed Da."""
    if p.simulation().member().memberIndex()==0 or p.ensembleDa():
        ivm.group( [convertToAtom(atom ,ivm.simulation) for atom in [p.oAtom(), p.xAtom(), p.p1Atom()]] )
        ivm.breakBond( [convertToAtom(atom, ivm.simulation) for atom in [p.p1Atom(), p.oAtom()]] )
        pass
    return

def configIVM_fixRh(p,ivm):
    """first argument is an RDCPot1 and the second is an IVM.
    this routine does topology setup for fixed Da."""
    if p.simulation().member().memberIndex()==0 or p.ensembleRh():
        ivm.group( [convertToAtom(atom ,ivm.simulation) for atom in [p.xAtom(), p.p2Atom(), p.o2Atom()]] )
        ivm.breakBond( [convertToAtom(atom, ivm.simulation) for atom in [p.p2Atom(), p.o2Atom()]] )
        pass
    return


def configIVM_varyDa(p,ivm):
    """first argument is an RDCPot1 and the second is an IVM.
    this routine does topology setup for varying Da.

    the o-p1 bond will rotate about an axis perpendicular to the
    y-o-p1 atom plane. since only the projection on the x-o-y plane is
    significant in the calculation of Da, o-p1 should always have
    zero-projection along the o-z axis.
    """
    if p.simulation().member().memberIndex()==0 or p.ensembleDa():
        ivm.hinge( "bend", p.oAtom(), p.yAtom(), p.p1Atom() )
        ivm.group(  (p.oAtom().index(), p.p1Atom().index()) )
        pass
    return

def configIVM_varyRh(p,ivm):
    """first argument is a VarTensor and the second is an IVM.
    this routine does topology setup for varying Rhombicity.

    the o-p2 bond will rotate about an axis perpendicular to the
    z-o2-p2 atom plane. Note that only the azimuthal angle p2-o2-z is
    considered in the calculation of Rh.
    """
    if p.simulation().member().memberIndex()==0 or p.ensembleRh():
        ivm.hinge( "bend", p.o2Atom(), p.zAtom(), p.p2Atom() )
        ivm.group(  (p.o2Atom().index(), p.p2Atom().index()) )
        pass
    return

def configIVM_fixRhToOther(p,ivm,pot2):
    """first argument is a VarTensor, the second an IVM, and the third another
    VarTensor.
    this routine does topology setup for Rhombicity fixed to that of another
    potential term.
    NOTE: for this to work, the axes must also be fixed together. """
    if p.ensembleDa() != pot2.ensembleDa():
        raise Exception("mismatch in ensembleDa() between %s and %s" % \
                        (p.instanceName(), pot2.instanceName() ))
    rh = pot2.Rh()
    p.modified.set()
    p.simulation().markAsModified()
    if p.simulation().member().memberIndex()==0 or p.ensembleRh():
        p.setRh( rh )
        ivm.group( [p.xAtom(), p.o2Atom()] )
        ivm.breakBond( [p.p2Atom(), p.o2Atom()] )
        groups=ivm.groupList()
        inGroup=0
        for i in range(0,len(groups)):
            if pot2.p2Atom().index() in groups[i]:
                inGroup=1
                groups[i].append(p.p2Atom().index())
                pass
            pass
        ivm.setGroupList( groups )
        if not inGroup:
            ivm.group( (p.p2Atom(), pot2.p2Atom()) )
            pass
        pass
    return

def configIVM_fixAxisToOther(p,ivm,pot2):
    """first argument is a VarTensor, the second an IVM, and the third another
    VarTensor.
    this routine does topology setup to fix the axis (tensor orientation)
    to that of another potential term."""
    Da = p.Da() ; Rh = p.Rh() #preserve Da, Rh values

    p.modified.set()
    p.simulation().markAsModified()
    if p.simulation().singleThread():

        #    ivm.breakBond( (p.p2Atom().index(), p.o2Atom().index()) ) ??
        p.oAtom().setPos( pot2.oAtom().pos() )
        p.o2Atom().setPos( pot2.o2Atom().pos() )
        p.xAtom().setPos( pot2.xAtom().pos() )
        p.yAtom().setPos( pot2.yAtom().pos() )
        p.zAtom().setPos( pot2.zAtom().pos() )
        p.setDa( Da) ; p.setRh( Rh )
        groups=ivm.groupList()
        inGroup=0
        for i in range(0,len(groups)):
            if pot2.xAtom().index() in groups[i]:
                inGroup=1
                groups[i].append(p.xAtom().index())
                groups[i].append(p.yAtom().index())
                groups[i].append(p.zAtom().index())
                pass
            pass
        ivm.setGroupList( groups )
        if not inGroup:
            ivm.group( (p.xAtom(), p.yAtom(), p.zAtom(), pot2.xAtom()) )
            pass
        pass
    p.simulation().multiThread()
    return

def massSetup(list=[],axisMass=300,pAtomFactor=3):
    """
    appropriately setup tensor atom masses.
    tensor orientational atom's masses are set to axisMass, and
    the parameter atom masses are set to axisMass * pAtomFactor.

    if list is not specified, then all registered VarTensor objects are
    configured. 
    """
    try:
        len(list)
    except TypeError:
        list = [list]
        pass
    from simulationTools import getRegisteredTopoTerms
    if not list: list = getRegisteredTopoTerms("VarTensor")

    for t in list:
        for a in [t.oAtom,t.o2Atom,t.xAtom,t.yAtom,t.zAtom]:
            a().setMass(axisMass)
            pass
        for a in [t.p1Atom,t.p2Atom]:
            a().setMass(pAtomFactor*axisMass)
            pass
        pass
    return
 
    

def topologySetup(ivm,list=[]):
    """
    configure the given IVM object's topology setup using the freedom string
    for each VarTensor in list. This function should be called prior to
    <m ivm>.IVM.autoTorsion() or <m protocol>.torsionTopology()

    The freedom language contains the following keywords:
       fix
       fixAxis
       fixDa, varyDa
       fixRh, varyRh
       fixAxisTo <name>
       fixRhTo <name>
       ignore

    More than one can be specified, using a comma as a seperator.
    """
    try:
        len(list)
    except TypeError:
        list = [list]
        pass
    
    from simulationTools import getRegisteredTopoTerms
    if len(list)==0: list = getRegisteredTopoTerms("VarTensor",
                                                   ivm.simulation)

    for t in list:
        fixAll=0; fixAxis=0; fixAxisTo=""; fixDa=0; fixRh=0;
        statements=[x.strip() for x in t.freedom().split(',')]
        doIgnore=False
        for statement in statements:
            if statement=='fix':
                fixAll=1;
            elif statement=='fixAxis':
                fixAxis=1;
            elif statement=='fixDa':
                fixDa=1
            elif statement=='varyDa':
                fixDa=0
            elif statement=='fixRh':
                fixRh=1
            elif statement=='varyRh':
                fixRh=0
            elif statement.startswith('fixAxisTo'):
                fixAxisTo=statement.split()[-1]
            elif statement.startswith('fixRhTo'):
                fixRh=statement.split()[-1]
            elif statement=='ignore':
                doIgnore=True
            else:
                raise Exception("topologySetup: " +
                                "term: %s: " % t.instanceName() +
                                "did not understand statement: %s"%
                                statement)
            pass
        if doIgnore:
            continue
        
        configIVM(t,ivm)
        if fixAll:
            configIVM_fix(t,ivm)
            continue
        if fixAxis==1:
            configIVM_fixAxis(t,ivm)
        elif fixAxisTo!="":
            try:
                ot = [x for x in list if x.instanceName()==fixAxisTo][0]
            except:
                raise Exception("topologySetup: " +
                                "could not find tensor for fixAxisTo %s " %
                                fixAxisTo)
            configIVM_fixAxisToOther(t,ivm,ot)
            pass
            
        if fixDa==0:
            configIVM_varyDa(t,ivm)
        else:
            configIVM_fixDa(t,ivm)
            pass
        if fixRh==0:
            configIVM_varyRh(t,ivm)
        elif fixRh==1:
            configIVM_fixRh(t,ivm)
        else:
            try:
                ot = [x for x in list if x.instanceName()==fixRh][0]
            except:
                raise Exception("topologySetup: " +
                      "could not find tensor for fixRhTo %s " % fixRh)
            if fixAxisTo=="":
                raise Exception("topologySetup: " +
                      "fixAxisTo must be used in conjunction with "+
                      "fixRhTo")
            configIVM_fixRhToOther(t,ivm,ot)
            pass
        pass
    return

def syncAxisAtoms(t):
    """copy atoms representing orientation to member simulations>0"""
    esim = t.simulation()
    if esim.size()==1: return
    # atoms p1 and p2 ( Da and Rh) are not synced, but their positions
    # must be reset because the other axis atoms are moved (for member>0).
    da = t.Da()
    rh = t.Rh()
    from atomSel import AtomSel
    for atom in [t.oAtom(), t.o2Atom(),
                 t.xAtom(), t.yAtom(),t.zAtom()]:
        a = AtomSel('atom "%s" %d  %s' % (atom.segmentName(),
                                          atom.residueNum(),
                                          atom.atomName()))[0]
        a.setPos(atom.pos())
        pass
    esim.barrier()
    t.setDa( da )
    t.setRh( rh )
    esim.barrier()
    return
    
    
def alignTensorToZ(oTensor,
                   tensorTrans=0):
    """First compute the best-fit alignment tensor orientation, and then
    rotate the entire system such that the z-axis of the alignment tensor
    actually lies in the Cartesian z direction. The tensorTrans
    argument specifies the x and y components of an optional translation to be
    applied to all pseudoatoms used to describe the alignment tensors. The z
    component is always zero.
    """
    from varTensorTools import calcTensorOrientation, VarTensor_analyze
    from atomSel import AtomSel
    from atomSelAction import Rotate, Translate
    from mat3 import Mat3, transpose
    from vec3 import unitVec,Vec3
    calcTensorOrientation(oTensor)
    vz=unitVec(oTensor.zAtom().pos() - oTensor.oAtom().pos())
    vy=unitVec(oTensor.yAtom().pos() - oTensor.oAtom().pos())
    vx=unitVec(oTensor.xAtom().pos() - oTensor.oAtom().pos())
    AtomSel("all").apply( Rotate(transpose(Mat3(vx[0],vy[0],vz[0],
                                                vx[1],vy[1],vz[1],
                                                vx[2],vy[2],vz[2]))) )
    AtomSel("resname ANI*").apply(Translate(Vec3(tensorTrans,tensorTrans,0)))
    return

def copyTensor(ten1,ten2):
    """copy positions of tensor atoms from ten2 to ten1"""
    ten1.oAtom().setPos( ten2.oAtom().pos() )
    ten1.xAtom().setPos( ten2.xAtom().pos() )
    ten1.yAtom().setPos( ten2.yAtom().pos() )
    ten1.zAtom().setPos( ten2.zAtom().pos() )
    ten1.o2Atom().setPos( ten2.o2Atom().pos() )
    ten1.p1Atom().setPos( ten2.p1Atom().pos() )
    ten1.p2Atom().setPos( ten2.p2Atom().pos() )

    ten1.setDa( ten2.Da() )
    ten1.setRh( ten2.Rh() )
    return

def saupeMatrix(t,
                eIndex=0):
    """given a VarTensor, calculate the associated Saupe matrix.
    The eIndex argument is used for ensembles in which members are allowed
    different alignment tensors."""

    esim = t.simulation()
    if eIndex<0 or eIndex>=esim.size():
        raise Exception("eIndex is out of range: %d"%eIndex)

    oPos, xPos, yPos, zPos = [convertToAtom(atom,esim.members(eIndex)).pos() for atom in [t.oAtom(), t.xAtom(), t.yAtom(), t.zAtom()]]

    Da = t.Da(eIndex)
    Rh = t.Rh(eIndex)

    from vec3 import unitVec, norm
    from cdsMatrix import RMat, transpose
    Ul = []
    Ul.append( unitVec( xPos - oPos ))
    Ul.append( unitVec( yPos - oPos ))
    Ul.append( unitVec( zPos - oPos ))
    U = transpose( RMat(3,3).fromList(Ul) )

    Szz = 2*Da
    Sxx = (1.5*Rh-1)*Da
    Syy = -(Sxx + Szz)
    diag = RMat(3,3,0)
    diag.set(0,0, Sxx )
    diag.set(1,1, Syy )
    diag.set(2,2, Szz )
    tensor = U * diag * transpose(U)

    return tensor
    
    
def scalarProd(t1,t2):
    """Calculate the (unnormalized) scalar product between tensors t1 and t2.

    For internal use.

    """
    ret=0
    for i in range(3):
        for j in range(3):
            ret += t1.get(i,j) * t2.get(i,j)
            pass
        pass
    return ret
    

def normalizedScalarProduct(t1,t2):
    """Return the normalized scalar product between tensors t1 and t2.

    t1 and t2 are varTensor.VarTensor instances.  The normalized scalar product
    is defined by Sass et al. JACS (1999) 121:2047-2055.
    
    """
    from math import sqrt
    
    t1 = saupeMatrix(t1)
    t2 = saupeMatrix(t2)

    p11 = scalarProd(t1,t1)
    p22 = scalarProd(t2,t2)

    ret=0
    if (p11>1e-20) and (p22>1e-20):
        ret = scalarProd(t1,t2) / sqrt(p11*p22)
        pass

    return ret


def generalizedDegreeOfOrder(t):
    """Return the generalized degree of order of tensor t.

    t is a varTensor.VarTensor instance.  The generalized degree of order is
    defined as sqrt((2/3)*sum(Sij**2)), where sqrt is the square root and the
    sum extends over all ij elements of the associated Saupe order matrix (Sij).

    Reference: Tolman et al. JACS (2001) 123:1416.

    """
    t = saupeMatrix(t)
    return ((2.0/3.0) * scalarProd(t, t))**(0.5)


def orthogonalize(t):
    """given a VarTensor t, check that the axis atoms are orthgonal. If not,
    print a warning message and then orthonalize.
    """
    from vec3 import Vec3, dot, unitVec, norm
    
    qO = t.oAtom().pos()
    qx = t.xAtom().pos()
    qy = t.yAtom().pos()
    qz = t.zAtom().pos()

    vx = qx - qO;
    vy = qy - qO;
    vz = qz - qO;

    acc=1e-3
    if ( abs(dot(vx,vy))>acc or
         abs(dot(vx,vz))>acc or
         abs(dot(vy,vz))>acc   ) :
        print("Warning: axis atoms do not form orthogonal axes.")
        print("\tOrthogonalizing...")

        zp = unitVec(vz);
        yp = unitVec(vy - dot(vy,zp)*zp);
        xp = unitVec(vx - dot(vx,yp)*yp - dot(vx,zp)*zp);

        DaOK=False
        try:
            Da = t.Da() ; Rh = t.Rh() #preserve Da, Rh values
            DaOK=True
        except SystemError:
            pass
        t.xAtom().setPos( qO+norm(vx)*xp )
        t.yAtom().setPos( qO+norm(vy)*yp )
        t.zAtom().setPos( qO+norm(vz)*zp )
        if DaOK:
            t.setDa( Da) ; t.setRh( Rh )
            pass

        print("axes:               old                          new")
        print("\to: (%7.3f ,%7.3f ,%7.3f ) " %tuple(qO), end=' ')
        print("(%7.3f ,%7.3f ,%7.3f ) " % tuple(t.oAtom().pos()))
        print("\tx: (%7.3f ,%7.3f ,%7.3f ) "%tuple(qx), end=' ')
        print("(%7.3f ,%7.3f ,%7.3f ) " % tuple(t.xAtom().pos()))
        print("\ty: (%7.3f ,%7.3f ,%7.3f ) "%tuple(qy), end=' ')
        print("(%7.3f ,%7.3f ,%7.3f ) " % tuple(t.yAtom().pos()))
        print("\tz: (%7.3f ,%7.3f ,%7.3f ) "%tuple(qz), end=' ')
        print("(%7.3f ,%7.3f ,%7.3f ) "% tuple(t.zAtom().pos()))
        pass

    d=t.oAtom().pos()-t.o2Atom().pos()
    if norm(d)>acc:
        print("o2 atom position differs from that of o: resetting o2 position")
        DaOK=False
        try:
            Da = t.Da() ; Rh = t.Rh() #preserve Da, Rh values
            DaOK=True
        except SystemError:
            pass
        t.o2Atom().setPos(t.oAtom().pos() )
        if DaOK:
            t.setDa( Da) ; t.setRh( Rh )
            pass
        pass
    return


current_axisResid=5000
current_axisSegid="AXIS"
def addAxisAtoms(resid=-1,segid="",sim=None):
    """ add atoms to the current structure. If the resid argument
    is omitted, a new number is automatically generated. If the
    segment name is omitted,  the default is used.  Initial coords
    for the various atoms are also generated.

    Returns a tuple of the segid and residue number.
    """
    if not sim:
        from simulation import currentSimulation
        sim = currentSimulation()
        pass
    from xplorSimulation import getXplorSimulation
    xSim = getXplorSimulation(sim)
    if sim.type()=="SymSimulation":
        print("addAxisAtoms: attempting to add atoms to a SymSimulation.")
        print("   Not Supported.")
        #this might be fixed by setting xSim=sim at this point: to be tested.
        raise Exception("addAxisAtoms: attempting to add atoms to a "
                        "SymSimulation.")
    global current_axisResid, current_axisSegid
    if segid=="":
        segid = current_axisSegid
        pass
    if resid<0:
        while len( AtomSel("segid %s and resid %d" %
                           (segid,current_axisResid)))>0:
            current_axisResid +=1 
            pass
        resid = current_axisResid
        pass
    import protocol
    protocol.initParams("axis",silent=True)
    if len( AtomSel("segid %s and resid %s" % (segid,resid))) > 0:
        print("WARNING: addAxisAtoms: " + \
              "segid %s and resid %d already exist." % (segid,resid))
        return (segid,resid)
    cmd = psfTemplate.replace('_n__','%-4d'%resid)
    cmd2 = cmd.replace('SGMT','%-4s'%segid)
    if sim.type()=="EnsembleMemberSimulation":
        esim=ensembleSimulation.EnsembleSimulation_currentSimulation()
        from ensembleSimulation import EnsembleMemberSimulationPtr
        msim = EnsembleMemberSimulationPtr(sim.this)
        if esim.member().memberIndex()!=msim.memberIndex():
            print("addAxisAtoms: attempting to load atoms in", end=' ')
            print("nonlocal EnsembleSimulationMember.")
            cmd2=""
            pass
        pass
    outputState = xSim.disableOutput()
    xSim.fastCommand(cmd2)
    xSim.syncFrom()
    xSim.enableOutput(outputState)
    sim.sync()
    # set initial coords for the new axis atoms
    resSel = "segid %s and resid %d" % (segid,resid)
    from math import sqrt
    AtomSel(resSel + " and name OO" ,sim)[0].setPos( (0,0,0) )
    AtomSel(resSel + " and name X"  ,sim)[0].setPos( (1,0,0) )
    AtomSel(resSel + " and name Y"  ,sim)[0].setPos( (0,1,0) )
    AtomSel(resSel + " and name Z"  ,sim)[0].setPos( (0,0,1) )
    AtomSel(resSel + " and name OO2",sim)[0].setPos( (0,0,0) )
    # start with Da = 1/2 DaMax     
    AtomSel(resSel + " and name PA1",sim)[0].setPos( (1./sqrt(2),1./sqrt(2),0) )
    # start with Rh = 1/3           
    AtomSel(resSel + " and name PA2",sim)[0].setPos( (0,1./sqrt(2),1./sqrt(2)) )
    sim.sync()
    return (segid,resid)    

psfTemplate = """
structure
PSF

       3 !NTITLE
 REMARKS   auto-generated structure parameters
 REMARKS the string _n__ is replaced by the residue number
 RENARKS the string SGMT is replaced by the segment name

       7 !NATOM
       1 SGMT _n__ ANI  X    XXX    0.000000E+00   10.0000           0
       2 SGMT _n__ ANI  Y    YYY    0.000000E+00   10.0000           0
       3 SGMT _n__ ANI  Z    ZZZ    0.000000E+00   10.0000           0
       4 SGMT _n__ ANI  OO   OOO    0.000000E+00   10.0000           0
       5 SGMT _n__ ANI  OO2  OOO    0.000000E+00   10.0000           0
       6 SGMT _n__ ANI  PA1  PAR    0.000000E+00   10.0000           0
       7 SGMT _n__ ANI  PA2  PAR    0.000000E+00   10.0000           0

       6 !NBOND: bonds
       4       1       4       2       4       3       4       6
       5       3       5       7

       3 !NTHETA: angles
       1       4       2       1       4       3       2       4       3

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

#
# given a known structure, and experimental RDCs, determine the dipolar
# coupling order matrix, Da, and rhombicity
#
def calcTensor(vTensor,expts=0,
               coords=None,
               weights=None,
               useErr=False,
               svdTolerance=0.1,
               svdNum=5,
               selection="all",
               suppressExceptions=False,
               ):
    """Given a known structure, and a <m varTensor>.VarTensor object with
    associated experimental RDC and/or CSA terms, determine the orientation
    tensor, Da, and rhombicity.  The Saupe matrix is returned (in the axis
    system of the input coordinates).  If there is an ensemble of
    structures, a single orientation tensor is calculated.  Note that when the
    varTensor instance (vTensor argument) has multiple associated experiments,
    each term in the solution matrix is premultiplied by 1/sqrt( term.scale() ),
    so that the appropriate scaling is maintained.

    The optional expts argument specifies which RDC and/or CSA experiments to
    use in the calculation of the tensor; if omitted, all experiments associated
    with the tensor will be used in the tensor calculation.

    If the coords argument is specified, it should be a sequence of sets of
    atomic coordinates (e.g., as obtained from <m simulation>.atomPosArr) to be 
    used in determination of the tensor.  If the sequence has more than one
    coordinate set, the weights argument specifies the relative importance of
    each set (uniform weight if weights is not specified).

    Set the useErr argument to True to appropriately weight the SVD calculation
    to take account of observed RDC and CSA errors.  This defaults to False.

    The selection argument can be used to select a subset of restraints to use
    for the SVD calculation.  The restraints used are those whose atom
    selections lie within the selection argument.

    svdTolerance specifies the size of the smallest singular value used in
    solving the linear equation.  The equation for the unique elements of the
    Saupe tensor is

        t = v * diag * uT * b 

    where v, diag and uT are results of SVD and b is a vector of observed
    measurements.  diag is a diagonal matrix, the nonzero elements of which are
    reciprocals of singular values of a matrix composed of geometrical bond
    vector information.  If the absolute value of a singular value is less than
    svdTolerance times the absolute value of the next largest singular value,
    the corresponding value of diag and those for all smaller singular values
    are set to zero.

    Alternatively, svdNum specifies the number of singular values to include
    (of a maximum of 5).

    The suppressExceptions argument can be set to True to avoid an exception
    being thrown if there are fewer than 5 observed RDCs. In that case, the
    alignment tensor is underdetermined, and the SVD solution will not be
    unique. 

    Reference: Losonczi, J.A., et al. J. Magn. Reson. (1999) 138: 334-342.

    """
    from vec3 import unitVec, norm, Vec3
    from rdcPotTools import DaScale
    from mat3 import SymMat3, eigen
    from cdsMatrix import RMat, transpose, svd
    from math import sqrt

    from selectTools import convertToAtomSel
    if not expts: expts = vTensor.expts()

    rdcs = [x for x in expts if x.potName()=='RDCPot1' and x.scale()>1e-15]
    csas = [x for x in expts if x.potName()=='CSAPot'  and x.scale()>1e-15]

    # first call to restraints() in multi-threaded region
    list(map(lambda pot: pot.rawRestraints(), rdcs+csas))

    if not rdcs and not csas: return

    if rdcs:
        ens = vTensor.simulation()
    else:
        ens = vTensor.simulation()
        pass
    
    Da = ens.sharedObj(0.)
    R = ens.sharedObj(0.)
    xPos= ens.sharedObj((0,0,0))
    yPos= ens.sharedObj((0,0,0))
    zPos= ens.sharedObj((0,0,0))
    savedX = ens.sharedObj(1)

    exptSim = rdcs[0].simulation() if rdcs else csas[0].simulation()
    selection = convertToAtomSel(selection,exptSim)

    if ens.singleThread():

        coordsSets=[]
        if coords:
            if not weights: weights = [1./len(coords)]*len(coords)
            for i in range(len(coords)):
                coordsSets.append( (coords[i],weights[i]) )
                pass                           
            pass
        else:
            for cnt in range(exptSim.size()):
                mem=exptSim.members(cnt)
                coordsSets.append( (mem.atomPosArr(), mem.weight()) )
                pass
            pass

        nRestraints=0
        for rdc in rdcs:
            for restraint in rdc.rawRestraints():
                include=False
                for selPair in restraint.selPairs():
                    for ai,aAtom in enumerate(selPair.a):
                        if selection.containsAtom(aAtom):
                            for bi in range(len(selPair.b)):
                                if rdc.aveType()=="pairwise":
                                    bi=ai
                                    pass
                                bAtom = selPair.b[bi]
                                if selection.containsAtom(bAtom):
                                    include=True
                                    break
                                if rdc.aveType()=="pairwise":
                                    break
                                pass
                            pass
                        pass
                    if include:
                        break
                    pass
                if include:
                    nRestraints+=1
                pass
            pass

        sim = exptSim.subSim()
        saveCoords = sim.atomPosArr()

        for csa in csas:
            for restraint in csa.rawRestraints():
                if selection.containsAtom(restraint.Selection1()[0]):
                    nRestraints += 1
                    pass
                pass
            pass

        if nRestraints<5:
            mess="insufficient restraints found: %d (<5)" % nRestraints
            if suppressExceptions:
                print("calcTensor: Warning:", mess)
                return
            else:
                raise Exception(mess)
            pass

        m = RMat( nRestraints , 5)
        b = CDSVector_double( nRestraints )

        row=0
        for rdc in rdcs:
            for restraint in rdc.rawRestraints():
                scale = sqrt( rdc.scale() )
                m[row,0] = m[row,1] = m[row,2] = m[row,3] = m[row,4] =0
                skip=False
                for selPair in restraint.selPairs():
                    for ai,aAtom in enumerate(selPair.a):
                        if not selection.containsAtom(aAtom):
                            skip=True
                            continue
                        for bi in range(len(selPair.b)):
                            if rdc.aveType()=="pairwise":
                                bi=ai
                                pass
                            bAtom = selPair.b[bi]
                            if not selection.containsAtom(bAtom):
                                skip=True
                                continue
                            for (coords, weight) in coordsSets:
                                aPos = coords[aAtom.index()]
                                bPos = coords[bAtom.index()]
                                
                                uVec = unitVec( aPos - bPos )

                                const = scale * weight * \
                                        DaScale(rdc,norm(aPos - bPos))
                            
                                if rdc.aveType()!="sum":
                                    const /= float(restraint.aveSize())
                                    #warning: aveSize is not correct if selection
                                    # selects a subset of selPairs...
                                    pass

                                if useErr:
                                    const /= 0.5*(restraint.plusErr()+
                                                  restraint.minusErr())
                                    pass
                            
                                m[row,0] += const * (uVec[1]**2 - uVec[0]**2)
                                m[row,1] += const * (uVec[2]**2 - uVec[0]**2)
                                m[row,2] += const * (2 * uVec[0] * uVec[1]  )
                                m[row,3] += const * (2 * uVec[0] * uVec[2]  )
                                m[row,4] += const * (2 * uVec[1] * uVec[2]  )
                                pass
                            if rdc.aveType()=="pairwise":
                                break
                            pass
                        pass
                    pass
                if not skip:
                    b[row] = scale * restraint.obs()
                    if useErr:
                        b[row] /= 0.5*(restraint.plusErr()+
                                       restraint.minusErr())
                        pass
                    row+=1
                    pass
                if row>=nRestraints: break
                pass
            pass

        for csa in csas:
            scale = sqrt( csa.scale() )
            for restraint in csa.rawRestraints():
                if selection.containsAtom(restraint.Selection1()[0]):
                    m[row,0] = m[row,1] = m[row,2] = m[row,3] = m[row,4] =0
                    for (coords, weight) in coordsSets:
                        sim.setAtomPosArr(coords)

                        delta = restraint.tensor()

                        const = scale * weight / csa.DaScale()
                        
                        if useErr:
                            const /= 0.5*(restraint.plusErr()+
                                       restraint.minusErr())
                            pass


                        m[row,0] += const * (-delta[0,0] + delta[1,1])
                        m[row,1] += const * (-delta[0,0] + delta[2,2])
                        m[row,2] += const * ( delta[0,1] + delta[1,0])
                        m[row,3] += const * ( delta[0,2] + delta[2,0])
                        m[row,4] += const * ( delta[1,2] + delta[2,1])

                        pass
                    b[row] = 3./2 * scale * restraint.obs() 
                    if useErr:
                        b[row] /= 0.5*(restraint.plusErr()+
                                       restraint.minusErr())
                        pass

                    row+=1
                    if row>=nRestraints: break
                    pass
                pass
            pass

        sim.setAtomPosArr( saveCoords )

        svdResults = svd(m,'S')

        v = transpose( svdResults.vT )
        uT = transpose( svdResults.u )
        sigma = svdResults.sigma 

        diag = RMat(5,5,0)
        for i in range(svdNum): # sigma sorted in decreasing order
            diag[i,i] = 1.0/sigma[i]
            if i<4 and (sigma[i+1]<1e-16 or
                        abs(sigma[i]/sigma[i+1])> 1./svdTolerance):
                print('calcTensor: large sigma range found:')
                print('  sigma: ',['%.3f'%s for s in sigma])
                print('  using elements 0..%d' % i)
                break
            pass

        # x is [ Syy Szz Sxy Sxz Syz ]
        x = v * diag * uT * b
        savedX.set( list(x) )

        S =SymMat3()
        S.set(0,0, -x[0]-x[1]); 
        S.set(1,0,x[2]);        S.set(1,1,x[0]) 
        S.set(2,0,x[3]);        S.set(2,1,x[4]);  S.set(2,2,x[1])

        eigenPairs = eigen(S)

        xInd = 0; yInd=1; zInd=2;
        
        if abs(eigenPairs[zInd].value()) < abs(eigenPairs[yInd].value()):
            zInd,yInd = yInd,zInd
            pass
        if abs(eigenPairs[zInd].value()) < abs(eigenPairs[xInd].value()):
            zInd,xInd = xInd,zInd
            pass
        if abs(eigenPairs[yInd].value()) < abs(eigenPairs[xInd].value()):
            yInd,xInd = xInd,yInd
            pass

        Azz = 2.0/3.0 * eigenPairs[zInd].value()  
        Ayy = 2.0/3.0 * eigenPairs[yInd].value()  
        Axx = 2.0/3.0 * eigenPairs[xInd].value()

        Aa = 3 * Azz / 2.
        Ar = (Axx-Ayy)
        
        Da.set( 0.5 * Aa )

        R.set( Ar/Aa )
        xPos.set( eigenPairs[xInd].vector() )
        yPos.set( eigenPairs[yInd].vector() )
        zPos.set( eigenPairs[zInd].vector() )
        
        pass

    ens.multiThread()

    if vTensor.oAtom().isValid():
        oPos = vTensor.oAtom().pos()
    else:
        oPos = Vec3(0,0,0)
        vTensor.oAtom().setPos( oPos )
        pass

    xPos = Vec3(xPos.get())
    yPos = Vec3(yPos.get())
    zPos = Vec3(zPos.get())

    from vec3 import dot, cross
    zPos *= dot(cross(xPos,yPos) , zPos) # make right-handed
        
    vTensor.o2Atom().setPos( oPos )
    vTensor.xAtom().setPos( xPos+oPos )
    vTensor.yAtom().setPos( yPos+oPos )
    vTensor.zAtom().setPos( zPos+oPos )
    vTensor.setDa( Da.get() )
    vTensor.setRh( R.get() )

    #update pot values
    for rdc in rdcs: rdc.calcEnergy()
    for csa in csas: csa.calcEnergy()

    x = savedX.get()

    S =SymMat3()
    S.set(0,0, -x[0]-x[1]); 
    S.set(1,0,x[2]);        S.set(1,1,x[0]) 
    S.set(2,0,x[3]);        S.set(2,1,x[4]);  S.set(2,2,x[1])
    
    return S


def calcTensor_ensemble(vTensor,expts=0,
                        svdTolerance=0.01,
                        useErr=False,
                        selection='known'):
    """given a known structure, and a VarTensor object with associated
    experimental RDC or CSA terms, determine the orientation tensor, Da, and
    rhombicity. The Saupe matrix is returned. If there is an ensemble of
    structures, multiple alignment tensors are calculated. Note when vTensor
    has multiple associated experiments: each term in the solution matrix is
    premultiplied by 1/sqrt( term.scale() ), so that the appropriate scaling
    is maintained.

    The optional expts argument specifies which RDC or CSA experiments to use
    in calculation of the tensor. If it is omitted, all expts associated with
    the tensor will be used in the tensor calculation.

    The selection argument can be used to select a subset of restraints to
    use for the SVD calculation. The restraints used are those whose atom
    selections lie within the selection argument.

    Set the useErr argument to True to appropriately weight the
    SVD calculation to take account of observed RDC and CSA errors. This
    defaults to False.

    svdTolerance specifies the size of the smallest singular value used in
    solving the linear equation. The equation for the unique elements of the
    Saupe tensor is
       t = v * diag * uT * b
    where v, diag and uT are results of SVD and b is a vector of observed
    measurements. diag is a diagonal matrix, the nonzero elements of which
    are reciprocals of singular values of a matrix composed of geometrical
    bond vector information. The elements of diag are set to zero if the
    absolute value of the corresponding singular value is less than
    svdTolerance times the average of the absoulte value of singular values.
    """
    from vec3 import unitVec, norm, Vec3
    from rdcPotTools import DaScale
    from mat3 import SymMat3, eigen
    from cdsMatrix import RMat, transpose, svd
    from math import sqrt

    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)
    if not expts: expts = vTensor.expts()

    rdcs = [x for x in expts if x.potName()=='RDCPot1']
    csas = [x for x in expts if x.potName()=='CSAPot']

    # first call to restraints() in multi-threaded region
    list(map(lambda pot: pot.rawRestraints(), rdcs+csas))

    if not rdcs and not csas: return

    if rdcs:
        ens = rdcs[0].simulation()
    else:
        ens = csas[0].simulation()
        pass

    eIndex = ens.member().memberIndex()
    msim = ens.member()
    
    Da = []; R  = []; xPos= []; yPos= []; zPos= []; savedX = []

    for cnt in range(ens.size()):
        Da.append(    ens.sharedObj(0.)      )
        R.append(     ens.sharedObj(0.)      )
        xPos.append(  ens.sharedObj((0,0,0)) )
        yPos.append(  ens.sharedObj((0,0,0)) )
        zPos.append(  ens.sharedObj((0,0,0)) )
        savedX.append(ens.sharedObj(1)       )
        pass

    

    mlT = [[] for i in range(5)]
    b = []

    for rdc in rdcs:
        for restraint in rdc.rawRestraints():
            include=False
            scale = sqrt( rdc.scale() )
            first=True
            for selPair in restraint.selPairs():
                for ai,aAtom in enumerate(selPair.a):
                    if not selection.containsAtom(aAtom):
                        continue
                    for bi in range(len(selPair.b)):
                        if rdc.aveType()=="pairwise":
                            bi=ai
                            pass
                        bAtom = selPair.b[bi]
                        if not selection.containsAtom(bAtom):
                            continue

                        include=True
                        if first:
                            first=False
                            for col in range(5):
                                mlT[col].append(0)
                                pass
                            pass
                            
                        aPos = aAtom.pos()
                        bPos = bAtom.pos()
                        
                        uVec = unitVec( aPos - bPos )
                                
                        const = scale * rdc.ensWeights()[eIndex] * \
                                DaScale(rdc,norm(aPos - bPos))
                                
                        if rdc.aveType()!="sum":
                            const /= float(restraint.aveSize())
                            pass

                        if useErr:
                            const /= 0.5*(restraint.plusErr()+
                                       restraint.minusErr())
                            pass

                        mlT[0][-1] +=const*(uVec[1]**2 - uVec[0]**2)
                        mlT[1][-1] +=const*(uVec[2]**2 - uVec[0]**2)
                        mlT[2][-1] +=const*(2 * uVec[0] * uVec[1]  )
                        mlT[3][-1] +=const*(2 * uVec[0] * uVec[2]  )
                        mlT[4][-1] +=const*(2 * uVec[1] * uVec[2]  )
                        
                        if rdc.aveType()=="pairwise":
                            break
                        pass
                    pass
                pass

            if include:
                b.append( scale * restraint.obs() )
                if useErr:
                    b[-1] /= 0.5*(restraint.plusErr()+
                                       restraint.minusErr())
                    pass
                pass

            pass
        pass

    for csa in csas:
        scale = sqrt( csa.scale() )
        for restraint in csa.rawRestraints():

            if not selection.containsAtom(restraint.Selection1()[0]):
                continue

                
            delta = restraint.tensor()

            const = scale * csa.ensWeights()[eIndex] / csa.DaScale()

            if useErr:
                const /= 0.5*(restraint.plusErr()+
                              restraint.minusErr())
                pass

            mlT[0].append( const * (-delta[0,0] + delta[1,1]) )
            mlT[1].append( const * (-delta[0,0] + delta[2,2]) )
            mlT[2].append( const * ( delta[0,1] + delta[1,0]) )
            mlT[3].append( const * ( delta[0,2] + delta[2,0]) )
            mlT[4].append( const * ( delta[1,2] + delta[2,1]) )
                
            b.append( 3./2 * scale * restraint.obs() )
            if useErr:
                b[-1] /= 0.5*(restraint.plusErr()+
                              restraint.minusErr())
                pass

            pass
        pass

    b = CDSVector_double(b)

    nRestraints = len(b)

    e_mlT=[]
    for cnt in range(ens.size()):
        e_mlT.append( ens.sharedObj( [list(l) for l in mlT] ) )
        pass

    e_mlT[eIndex].set( [list(l) for l in mlT] )

    if ens.singleThread():

        m = RMat( nRestraints , 5*ens.size())

        for cnt in range(ens.size()):
            mlT = e_mlT[cnt].get()
            for i in range(nRestraints):
                for j in range(5):
                    m[i,j+cnt*5] = mlT[j][i]
                    pass
                pass
            pass


        numParams = 5*ens.size()
        if len(b)<numParams:
            raise Exception("insufficient restraints found: %d" % len(b))


        svdResults = svd(m,'S')
    
        v = transpose( svdResults.vT )
        uT = transpose( svdResults.u )
        diag = RMat(numParams,numParams,0)
        aveVal=0
        for i in range(numParams):
            aveVal += abs(svdResults.sigma[i])
            pass
        aveVal /= numParams
        
        for i in range(numParams):
            if abs(svdResults.sigma[i]/aveVal)>svdTolerance:
                diag[i,i] = 1.0/svdResults.sigma[i]
            else:
                print('calcTensor_ensemble: discarding singular value:', \
                      svdResults.sigma[i])
            pass
    
        # x is [ Syy Szz Sxy Sxz Syz ]
        fullX = v * diag * uT * b

        for cnt in range(ens.size()):
            x = fullX[cnt*5:(cnt+1)*5]
            savedX[cnt].set( list(x) )
    
            S =SymMat3()
            S.set(0,0, -x[0]-x[1]); 
            S.set(1,0,x[2]);        S.set(1,1,x[0]) 
            S.set(2,0,x[3]);        S.set(2,1,x[4]);  S.set(2,2,x[1])
    
        

            eigenPairs = eigen(S)
    
            xInd = 0; yInd=1; zInd=2;
            
            if abs(eigenPairs[zInd].value()) < abs(eigenPairs[yInd].value()):
                zInd,yInd = yInd,zInd
                pass
            if abs(eigenPairs[zInd].value()) < abs(eigenPairs[xInd].value()):
                zInd,xInd = xInd,zInd
                pass
            if abs(eigenPairs[yInd].value()) < abs(eigenPairs[xInd].value()):
                yInd,xInd = xInd,yInd
                pass
    
            Azz = 2.0/3.0 * eigenPairs[zInd].value()
            Ayy = 2.0/3.0 * eigenPairs[yInd].value()
            Axx = 2.0/3.0 * eigenPairs[xInd].value()
    
            Aa = 3 * Azz / 2.
            Ar = (Axx-Ayy)

            Da[cnt].set( 0.5*Aa )
            

            R[cnt].set( Ar/Aa )
            xPos[cnt].set( eigenPairs[xInd].vector() )
            yPos[cnt].set( eigenPairs[yInd].vector() )
            zPos[cnt].set( eigenPairs[zInd].vector() )
            pass
        
        pass
    
    ens.multiThread()

    vTensor.setEnsembleAxis(1)
    vTensor.setEnsembleDa(1)
    vTensor.setEnsembleRh(1)

    oPos = vTensor.oAtom().pos()
    vTensor.o2Atom().setPos( oPos )

    xPos = Vec3(xPos[eIndex].get())
    yPos = Vec3(yPos[eIndex].get())
    zPos = Vec3(zPos[eIndex].get())

    from vec3 import dot, cross
    zPos *= dot(cross(xPos,yPos) , zPos) # make right-handed
         
    vTensor.xAtom().setPos( xPos+oPos )
    vTensor.yAtom().setPos( yPos+oPos )
    vTensor.zAtom().setPos( zPos+oPos )

    vTensor.setDa( Da[eIndex].get() )
    vTensor.setRh( R[eIndex].get() )

    #update pot values
    for rdc in rdcs: rdc.calcEnergy()
    for csa in csas: csa.calcEnergy()

    Sret=[]
    for cnt in range(ens.size()):
        x = savedX[cnt].get()

        S =SymMat3()
        S.set(0,0, -x[0]-x[1]); 
        S.set(1,0,x[2]);        S.set(1,1,x[0]) 
        S.set(2,0,x[3]);        S.set(2,1,x[4]);  S.set(2,2,x[1])

        Sret.append(S)
        pass
        
    return Sret


def calcTensorOrientation(oTensor):
    """
    given a known structure, and experimental RDCs, and 
    Da and rhombicity, determine the tensor orientation.

    All rdcs in oTensor.expts are used in determining the tensor.

    The Saupe matrix is returned.
    """
    
    from ivm import IVM
    ivm = IVM()
    for rdc in oTensor.expts():
        ivm.potList().add( rdc )
        pass

    ivm.fix(' not (segid "%s" and resid %d)' %
            (oTensor.oAtom().segmentName(),oTensor.oAtom().residueNum()))
    oldFreedom = oTensor.freedom()
    oTensor.setFreedom("fixDa, fixRh")
    topologySetup(ivm,oTensor)

    import protocol
    protocol.initMinimize(ivm)
    ivm.setDEpred( 10 )
    ivm.setNumSteps(500)
    ivm.setPrintInterval(10)
    ivm.run()

    oTensor.setFreedom(oldFreedom)
    return saupeMatrix(oTensor)

def getVarTensors(potList):
    """ given a list of potential terms, return a list of unique VarTensor
    objects either in the list or refered to by another potential type.
    """
    uniqueList=[]
    for pot in potList:
        try:
            if hasattr(pot,"oTensor"): pot = pot.oTensor
        except:
            pass
        if pot.potName()=="VarTensor":
            if not pot.instanceName() in [p.instanceName() for p in uniqueList]:
                uniqueList.append( pot )
                pass
            pass
        pass
    return uniqueList
            
                     
def VarTensor_analyze(potList):
    "perform analysis of VarTensor terms and return nicely formatted summary"
    from math import atan2, acos, pi, sqrt, sin
    from vec3 import unitVec, norm

    ret = ""

    from simulationTools import getPotTerms
    varTensors = getPotTerms(potList,'VarTensor')

    rdccsa = getPotTerms(potList,('RDCPot1','CSAPot'))

    varTensors += [x.oTensor for x in rdccsa]

    if not varTensors: return ''

    varTensors.sort(key=lambda x:x.instanceName())
    rVarTensors = varTensors
    varTensors=[]
    for t in rVarTensors:
        if not t.instanceName() in [x.instanceName() for x in varTensors]:
            varTensors.append(t)
            pass
        pass

    for t in varTensors: 
        if not t.ensembleAxis():
            syncAxisAtoms( t )
            pass
        pass

    esim = varTensors[0].simulation()
    if esim.size()>1:
        ret += ' %-9s  %7s %7s %7s %7s %7s %7s\n' %\
                ('','aveDa', 'Da', 'rms(Da)', 'aveRh', 'Rh', 'rms(Rh)')
        for t in varTensors:
            rmsDa = 0
            rmsRh = 0
            for cnt in range(esim.size()):
                rmsDa += esim.weight(cnt) * (t.Da(cnt) - t.aveDa())**2
                rmsRh += esim.weight(cnt) * (t.Rh(cnt) - t.aveRh())**2
                pass
            rmsDa = sqrt(rmsDa)
            rmsRh = sqrt(rmsRh)
            ret += ' %-9s  %7.3f %7.3f %7.3f %7.3f %7.3f %7.3f\n' %\
                   (t.instanceName(),
                    t.aveDa(),t.Da(),rmsDa,
                    t.aveRh(),t.Rh(),rmsRh)
            pass
        pass
    
    # print out info on each tensor
    for t in varTensors:
        if t.ensembleAxis():
            ret += ensembleTensorInfo(t)
            continue

        ret += "\nVarTensor: %s" % t.instanceName()
        ret += "      [resid %d]\n" % t.oAtom().residueNum()
        ret += "  Da: %6.2f  Rhombicity: %6.3f  (GDO: %6.2f)\n" % (t.Da(),
                                                                   t.Rh(),
                                                generalizedDegreeOfOrder(t))
        z = unitVec(t.zAtom().pos() - t.oAtom().pos())
        y = unitVec(t.yAtom().pos() - t.oAtom().pos())
        x = unitVec(t.xAtom().pos() - t.oAtom().pos())
        ret += "  principal axes: z: %6.3f %6.3f %6.3f\n" % tuple(z)
               
        ret += "                  y: %6.3f %6.3f %6.3f\n" % tuple(y)
        ret += "                  x: %6.3f %6.3f %6.3f\n" % tuple(x)
        S = saupeMatrix(t)
        ret += '  Saupe matrix:     %8.3f %8.3f %8.3f\n' %\
               ( S.get(0,0) ,S.get(0,1), S.get(0,2) )
        ret += '                    %8.3f %8.3f %8.3f\n' %\
               ( S.get(1,0) ,S.get(1,1), S.get(1,2) )
        ret += '                    %8.3f %8.3f %8.3f\n' %\
               ( S.get(2,0) ,S.get(2,1), S.get(2,2) )
        #
        # Euler angles computed using Goldstein's x-convention
        # p. 147
        #
        #
        # Euler Angles are for rotation matrix which converts the
        #  matrix R whose columns are composed of the vectors [ijk]
        #  into a matrix whose columns comprise the current [xyz]
        #  vectors:
        #           [xyz] = R [ijk]
        #
        # in this def. R is simply [xyz]^T
        #
        ret += '\n  Euler Angles for 4 equiv. orientations\n'
        line1 = "   phi    :"
        line2 = "   theta  :"
        line3 = "   psi    :"
        z0,z1,z2 = tuple(z)
        x2 = x[2]
        y2 = y[2]
        def sign(x):
            if x<0: return -1
            return 1
        for (s1,s2,s3) in ((1,1,1),
                           #(-1,1,1),     #these involve inversion
                           #(1,-1,1),
                           #(1,1,-1),
                           (1,-1,-1),
                           (-1,1,-1),
                           (-1,-1,1),
                           #(-1,-1,-1),
                           ):
            if (x2==0. and y2==0.) or (z0==0. and z1==0.):
                #FIX: should properly handle case of being exactly on axis
                phi,theta,psi = (0,0,0)
            else:
                psi,theta,phi = (atan2(s1*x2,s2*y2),
                                 acos(s3*z2)           ,
                                 atan2(s3*z0,-z1*s3)   )
                theta *= sign( sin(psi) * z2 * s3)
                psi = atan2(s1*x2*sign(theta),s2*y2*sign(theta))
                phi = atan2(s3*z0*sign(theta),-z1*s3*sign(theta)) 
                pass
            # get in range 0..2*PI
            phi,theta,psi = [x-2*pi*(-1+(x>0)) for x in (phi,theta,psi)]
            # convert to degrees
            phi,theta,psi = [x*180./pi for x in (phi,theta,psi)]
            line1 += " %7.2f " % phi
            line2 += " %7.2f " % theta
            line3 += " %7.2f " % psi
            pass
        ret += line1 + '\n' + line2 + '\n' + line3 + '\n'

        pass

    if len(varTensors)<2:
        return ret

        
    ret += '-'*60 + '\n'
    ret += "  Normalized scalar product between tensors \n\n"
    ret += "%8s "%' '
    for t1 in varTensors:
        ret += "%8s "%t1.instanceName()
        pass
    ret += '\n'
    for t1 in varTensors:
        ret += "%8s "%t1.instanceName()
        for t2 in varTensors:
            ret += "%8.3f " % normalizedScalarProduct(t1,t2)
            pass
        ret += '\n'
        pass

    return ret

def eulerAngles(t,eIndex=-1):
    """
    For the given <m varTensor>.VarTensor object, return the four sets of
    equivalent euler angles.

    
     Euler angles computed using Goldstein's x-convention
     p. 147
    
    
     Euler Angles are for rotation matrix which converts the
      matrix R whose columns are composed of the vectors [ijk]
      into a matrix whose columns comprise the current [xyz]
      vectors:
               [xyz] = R [ijk]
    
     in this def. R is simply [xyz]^T
    
    """
    esim=t.simulation()
    if eIndex<0:
        eIndex=esim.member().memberIndex()
        pass
    
    
    oPos, xPos, yPos, zPos = [AtomSel('atom "%s" %d  %s' %
                                         (atom.segmentName(),
                                          atom.residueNum(),
                                          atom.atomName()),
                                         esim.members(eIndex))[0].pos() for atom in (t.oAtom(), t.xAtom(),
                                  t.yAtom() , t.zAtom())]

    from vec3 import unitVec
    from math import sin, acos, atan2, pi
    
    z = unitVec(zPos - oPos)
    y = unitVec(yPos - oPos)
    x = unitVec(xPos - oPos)
    S = saupeMatrix(t,eIndex)

    z0,z1,z2 = tuple(z)
    x2 = x[2]
    y2 = y[2]
    def sign(x):
        if x<0: return -1
        return 1
    ret=[]
    for (s1,s2,s3) in ((1,1,1),
                       #(-1,1,1),     #these involve inversion
                       #(1,-1,1),
                       #(1,1,-1),
                       (1,-1,-1),
                       (-1,1,-1),
                       (-1,-1,1),
                       #(-1,-1,-1),
                       ):
        if (x2==0. and y2==0.) or (z0==0. and z1==0.):
            #FIX: should properly handle case of being exactly on axis
            phi,theta,psi = (0,0,0)
        else:
            psi,theta,phi = (atan2(s1*x2,s2*y2),
                             acos(s3*z2)           ,
                             atan2(s3*z0,-z1*s3)   )
            theta *= sign( sin(psi) * z2 * s3)
            psi = atan2(s1*x2*sign(theta),s2*y2*sign(theta))
            phi = atan2(s3*z0*sign(theta),-z1*s3*sign(theta)) 
            pass
        # get in range 0..2*PI
        phi,theta,psi = [x-2*pi*(-1+(x>0)) for x in (phi,theta,psi)]
        # convert to degrees
        phi,theta,psi = [x*180./pi for x in (phi,theta,psi)]
        ret.append( (phi,theta,psi) )
        pass
    return ret
    

def ensembleTensorInfo(t):
    """ get info for all tensors of the ensemble
    """
    from math import atan2, acos, pi, sqrt, sin
    from vec3 import unitVec, norm

    ret=''
    ret += "\nVarTensor: %s" % t.instanceName()
    ret += "      [resid %d]\n" % t.oAtom().residueNum()

    esim = t.simulation()

    for eIndex in range(esim.size()):
        ret += " [Tensor for ensemble member %d]\n" % eIndex
        ret += "  Da: %6.2f  Rhombicity: %6.3f\n" % (t.Da(eIndex),
                                                     t.Rh(eIndex))
        oPos, xPos, yPos, zPos = [AtomSel('atom "%s" %d  %s' %
                                             (atom.segmentName(),
                                              atom.residueNum(),
                                              atom.atomName()),
                                             esim.members(eIndex))[0].pos() for atom in (t.oAtom(), t.xAtom(),
                                              t.yAtom() , t.zAtom())]

        z = unitVec(zPos - oPos)
        y = unitVec(yPos - oPos)
        x = unitVec(xPos - oPos)
        ret += "  principal axes: z: %6.3f %6.3f %6.3f\n" % tuple(z)
               
        ret += "                  y: %6.3f %6.3f %6.3f\n" % tuple(y)
        ret += "                  x: %6.3f %6.3f %6.3f\n" % tuple(x)
        S = saupeMatrix(t,eIndex)
        ret += '  raw matrix form:  %8.3f %8.3f %8.3f\n' %\
               ( S.get(0,0) ,S.get(0,1), S.get(0,2) )
        ret += '                    %8.3f %8.3f %8.3f\n' %\
               ( S.get(1,0) ,S.get(1,1), S.get(1,2) )
        ret += '                    %8.3f %8.3f %8.3f\n' %\
               ( S.get(2,0) ,S.get(2,1), S.get(2,2) )
        #
        # Euler angles computed using Goldstein's x-convention
        # p. 147
        #
        #
        # Euler Angles are for rotation matrix which converts the
        #  matrix R whose columns are composed of the vectors [ijk]
        #  into a matrix whose columns comprise the current [xyz]
        #  vectors:
        #           [xyz] = R [ijk]
        #
        # in this def. R is simply [xyz]^T
        #
        ret += '\n  Euler Angles for 4 equiv. orientations\n'
        line1 = "   phi    :"
        line2 = "   theta  :"
        line3 = "   psi    :"
        for (phi,theta,psi) in eulerAngles(t,eIndex):
            line1 += " %7.2f " % phi
            line2 += " %7.2f " % theta
            line3 += " %7.2f " % psi
            pass
        ret += line1 + '\n' + line2 + '\n' + line3 + '\n'
        ret += '\n'
        pass
    return ret



def composite_RDCRfactor(term):
    """For the RDC terms associated with the argument VarTensor object,
    return composite R-factor by calling
    <m rdcPotTools>.composite_Rfactor_infinite."""
    from simulationTools import getPotTerms
    rdcTerms = getPotTerms(term.expts(),'RDCPot1')
    
    from rdcPotTools import composite_Rfactor_infinite
    return composite_Rfactor_infinite(rdcTerms)

def composite_RDCchi2(term):
    """For the RDC terms associated with the argument VarTensor object,
    return composite chi^2 by calling
    <m rdcPotTools>.composite_chi62."""
    from simulationTools import getPotTerms
    rdcTerms = getPotTerms(term.expts(),'RDCPot1')
    
    from rdcPotTools import composite_chi2
    return composite_chi2(rdcTerms)

import simulationTools
simulationTools.registerTerm(VarTensor_analyze,
                             "Variable Tensor Analysis","VarTensor",
                             r"""
Printed are values of the Da and rhombicity (Rh) for one or more alignment
tensors averaged over the ensemble.

For calculations in which each ensemble member has its own alignment tensor,
details are printed for each.

resid: the residue number for the pseudo atoms describing the alignment tensor.

GDO: generalized degree of order. [Ref: Tolman et al. JACS 123, 1416 (2001)]:

.. math::

  \sqrt{\frac{2}{3} \sum_{i,j} S_{i,j}^2},

where :math:`S_{i,j}` are the elements of the Saupe matrix.

Euler angles for four equivalent orientations are computed using Goldstein's
x-convention [Goldstein, "Classical Mechanics" p. 147 ]

If more than two tensors are present, the normalized scalar products are
reported [Sass et al. JACS 121, 2047 (1999).]


                             """)

simulationTools.registerExtraStats("VarTensor","RDC R-fac",
                                   composite_RDCRfactor)
simulationTools.registerExtraStats("VarTensor","RDC Chi^2",
                                   composite_RDCchi2)

def xyzFromEuler(tensor,phi,theta,psi):
    """set the x, y, and z axes of tensor to correspond to the
    given Euler phi, theta and psi rotations (in radians)
    """

    Da = tensor.Da()
    Rh = tensor.Rh()
    
    from vec3 import Vec3
    from mat3 import rotVector
    D = rotVector([0,0,1],phi)
    C = rotVector([1,0,0],theta)
    B = rotVector([0,0,1],psi)
    xyz = B * C * D

    from mat3 import transpose
    xyz = transpose(xyz)

    x = xyz * Vec3(1,0,0) + tensor.oAtom().pos()
    y = xyz * Vec3(0,1,0) + tensor.oAtom().pos()
    z = xyz * Vec3(0,0,1) + tensor.oAtom().pos()

    tensor.xAtom().setPos(x)
    tensor.yAtom().setPos(y)
    tensor.zAtom().setPos(z)

    tensor.setDa(Da)
    tensor.setRh(Rh)

    return

