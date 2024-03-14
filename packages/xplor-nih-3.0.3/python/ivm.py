#!/usr/bin/env python

from publicIVM import PublicIVM
from publicIVM import PublicIVM
from publicIVM import IVMError
from atomSel import AtomSel
from selectTools import convertToSelectionIndices
import re
import string
import sys
import simulation
#from trajFile import TrajFile

from os import environ as env

def help():
    return open(env["XPLOR_DIR"]+"/helplib/nih-py-ivm").read()
pyXplorHelp = help

class IVM(PublicIVM):
    """ Internal Variable Module class
          members:
            reuse: reuse coordinates and topology from previous call to 
                    run(), init(), or autoTorsion()
                   this can be unset explicitly, or whenever a topology
                   element is modified."""
    def __init__(s,sim=0):
        """
        Constructor. The optional argument is a <m simulation>.Simulation. By
        default, the current simulation is used.
        """
        if not sim: sim = simulation.currentSimulation()
        PublicIVM.__init__(s,sim)
        s.simulation        = sim
        s.write             = sys.stdout.write
        s.stepsize_         = 0.01
        s.setMaxDeltaE(10000)
        s.setResponseTime(20)
        s.numSteps_         = 0     #number of dynamics/minimization steps
        s.resetCMInterval_  = 0     #how often to zero CM velocity/AM
        s.reuse_            = 0     #boolean: reuse topology and traj. info
        s.finalTime_        = 0     #time duration for integration
        s.printInterval_    = 1  #how often to print energy info
        s.stepsizeThreshold_= None  #if stepsize falls below this value
                                    # stepsizeEscape is called
                                    # try setting to 1e-6
        s.time_=0                 
        s.iter_=0                 
        s.stepsizeEscape = minimizeEscape
        s.setStepType("PC6")
        s.verbosePairs = []
        for f in ("printCoords",     	 
                  "printResetCM",    	 
                  "printVelFromCartCost", 
                  "printTemperature",     
                  "printEnergy",    	 
                  "printCMVel",     	 
                  "printNodeForce", 	 
                  "printNodePos",   	 
                  "printNodeTheta", 	 
                  "printStepDebug", 	 
                  "printStepInfo",  	 
                  "printNodeDef",   	 
                  "printLoopDebug", 	 
                  "printLoopInfo"):
            exec ( "s.verbosePairs.append( (s.%s, '%s') )" % (f,f) )
            pass
        
        s.configActions=[]

        #FIX: deal with the next items
        s.trajectory_        = None
        #s.velTrajFile       = TrajFile()
        #s.velTrajFile.type  = "velocity"
        return

    #def __del__(s):
    #    #s.posTrajFile.close() : close is called in trajFile destructor
    #    #s.velTrajFile.close()
    #    PublicIVM.__del__(s)
    #    return

    def addConfigAction(s,action):
        """
        Add a command to be run before topology configuration via
        protocol.torsionTopology or protocol.cartesianTopology.
        """
        from simulationTools import IVMAction
        if not isinstance(action,IVMAction):
            raise Exception("action must be an instance of " +
                            "simulationTools.IVMAction")
        
        s.configActions.append(action)
        return

    def help(s):
        return open(env["XPLOR_DIR"]+"/helplib/nih-py-ivm").read()

    def idAtom(s,i):
        sim = s.simulation
        ret = "(%d)" % i
        if i<0 or i>=sim.numAtoms():
            ret += " [unknown]"
        else:
            ret += " " + sim.atomName(i)
            ret += " " + sim.residueName(i)
            ret += " " + repr(sim.residueNum(i))
            ret += " " + sim.segmentName(i)
            pass
        return ret

    def fixedAtoms(s):
        """Return an <m atomSel>.AtomSel containing all atoms fixed in space.
        """
        indices=[]
        for g in s.groupList():
            if -1 in g:
                for i in g:
                    if i>=0: indices.append(i)
                    pass
                pass
            pass
        from atomSel import AtomSel
        return AtomSel(indices,s.simulation)
        

    def time(s):
        """
        Actual time elapsed during dynamics. It is reset and updated by
        run().
        """
        return s.time_

    def iters(s):
        """
        Actual number of steps taken during dynamics or minimization. It is
        reset and updated by run().
        """
        return s.iter_

    

    def info(s,what="default"):
        """return a string with info on the IVM settings.
        An optional argument specifies information to return:
           default
           integrate    - integration specific info
           minimization - minimization specific info
           topology     - information about the topology settings"""
        
        ret = "current settings for the Internal Variable Module:\n"

        if what=="default" or what=="all":
            ret += "\tStep type: " + s.stepType()
            try:
                if s.minimization(): ret += "[minimization]"
                else: ret += "[integration]"
            except:
                pass
            ret += '\n'
            ret += "\tstepsize: " + repr(s.stepsize()) + "\n"
            ret += "\tnumSteps: " + repr(s.numSteps()) + '\n'
            ret += "\teTolerance: " +  repr(s.eTolerance()) + '\n'
            ret += "\tconstrainLengths: " + repr(s.constrainLengths())+'\n'
            ret += "\tprintInterval: " + repr(s.printInterval()) + '\n'
            ret += "\tverbose info to print:\n"
            for v in s.verbosePairs:
                if v[0]&s.verbose(): ret += "\t\t" + v[1] + '\n'
                pass
            ret += "\ttrajectory: " + repr(s.trajectory()) + '\n'
            # ostr << "VELocity-file= " << velFilename << '\n';
            ret += "\tnumber of internal degrees of freedom: %d\n" % s.dof()
            pass
        

        if what.find("integ")==0 or what=="all":
            ret += "\n\tIntegration parameters:\n"
            ret += "\t\tfinalTime: " + repr(s.finalTime()) + "\n"
            ret += "\t\tmaxTSFactor: " + repr(s.maxTSFactor()) + '\n'
            ret += "\t\tmaxDeltaE: " + repr(s.maxDeltaE()) + '\n'
            ret += "\t\tresponseTime: " + repr(s.responseTime()) + '\n';
            ret += "\t\tresetCMInterval: " + repr(s.resetCMInterval()) + '\n'
            ret += "\t\tscaleVel: " + repr(s.scaleVel()) + '\n'
            ret += "\t\tadjustStepsize: " + repr(s.adjustStepsize()) + '\n'
            ret += "\t\tbathTemp: " + repr(s.bathTemp()) + '\n'
            ret += "\t\tfrictionCoeff: " + repr(s.frictionCoeff()) + '\n'
            pass

        if what.find("min")==0 or what=="all":
            ret += "\n\tMinimization parameters:\n"
            ret += "\t\tdEpred: " + repr(s.dEpred()) + '\n'
            ret += "\t\tgTolerance: " + repr(s.gTolerance()) + '\n'
            ret += "\t\tmaxCalls: " + repr(s.maxCalls()) + '\n'
            pass

        if what.find("topo")==0 or what=="all":
            ret += "\n\tTopology parameters:\n"
            ret += "\t\treuse: " + repr(s.reuse()) + '\n'
            ret += "\t\tBonds to ignore:\n"
            for b in s.bondExclude():
                ret += "\t\t\t" + s.idAtom( b[0] ) + " " + s.idAtom( b[1] ) + '\n'
                pass
            ret += "\t\tBonds to constrain fixed:\n"
            for b in s.constraintList():
                ret += "\t\t\t" + s.idAtom( b[0] ) + " " + s.idAtom( b[1] ) + '\n'
                pass
            ret += "\t\tAtoms fixed in space:\n"
            for g in s.groupList():
                if -1 in g:
                    for a in g: ret += "\t\t\t" + s.idAtom(a) + '\n'
                    pass
                pass
            ret += "\t\tAtoms in base nodes:\n"
            for a in s.baseAtoms():
                for a in g: ret += "\t\t\t" + s.idAtom(a) + '\n'
                pass

            ret += "\t\tatoms held rigid wrt each other:\n"
            cnt=0
            for g in s.groupList():
                if not -1 in g:
                    ret += "\t\t\tgroup %d (%d atoms):\n" % (cnt,len(g))
                    for a in g: ret += "\t\t\t\t" + s.idAtom(a) + '\n'
                    cnt += 1
                    pass
                pass
            ret += "\t\thinge specifications:\n"
            for h in s.hingeList():
                ret += "\t\t\t" + h[0]
                if len(h)>2: ret += " [" + s.idAtom(h[2]) + "] "
                if len(h)>3: ret += " [" + s.idAtom(h[3]) + "]"
                ret += '\n'
                for a in h[1]: ret += "\t\t\t\t" + s.idAtom(a) + '\n'
                pass
            ret += "\t\tnumber of internal degrees of freedom: %d\n" % s.dof()
            ret += "\t\tdimension of the state vector: %d\n" % s.dim()
            pass
        return ret
        
    #accessors

    def baseAtoms(s):            return s.oldBaseAtoms()
    def setBaseAtoms(s,l):
        s.reuse_=0
        oIndices=convertToSelectionIndices(l,s.simulation)
        indices=set( oIndices )
        s.setOldBaseAtoms( list(indices) )
        return
    def stepsize(s):             return s.stepsize_
    def setStepsize(s,v):        s.stepsize_ = v
    def numSteps(s):             return s.numSteps_
    def setNumSteps(s,v):        s.numSteps_ = v
    def resetCMInterval(s):      return s.resetCMInterval_
    def setResetCMInterval(s,v): s.resetCMInterval_ = v
    def reuse(s):                return s.reuse_
    def resetReuse(s):           s.reuse_ = 0
    def finalTime(s):            return s.finalTime_
    def setFinalTime(s,v):       s.finalTime_ = v
    def printInterval(s):        return s.printInterval_
    def setPrintInterval(s,v):   s.printInterval_=v
    def trajectory(s):           return s.trajectory_
    def setTrajectory(s,v):      s.trajectory_ = v
    def stepsizeThreshold(s):      return s.stepsizeThreshold_
    def setStepsizeThreshold(s,v): s.stepsizeThreshold_ = v
    def minStepsize(s):            return s.minStepSize()
    def setMinStepsize(s,v):       s.setMinStepSize(v)

    def setBondExclude(s,v):    s.reuse_=0; PublicIVM.setBondExclude(s,v)
    def setGroupList(s,v):      s.reuse_=0; PublicIVM.setGroupList(s,v)
    def setConstraintList(s,v): s.reuse_=0; PublicIVM.setConstraintList(s,v)
    def setHingeList(s,v):      s.reuse_=0; PublicIVM.setHingeList(s,v)
    def setStepType(s,v):       s.reuse_=0; PublicIVM.setStepType(s,v)

    def reset(s):
        s.setBondExclude( [] )
        s.setGroupList( [] )
        s.setHingeList( [] )
        s.setConstraintList( [] )
        s.setBaseAtoms( [] )
        #s.fixAtomsFromSim()
        #fixAtomsFromExternal(natom,imove); # read atoms fixed externally
        return

    def fix(s,sel):
        sel = convertToSelectionIndices(sel,s.simulation)
        aList = [-1]
        aList += sel
        s.setGroupList( s.groupList() + [aList] )
        return

    def group(s,sel):
        sel = convertToSelectionIndices(sel,s.simulation)
        l = s.groupList()
        l += [sel,]
        s.setGroupList( l )
        return

    def groups(s):
        """Return a list of <m atomSel>.AtomSel objects corresponding to
        atoms grouped in rigid bodies in this IVM object.

        [This method converts the indices from groupList() into AtomSel
        objects, discarding negative indices which are used for IVM internal
        purposes.]
        """
        from atomSel import AtomSel
        sim = s.simulation
        ret = [AtomSel(indices,sim) for indices in s.simulationGroupList()]
        return ret

    def hinge(s,hingeType, sel=None, *sels):
        """add a hinge definition of type hingeType for atoms specified by
        sel, an <m atomSel>.AtomSel, string or list of indices. Some
        hinge types require additional selections specified as addtional
        arguments.        
        """

        if sel is None:
            # this clause for backwards compatibility
            h = hingeType
            hingeType, sel = (h[0], h[1])
            if len(h)>2:
                sels += h[2:]
                pass
            pass
                
        entry = [hingeType, convertToSelectionIndices(sel,s.simulation)]
        for sel in sels:
            entry.append(convertToSelectionIndices(sel,s.simulation))
            pass

        new = s.hingeList() + [entry]
        s.setHingeList( new )
        return

    def breakBond(s,*args):
        """ Argument is either two selections specifying one atom each or
            a single selection specifying two atoms corresponding to a 
            bond to break (topologically). 
            The BOND potential is not altered."""
        if len(args)==1:
            atoms=convertToSelectionIndices(args[0],s.simulation)
        elif len(args)==2:
            atoms=(convertToSelectionIndices(args[0],s.simulation),
                   convertToSelectionIndices(args[1],s.simulation))
        else:
            raise Exception("args must specify one or two atom selections")
        s.setBondExclude( s.bondExclude() + [atoms] )
        return
    def constrainBond(s,el):
        """ specify selection consisting of two atoms corresponding to a
            bond to constrain"""
        s.breakBond(el)
        s.setConstraintList( s.constraintList() +
                          [convertToSelectionIndices(el,s.simulation)] )
        return

    def breakAllBondsTo(s,sel):
        """break all bonds to bonded atoms in selection"""
        sel = convertToSelectionIndices(sel,s.simulation)
        for id in sel:
            for id2 in AtomSel("bondedto id %d" % (id+1),
                               s.simulation).indices():
                s.breakBond( (id, id2) )
                pass
            pass
        return

    def breakAllBondsIn(s,sel):
        """break all bonds between atoms in selection"""
        sel = convertToSelectionIndices(sel,s.simulation)
        PublicIVM.breakAllBondsIn(s,sel)
        return

    #def fixAtomsFromSim(s):
    #    gList = [-1]
    #    for i in s.sim.constraintList:
    #        gList.append(i)
    #        pass
    #    s.groupList.append( gList )
    #    return

    def autoTorsion(s):
        "set up hinges for torsion angle dynamics/minimization"

        #fix: this could support an atom selection argument if the next
        # function did. Or- it's possible that those atoms not in the sel
        # and whose group is not set could be sent to a temporary group
        # and then unset afterwards.
        s.groupTorsion()
        
        hingeList=s.hingeList()

        #get all grouped atoms into a hingeList
        # if one a hinge is specified for any atom in the group
        for group in s.groupList():
            for g in group:
                for (h,hinge) in enumerate(hingeList):
                    if hinge[0].startswith('bend'):
                        # don't mess with these hinges-
                        # they are modified by the ivm
                        continue
                    hingeAtoms=set(hinge[1])
                    if g in hingeAtoms:
                        hingeAtoms=hingeAtoms.union(set(group))
                        hingeList[h]=(hingeList[h][0],sorted(list(hingeAtoms)))
                        break
                    pass
                pass
            pass

        # any atoms not now in a hinge should be specified as
        # torsion hinges
        ids = []
        for i in range(0, s.simulation.numAtoms()):
            inHinge=0;
            for hinge in hingeList:
                if i in hinge[1]:
                    inHinge=1;
                    break
                pass
                
            if not inHinge:
                ids.append(i)
                pass
            pass
        
        if len(ids):
            hingeList.append( ("torsion" , ids) )
            pass

        s.setHingeList( hingeList )

        s.init() # FIX: necessary?
        return

    def outputStepInfo(s, step, stepsize, type, totTime):
        #print "timefactor: " + `s.simulation.timeFactor`
        #print "stepsize: " + `stepsize`
        if s.minimization():
            s.write("*-- %s" % type.upper() + " " + "-"*(10-len(type)))
            s.write("-- step=%7d" % step)
            s.write(" --- stepsize=%10.5f" % s.stepsize())
            s.write(" --- energy evals=%5d -*\n" % s.eCount())
        else:
            s.write("*--- Dynamics ---- step=%7d" % step)
            s.write(" ---- time=%10.5g" % totTime)
            s.write(" ---- delta_t=%10.5g --*\n" % stepsize)

        if not s.minimization():
            s.write("| E(kin)+E(poten)=%14.3f" % s.Etotal())
            s.write(" E(kin)=%14.3f" % s.Ekinetic())
            s.write(" temperature=%11.3f |\n" % s.currentTemp())
            pass
        
        terms = s.potList().energyReports()
        for i in range(0,len(terms)):
            if terms[i][0] == "SCRI":
                del terms[i]
                break
            pass

        terms = [("E(poten)",s.Epotential(),0),
                 ("grad",s.gradMagnitude(),0)] + terms

        fieldFormat = r"%8s=%14.7f"
        for i in range(0,len(terms),3):
            s0 =  fieldFormat% (terms[i][0], terms[i][1])
            s1=""
            if i+1 < len(terms): s1 = fieldFormat % (terms[i+1][0],
                                                      terms[i+1][1])
            s2 = ""
            if i+2 < len(terms): s2 = fieldFormat % (terms[i+2][0],
                                                     terms[i+2][1])
            s.write("| %25s %25s %24s |\n" % (s0,s1,s2) )
            pass
        s.write("*" + "-"*78 + "*\n")
        s.printStepDetails();
        s.eCountReset() #reset
        return

    #def initDynamics(s,bool):
    #    s.be.initDynamics(1)
    #    return

    def init(s):
        #s.be.setGroupList( s.groupList )
        #s.be.setConstraintList( s.constraintList )
        #s.be.setBondExclude( s.bondExclude )
        #s.be.setHingeList( s.hingeList )
        #s.be.setOldBaseAtoms( s.baseList )
        PublicIVM.init(s)
        s.minimizationFlag   = s.minimization()
        s.reuse_ = 1
        #s.groupList      = s.be.groupList()
        #s.constraintList = s.be.constraintList()
        #s.bondExclude    = s.be.bondExclude()
        #s.hingeList      = s.be.hingeList()
        #s.baseList       = s.be.oldBaseAtoms()
        return

    def calcEnergy(s):
        """calculate energy in internal coordinates, and return the total
        energy. This is a safe wrapper around the underlying calcEnergy"""
        s.init()
        PublicIVM.calcEnergy(s)
        return s.Etotal()

    def run(s):
        """ Perform dynamics or minimization. Stop when numSteps have been
        taken, or finalTime has been reached (whichever is first).

        Returns an integer which is 1 on successfull completion of dynamics or
        minimization and negative on failure.
        """

        
        for atom in AtomSel("not known",s.simulation):
            print("WARNING: large displacement for atom %s:" %atom.string(),\
                                                              atom.pos())
            pass

        if s.numSteps()<1 and s.finalTime()<=0.:
            return
        if s.stepsize() <= 0.:
            print("stepsize=0. No steps taken.")
            return;
        s.initStepsize = s.stepsize()
        if s.resetCMInterval()>0: s.resetCM()

        from coordComparer import CoordComparer
        #            s.syncPosVelFrom()
        if not s.reuse():
            s.init()
            pass

        if s.bathTemp()==0. and s.scaleVel() and not s.minimization():
            print("bathTemp=0:")
            print("  This is not a valid option when scaleVel is set.")
            return

        coordComparer=CoordComparer(s.simulation)
        s.potList().updateDelta(0,coordComparer.update())

        s.initDynamics(1)

        if s.printInterval():
            s.outputStepInfo(0,s.stepsize(),s.stepType(),0.);

        if s.trajectory():
            s.trajectory().write(0,0.)
            pass

        done=None
        s.iter_=0
        s.time_=0.0
        while not done:
            s.iter_ += 1

            if  s.resetCMInterval()!=0 and s.iters()%s.resetCMInterval()==0:
                s.resetCM();
                pass

            s.time_ += s.stepsize();

            s.potList().updateDelta(s.iters(),coordComparer.update())

            try:
                (done,s.stepsize_) = PublicIVM.step(s,s.stepsize());
            except:

                if (s.stepsizeThreshold()!= None and
                    s.stepsize_<s.stepsizeThreshold()):
                    print("IVM: attempting to avoid small stepsize...")
                    s.stepsizeEscape(s)
                    s.initDynamics(1) #reset velocities
                    s.setStepsize(s.initStepsize)
                    pass
                else:
                    raise
                    pass
                pass

            if not s.minimization():
                if s.finalTime()>0.0 and s.time() >= s.finalTime():
                    done=1
                    pass
                pass
            if  s.numSteps()>0 and s.iters()>=s.numSteps():
                done=1
                pass

            if s.trajectory():
                s.trajectory().write(s.iters(),s.time())
                pass

            #if s.trajFile and nsavc>0 and iter%nsavc==0:
            #    s.trajFile.writeCoords()
            #    pass
            #
            #if s.velFile and nsavv>0 and iter%nsavv==0:
            #    s.velFile.writeVel()
            #    pass

            if s.printInterval() and s.iters()%s.printInterval()==0:
                s.outputStepInfo(s.iters(),s.stepsize(),s.stepType(),s.time())
                pass

            if done<0: break

            pass

        #FIX: update xplor variables ???
        #   const char* dum = "";
        #   FORTRAN(declar)("DINT_TIME" ,"DP",
        #                   dum,dum,time*xplorVars->timeFac,9,2,0,0);
        #   FORTRAN(declar)("DINT_STEPS","DP",dum,dum,(double)iter,10,2,0,0);
        #   if ( coordUnit>0 ) FORTRAN(vclose)(coordUnit,"KEEP",error,4);


        if s.printInterval() and s.iters()%s.printInterval()!=0:
            s.outputStepInfo(s.iters(),s.stepsize(),s.stepType(),s.time())
            pass
        #        except ivmBE.error:
        #            print "error in ivmBE"
        #            raise
        #except:
        #    import traceback
        #    print 'caught exception: '
        #    traceback.print_exc()
        #    pass
        return done
    pass

#
# these used to correct for small stepsize
#

minimizeEscapeCnt=0
def minimizeEscape(ivm):
    escMin=0
    if not escMin:
        escMin = IVM(ivm.simulation)
        escMin.setBaseAtoms( ivm.baseAtoms() )
        escMin.setBondExclude( ivm.bondExclude() )
        escMin.setGroupList( ivm.groupList() )
        escMin.setConstraintList( ivm.constraintList() )
        escMin.setHingeList( ivm.hingeList() )
        escMin.setConstrainLengths( ivm.constrainLengths() )
        escMin.setPotList( ivm.potList() )
        pass
        
    import protocol
    protocol.initMinimize(escMin,numSteps=20)
    global minimizeEscapeCnt
#    if minimizeEscapeCnt>6:
#        from simulationTools import testGradient
#        if not testGradient(escMin.potList(),eachTerm=1):
#            raise Exception("bad gradient")
#        pass
    minimizeEscapeCnt+=1
    
    for i in range(1):
        escMin.run()
        pass
    return
    

