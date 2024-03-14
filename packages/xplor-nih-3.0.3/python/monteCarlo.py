"""perform some Monte-Carlo operations.

see randomizeTorsions below.
"""

def randSign():
    from random import choice
    return choice([1,-1])

from ivm import IVM
class MC( IVM ):  #doesn't pass the is-a test
    """analog of the <m ivm>.IVM class which performs raw Monte Carlo
    instead of dynamics/minimization.

    NOTE: untested!
    """
    def __init__(s,simulation):
        IVM.__init__(s,simulation)

        s.setIntegrateType( "powell" ) # so it's it minimization mode
        s.type_ = "MonteCarlo"
        s.setStepsize( 0.01 )

        s.stepsizeA_        = None
        s.stepsizeInfo_     = None
        s.Ra0_              = 0.5   # target acceptance ratio
        s.tau_              = 20    # response time of stepsize adjuster
        return

    #accessors
    def setStepsizeA(s,v): s.stepsizeA_ = v
    def stepsizeA(s): return s.stepsizeA_
    def setRa0(s,v): s.Ra0_ = v
    def Ra0(s): return s.Ra0_
    def setTau(s,v): s.tau_ = v
    def tau(s): return s.tau_
    def type(s): return s.type_
    

    def __del__(s):
        #s.posTrajFile.close()
        #s.velTrajFile.close()
        return

    def help(s):
        print("FIX: help")
        return

    def info(s):
        print("FIX: current settings")
        return

    def run(s):
        import random
        if s.numSteps()<1:
            return
        if s.stepsize() <= 0.:
            print("stepsize=0. No steps taken.")
            return;
        try:
            #for eaach coordinate:
            # 1) update it (random size, random magnitude w/in its range)
            # 2) evaluate energy accept if it passes MC criteria
            # 3) update range s.t. the acceptance ratio is %50

            if s.reuse():
                s.initDynamics(1)
            else:
                s.init()
                if (type(s.stepsizeA()).__name__ != "list" or
                    len(s.stepsizeA()) != len(s.pos())):
                    s.setStepsizeA( [0] * len(s.pos()) )
                    pass
                for i in range(len(s.stepsizeA())):
                    if s.stepsizeA_[i]==0:
                        s.stepsizeA_[i] = s.stepsize()
                        pass
                    pass
                pass
   
            if s.resetCMInterval()>0: s.resetCM()

            s.outputStepInfo(0,s.stepsize(),s.type());
            eCount=0


            done=None
            iter=0
            time=0.0
            beta=1.0 / (s.bathTemp() * s.simulation.kBoltzmann())
            coords = s.pos()
            s.calcEnergy()
            energy = s.Epotential()
            print("energy: " + repr(energy))
            minEnergy = energy
            print(coords)
            minCoords = coords
            s.stepsizeInfo_ = coords[:]
            for i in range(len(coords)):
                s.stepsizeInfo_[i] = [0,0]
                pass
            
            while not done:
                iter += 1

                if  s.resetCMInterval()!=0 and iter%s.resetCMInterval()==0:
                    s.resetCM();
                    pass
     
                for i in range(len(coords)):
                    testCoords = coords[:]
                    testCoords[i] += randSign() * s.stepsizeA_[i]
                    s.setPos( testCoords )
                    s.calcEnergy()
                    testEnergy = s.Epotential()
                    accepted=0
                    if (testEnergy < energy or 
                        exp(-beta*(testEnergy-energy)) > random.uniform(0,1)):
                        accepted=1
                        energy = testEnergy
                        coords = testCoords
                        pass
                    s.stepsizeA_[i] *= s.updateStepSize(i,accepted)
                    if energy < minEnergy:
                        minEnergy = energy
                        minCoords = coords
                        pass
                    pass
                pass
            
            
                if  s.numSteps()>0 and iter>=s.numSteps():
                    done=1
                    pass
     
                if iter%s.printInterval()==0:
                    s.outputStepInfo(iter,s.stepsize(),s.type())
                    pass
                pass
            s.setPos( minCoords )
            s.calcEnergy()
            print("MonteCarlo: minimum energy: %f %f" % (minEnergy,
                                                         s.Epotential()))
        except:
            import traceback
            print('caught exception: ')
            traceback.print_exc()
            raise
            pass
        return
    def updateStepSize(s,dof,accepted):
        #FIX: if the initial stepsize is way off, stepsize equilibration
        # will take a long time. Better: acceptance ratio should use only
        # most recent data.
        from math import sqrt
        s.stepsizeInfo_[dof][0] += 1
        s.stepsizeInfo_[dof][1] += accepted
        fac = 1
        if s.stepsizeInfo_[dof][0] >= 10:
            Ra = float(s.stepsizeInfo_[dof][1]) / s.stepsizeInfo_[dof][0]
            fac = sqrt(1 - (s.Ra0()-Ra) / (s.tau() * s.Ra0()))
            print("updateStepSize: %d: %d: %f Ra: %f" % \
                  (dof, s.stepsizeInfo_[dof][0], s.stepsizeA_[dof], Ra))
            pass
        return fac
    pass


from atomSel import AtomSel
import random

def randomizeTorsions(ivm,
                      sel="all",
                      range=180):
    """Randomize torsion angles in the input ivm (an <m ivm>.IVM instance).

    Applicable torsion angles are those whose two center atoms are contained in
    the sel argument (an XPLOR selection string).  Note that the ivm
    configuration is respected so that rigid regions are not affected.  The
    argument range (a number; degrees) specifies the maximum +/- change in value
    of each torsion angle.
    
    """

    from math import pi
    pmAngle = range * pi / 180
    if type(sel).__name__ == type("s").__name__:
        sel = AtomSel(sel)
        pass

    pos = ivm.pos()
    for node in ivm.nodeList():
        if node.type() == "torsion" and \
           sel.containsIndex( node.atoms()[0] ) and \
           sel.containsIndex( node.parentAtom() ):
            pos[node.startIndex()] = random.uniform(-pmAngle,pmAngle)
            pass
        pass
    ivm.setPos( pos )
    return

        
        
