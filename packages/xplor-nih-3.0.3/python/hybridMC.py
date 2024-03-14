"""mixed Monte-Carlo molecular dynamics minimization

Note: untested!
"""
import random
import xplor

from math import sqrt
from math import exp


sim = xplor.simulation

debug=1


def printMem(label=""):
    #print "musage hmc:" + label + ' ' + \
    #      ` open("/proc/self/statm").read().split()[0]`
    return

def maxwellDist(temp,mass):
    return (random.gauss( 0, sqrt(temp*sim.kBoltzmann()/mass) ),
            random.gauss( 0, sqrt(temp*sim.kBoltzmann()/mass) ),
            random.gauss( 0, sqrt(temp*sim.kBoltzmann()/mass) ))

def mc(ivm,temp,numMCSteps,targetAcceptRatio=0.5):
    """given an <m ivm>.IVM setup for dynamics, perform the specified number
    of Monte-Carlo steps.

    The targetAcceptRatio argument is used to update the ivm's finalTime
    value for a possible subsequent call to this function.
    """
    from xplorSimulation import getXplorSimulation
    sim = getXplorSimulation()
    
    beta = 1 / (temp * sim.kBoltzmann())
    acceptCnt=0
    ivm.init()
    printMem()
    for cnt in range(0,numMCSteps):
        #  a) generate random momenta consistent w/ temp.
        for i in range(0,sim.numAtoms()):
            sim.setAtomVel(i, maxwellDist(temp,sim.atomMass(i)))
            pass
        printMem("0")
        ivm.velFromCartesian()
        printMem("1")
        ivm.initDynamics(1)
        printMem("2")
        energy = ivm.Epotential()
        printMem("3")
        oldAtomPosArr = sim.atomPosArr()
        print("musage: " + repr(oldAtomPosArr))
        printMem("4")
    
        #  b) take (random??) number of VV MD steps
        #ivm.setNumSteps(numMCSteps)
        printMem("5")
        ivm.run()
        printMem("6")
    
        #  c) accept or reject this configuration based on the metropolis
        #     criterion, based on the total energy (pot+kinetic).
        testEnergy = ivm.Epotential()
        if debug: print("mc: testEnergy: %f   energy: %f   temp: %f" % \
           (testEnergy,energy,temp))
        if (testEnergy < energy or 
            exp(-beta*(testEnergy-energy)) > sim.uniformRandom() ):
            acceptCnt += 1
        else:
            sim.setAtomPosArr( oldAtomPosArr )
        printMem("7")
        pass
    acceptRatio = float(acceptCnt) / numMCSteps
    print("mc: acceptance ratio: %f" % acceptRatio)
    ivm.setFinalTime( sqrt(1 + max(-1,(acceptRatio-targetAcceptRatio)/0.6)) *
                      ivm.finalTime() )
    print("mc: trajectory duration: %f" % ivm.finalTime())
    return
