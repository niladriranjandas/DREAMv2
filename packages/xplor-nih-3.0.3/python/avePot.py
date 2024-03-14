
import avePot1

from os import environ as env

def help():
    return open(env["XPLOR_DIR"]+"/helplib/nih-py-avePot").read()


class AvePot(avePot1.realAvePot):
    """
    wrapper around avePot1.AvePot:
      construct from either an already created Pot:
        pot = XplorPot("BOND")
        AvePot(pot)
      or internally call the averaged pot's constructor
        AvePot(XplorPot,"BOND")
      where the first argument is the class constructor and the remaining
      are parameters for that constructor.

      __getattr__ is defined such that first the avePot1.AvePot attributes are
      used, and then those of the subPot are attempted.
    """
    def __init__(self,*args):
        if len(args)>1:
            from ensembleSimulation import EnsembleSimulation_currentSimulation
            sim = EnsembleSimulation_currentSimulation()
            if not sim:
                from simulation import Simulation_currentSimulation
                sim = Simulation_currentSimulation()
                subSim = sim
            else:
                subSim = sim.subSim()
                pass
            pot=args[0](*args[1:]+(subSim,))
            pass
        else:
            pot = args[0]
            pass
                 
        try:
            avePot1.realAvePot.__init__(self,pot.potObj)
        except AttributeError:
            avePot1.realAvePot.__init__(self,pot)
            pass

        return
    def __getattr__(self,name):
        if name=="this":
            if hasattr(avePot1.realAvePot,"this"):
                return avePot1.realAvePot.this
            else:
                return None
            pass
        try:
            return avePot1.realAvePot.__getattr__(self,name)
        except AttributeError:
            return eval("self.subPot().%s" % name)
        return
    pass

realAvePot = AvePot
def AvePot(*args):
    from potProxy import PotProxy
    return PotProxy( realAvePot(*args) )


class AvePotPtr(avePot1.AvePotPtr):
    """
    this so that potDerive.createDerivedPot works correctly.
    wrapper around avePot1.AvePotPtr:
      construct from a pointer to avePot1.AvePot

      __getattr__ is defined such that first the avePot1.AvePot attributes are
      used, and then those of the subPot are attempted.
    """
    def __init__(self,pointer):
        avePot1.AvePotPtr.__init__(self,pointer)
        return
    def __getattr__(self,name):
        try:
            return avePot1.AvePotPtr.__getattr__(self,name)
        except AttributeError:
            return eval("self.subPot().%s" % name)
        return
    pass

