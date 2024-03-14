# This file was automatically generated by SWIG (http://www.swig.org).
# Version 3.0.12
#
# Do not make changes to this file unless you know what you are doing--modify
# the SWIG interface file instead.

from sys import version_info as _swig_python_version_info
if _swig_python_version_info >= (2, 7, 0):
    def swig_import_helper():
        import importlib
        pkg = __name__.rpartition('.')[0]
        mname = '.'.join((pkg, '_ensembleSimulation')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_ensembleSimulation')
    _ensembleSimulation = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_ensembleSimulation', [dirname(__file__)])
        except ImportError:
            import _ensembleSimulation
            return _ensembleSimulation
        try:
            _mod = imp.load_module('_ensembleSimulation', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _ensembleSimulation = swig_import_helper()
    del swig_import_helper
else:
    import _ensembleSimulation
del _swig_python_version_info

try:
    _swig_property = property
except NameError:
    pass  # Python < 2.2 doesn't have 'property'.

try:
    import builtins as __builtin__
except ImportError:
    import __builtin__

def _swig_setattr_nondynamic(self, class_type, name, value, static=1):
    if (name == "thisown"):
        return self.this.own(value)
    if (name == "this"):
        if type(value).__name__ == 'SwigPyObject':
            self.__dict__[name] = value
            return
    method = class_type.__swig_setmethods__.get(name, None)
    if method:
        return method(self, value)
    if (not static):
        object.__setattr__(self, name, value)
    else:
        raise AttributeError("You cannot add attributes to %s" % self)


def _swig_setattr(self, class_type, name, value):
    return _swig_setattr_nondynamic(self, class_type, name, value, 0)


def _swig_getattr(self, class_type, name):
    if (name == "thisown"):
        return self.this.own()
    method = class_type.__swig_getmethods__.get(name, None)
    if method:
        return method(self)
    raise AttributeError("'%s' object has no attribute '%s'" % (class_type.__name__, name))


def _swig_repr(self):
    try:
        strthis = "proxy of " + self.this.__repr__()
    except __builtin__.Exception:
        strthis = ""
    return "<%s.%s; %s >" % (self.__class__.__module__, self.__class__.__name__, strthis,)


def _swig_setattr_nondynamic_method(set):
    def set_attr(self, name, value):
        if (name == "thisown"):
            return self.this.own(value)
        if hasattr(self, name) or (name == "this"):
            set(self, name, value)
        else:
            raise AttributeError("You cannot add attributes to %s" % self)
    return set_attr



def fromSimulation(*args, **kwargs):
    return _ensembleSimulation.fromSimulation(*args, **kwargs)
fromSimulation = _ensembleSimulation.fromSimulation

def memberFromSimulation(*args, **kwargs):
    return _ensembleSimulation.memberFromSimulation(*args, **kwargs)
memberFromSimulation = _ensembleSimulation.memberFromSimulation

def singleThread(index=0):
    """
    Wrap around <m ensembleSimulation>.EnsembleSimulation.singleThread
    This returns true for the ensemble member of an EmsembleSimulation 
    specified by the index argument or if no EnsembleSimulation is used.
    """
    esim = EnsembleSimulation_currentSimulation()
    import simulation
    if esim and simulation.currentSimulation().name() == esim.name():
        return esim.singleThread(index)
    else:
        return 1
    pass

def multiThread():
    """
    Wrap around <m ensembleSimulation>.EnsembleSimulation.multiThread
    which also works if no EnsembleSimulation is used.
    """
    esim = EnsembleSimulation_currentSimulation()
    import simulation
    if esim and simulation.currentSimulation().name() == esim.name():
        esim.multiThread()
        pass
    return

def commBarrier(comm):
    """ perform EnsembleSimulation-safe <m socketComm>.SocketComm barrier

	It's probably a better idea to instead create and use an instance of
	the local Comm class.
    """
    esim = EnsembleSimulation_currentSimulation()

    if not esim:
        procs = comm.barrier()
        return procs

    from ensembleSharedObj import SharedObj
    sharedProcs = SharedObj()
    if singleThread():
        procs = comm.barrier()
        sharedProcs.set(procs)
        pass
    multiThread()
    procs = sharedProcs.barrierGet()
    return procs

class Comm:
    """
    An ensemble-safe version of <m socketComm>.Comm. The members of this class
    block non-ensemble-member-0 process from communicating: only the process
    associated with member 0 communicates.
    """
    def __init__(s,comm,esim):
        """
        The contructor takes an existing <m socketComm>.Comm and an
        EnsembleSimulation as arguments.
        """
        if not esim:
            raise Exception("esim argument is not valid")

        s.comm= comm
        s.procNum = comm.procNum
        s.esim=esim

        return

    def barrier(s,timeout=-1):
        from ensembleSharedObj import SharedObj
        sharedProcs = SharedObj()
        if singleThread():
            procs = s.comm.barrier(timeout)
            sharedProcs.set(procs)
            pass
        multiThread()
        procs = sharedProcs.barrierGet()
        return procs
    def collect(s,msg):
        from ensembleSharedObj import SharedObj
        sharedProcs = SharedObj()
        if singleThread():
            data = s.comm.collect(msg)
            sharedProcs.set(data)
            pass
        multiThread()
        data = sharedProcs.barrierGet()
        return data
    def distribute(s,msg):
        from ensembleSharedObj import SharedObj
        sharedProcs = SharedObj()
        if singleThread():
            msg = s.comm.distribute(msg)
            sharedProcs.set(msg)
            pass
        multiThread()
        msg = sharedProcs.barrierGet()
        return msg
    def writeDataTo(s,proc,msg):
        from ensembleSharedObj import SharedObj
        sharedProcs = SharedObj()
        if singleThread():
            s.comm.writeDataTo(proc,msg)
            sharedProcs.set(msg)
            pass
        multiThread()
        msg = sharedProcs.barrierGet()
        return msg
    def readDataFrom(s,proc):
        from ensembleSharedObj import SharedObj
        sharedProcs = SharedObj()
        if singleThread():
            msg=s.comm.readDataFrom(proc)
            sharedProcs.set(msg)
            pass
        multiThread()
        msg = sharedProcs.barrierGet()
        return msg
    def procs(s):
        """
        Return the a sorted list of the current proc numbers. Only to be
        called by proc 0. 
        """
        from ensembleSharedObj import SharedObj
        sharedProcs = SharedObj()
        if singleThread():
            procs = s.comm.procs()
            sharedProcs.set(procs)
            pass
        multiThread()
        procs = sharedProcs.barrierGet()
        return procs
    def info(s,procNum):
        return s.comm.info(procNum)



import simulation
class EnsembleSimulation(simulation.Simulation):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _ensembleSimulation.new_EnsembleSimulation(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _ensembleSimulation.delete_EnsembleSimulation
    __del__ = lambda self: None

    def type(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_type(self, *args, **kwargs)

    def lookupID(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_lookupID(self, *args, **kwargs)

    def atomID(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_atomID(self, *args, **kwargs)

    def barrier(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_barrier(self, *args, **kwargs)

    def shutdown(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_shutdown(self, *args, **kwargs)

    def abort(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_abort(self, *args, **kwargs)

    def numThreads(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_numThreads(self, *args, **kwargs)

    def singleThreaded(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_singleThreaded(self, *args, **kwargs)

    def singleThread(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_singleThread(self, *args, **kwargs)

    def multiThread(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_multiThread(self, *args, **kwargs)

    def size(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_size(self, *args, **kwargs)

    def subSim(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_subSim(self, *args, **kwargs)

    def deleteAtoms_byIndex(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_deleteAtoms_byIndex(self, *args, **kwargs)

    def members(self, *args):
        return _ensembleSimulation.EnsembleSimulation_members(self, *args)

    def member(self, *args):
        return _ensembleSimulation.EnsembleSimulation_member(self, *args)

    def weight(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_weight(self, *args, **kwargs)

    def setWeights(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setWeights(self, *args, **kwargs)
    AVETYPE_AVE = _ensembleSimulation.EnsembleSimulation_AVETYPE_AVE
    AVETYPE_SUM = _ensembleSimulation.EnsembleSimulation_AVETYPE_SUM

    def aveType(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_aveType(self, *args, **kwargs)

    def setAveType(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setAveType(self, *args, **kwargs)

    def kineticEnergy(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_kineticEnergy(self, *args, **kwargs)

    def registerCallbacks(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_registerCallbacks(self, *args, **kwargs)

    def bondPairByID(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_bondPairByID(self, *args, **kwargs)

    def atomString(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_atomString(self, *args, **kwargs)

    def addDependent(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_addDependent(self, *args, **kwargs)

    def removeDependent(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_removeDependent(self, *args, **kwargs)

    def markAsModified(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_markAsModified(self, *args, **kwargs)

    def modifiedID(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_modifiedID(self, *args, **kwargs)

    def select(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_select(self, *args, **kwargs)

    def meanAtomPosArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_meanAtomPosArr(self, *args, **kwargs)

    def atomPosArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_atomPosArr(self, *args, **kwargs)

    def atomVelArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_atomVelArr(self, *args, **kwargs)

    def atomMassArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_atomMassArr(self, *args, **kwargs)

    def atomNameArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_atomNameArr(self, *args, **kwargs)

    def residueNameArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_residueNameArr(self, *args, **kwargs)

    def segmentNameArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_segmentNameArr(self, *args, **kwargs)

    def chemTypeArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_chemTypeArr(self, *args, **kwargs)

    def residueNumArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_residueNumArr(self, *args, **kwargs)

    def setAtomPosArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setAtomPosArr(self, *args, **kwargs)

    def setAtomVelArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setAtomVelArr(self, *args, **kwargs)

    def setAtomMassArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setAtomMassArr(self, *args, **kwargs)

    def setAtomPos(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setAtomPos(self, *args, **kwargs)

    def setAtomVel(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setAtomVel(self, *args, **kwargs)

    def setAtomMass(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setAtomMass(self, *args, **kwargs)

    def setAtomFric(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setAtomFric(self, *args, **kwargs)

    def setAtomCharge(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setAtomCharge(self, *args, **kwargs)

    def setSegmentName(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setSegmentName(self, *args, **kwargs)

    def setResidueName(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setResidueName(self, *args, **kwargs)

    def setResidueNum(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setResidueNum(self, *args, **kwargs)

    def setAtomName(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setAtomName(self, *args, **kwargs)

    def setChemType(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_setChemType(self, *args, **kwargs)

    def atomPos(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_atomPos(self, *args, **kwargs)

    def atomVel(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_atomVel(self, *args, **kwargs)

    def atomMass(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_atomMass(self, *args, **kwargs)

    def atomFric(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_atomFric(self, *args, **kwargs)

    def atomCharge(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_atomCharge(self, *args, **kwargs)

    def segmentName(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_segmentName(self, *args, **kwargs)

    def residueName(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_residueName(self, *args, **kwargs)

    def residueNum(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_residueNum(self, *args, **kwargs)

    def atomName(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_atomName(self, *args, **kwargs)

    def chemType(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_chemType(self, *args, **kwargs)

    def sync(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_sync(self, *args, **kwargs)

    def resize(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_resize(self, *args, **kwargs)

    def sharedString(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_sharedString(self, *args, **kwargs)
    currentSimulation = staticmethod(_ensembleSimulation.EnsembleSimulation_currentSimulation)
    makeCurrent = staticmethod(_ensembleSimulation.EnsembleSimulation_makeCurrent)
    sizeOneSimulation = staticmethod(_ensembleSimulation.EnsembleSimulation_sizeOneSimulation)
    getEnsembleSimulation = staticmethod(_ensembleSimulation.EnsembleSimulation_getEnsembleSimulation)


    #somewhat fragile...
    def __eq__(self,other):
        if other==None: return False
        return self.name()==other.name()
    def __ne__(self,other):
        if other==None: return True
        return self.name()!=other.name()

    def deleteAtoms(self,arg,noSync=False):
        from atomSel import AtomSel
        if str(type(arg)).find('AtomSel')>=0:
          self.deleteAtoms_byIndex( arg.indices(), noSync )
        else:
          self.deleteAtoms_byIndex( AtomSel(arg,self).indices(), noSync )
        import protocol
        protocol.updatePseudoAtoms(self)
        pass

    def sharedObj(self,obj=0):
        from ensembleSharedObj import SharedObj
        return SharedObj(self,obj)
    def collect(self,val):
        from ensembleSharedObj import collect
        return collect(self,val)
    def __del__(self):
      try:
        if self.thisown: 
          import simulation
          simulation.Simulation_deleteSimulation(self)
          self.thisown = False #don't delete a second time!
      except:
        pass


    def help(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_help(self, *args, **kwargs)

class EnsembleSimulationPtr(EnsembleSimulation):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = EnsembleSimulation

EnsembleSimulation_swigregister = _ensembleSimulation.EnsembleSimulation_swigregister
EnsembleSimulation_swigregister(EnsembleSimulation)
cvar = _ensembleSimulation.cvar
maxMessageSize = cvar.maxMessageSize

def EnsembleSimulation_currentSimulation(*args):
    return _ensembleSimulation.EnsembleSimulation_currentSimulation(*args)
EnsembleSimulation_currentSimulation = _ensembleSimulation.EnsembleSimulation_currentSimulation

def EnsembleSimulation_makeCurrent(*args, **kwargs):
    return _ensembleSimulation.EnsembleSimulation_makeCurrent(*args, **kwargs)
EnsembleSimulation_makeCurrent = _ensembleSimulation.EnsembleSimulation_makeCurrent

def EnsembleSimulation_sizeOneSimulation(*args, **kwargs):
    return _ensembleSimulation.EnsembleSimulation_sizeOneSimulation(*args, **kwargs)
EnsembleSimulation_sizeOneSimulation = _ensembleSimulation.EnsembleSimulation_sizeOneSimulation

def EnsembleSimulation_getEnsembleSimulation(*args, **kwargs):
    return _ensembleSimulation.EnsembleSimulation_getEnsembleSimulation(*args, **kwargs)
EnsembleSimulation_getEnsembleSimulation = _ensembleSimulation.EnsembleSimulation_getEnsembleSimulation

class EnsembleSimulation_SharedString(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _ensembleSimulation.new_EnsembleSimulation_SharedString(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _ensembleSimulation.delete_EnsembleSimulation_SharedString
    __del__ = lambda self: None

    def get(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_SharedString_get(self, *args, **kwargs)

    def set(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleSimulation_SharedString_set(self, *args, **kwargs)

class EnsembleSimulation_SharedStringPtr(EnsembleSimulation_SharedString):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = EnsembleSimulation_SharedString

EnsembleSimulation_SharedString_swigregister = _ensembleSimulation.EnsembleSimulation_SharedString_swigregister
EnsembleSimulation_SharedString_swigregister(EnsembleSimulation_SharedString)

class EnsembleMemberSimulation(simulation.Simulation):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined")
    __repr__ = _swig_repr
    __swig_destroy__ = _ensembleSimulation.delete_EnsembleMemberSimulation
    __del__ = lambda self: None

    def type(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_type(self, *args, **kwargs)

    def lookupID(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_lookupID(self, *args, **kwargs)

    def atomID(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_atomID(self, *args, **kwargs)

    def deleteAtoms_byIndex(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_deleteAtoms_byIndex(self, *args, **kwargs)

    def pid(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_pid(self, *args, **kwargs)

    def barrierID(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_barrierID(self, *args, **kwargs)

    def barrierCnt(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_barrierCnt(self, *args, **kwargs)

    def barrierIncr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_barrierIncr(self, *args, **kwargs)

    def sleep(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_sleep(self, *args, **kwargs)

    def wake(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_wake(self, *args, **kwargs)

    def memberIndex(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_memberIndex(self, *args, **kwargs)

    def weight(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_weight(self, *args, **kwargs)

    def ensembleSim(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_ensembleSim(self, *args, **kwargs)

    def subSim(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_subSim(self, *args, **kwargs)

    def resize(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_resize(self, *args, **kwargs)

    def bondPairByID(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_bondPairByID(self, *args, **kwargs)

    def atomPosArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_atomPosArr(self, *args, **kwargs)

    def atomVelArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_atomVelArr(self, *args, **kwargs)

    def atomMassArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_atomMassArr(self, *args, **kwargs)

    def atomNameArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_atomNameArr(self, *args, **kwargs)

    def residueNameArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_residueNameArr(self, *args, **kwargs)

    def segmentNameArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_segmentNameArr(self, *args, **kwargs)

    def chemTypeArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_chemTypeArr(self, *args, **kwargs)

    def residueNumArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_residueNumArr(self, *args, **kwargs)

    def setAtomPosArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setAtomPosArr(self, *args, **kwargs)

    def setAtomVelArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setAtomVelArr(self, *args, **kwargs)

    def setAtomMassArr(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setAtomMassArr(self, *args, **kwargs)

    def setAtomPos(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setAtomPos(self, *args, **kwargs)

    def setAtomVel(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setAtomVel(self, *args, **kwargs)

    def setAtomMass(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setAtomMass(self, *args, **kwargs)

    def setAtomFric(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setAtomFric(self, *args, **kwargs)

    def setAtomCharge(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setAtomCharge(self, *args, **kwargs)

    def setSegmentName(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setSegmentName(self, *args, **kwargs)

    def setResidueName(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setResidueName(self, *args, **kwargs)

    def setResidueNum(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setResidueNum(self, *args, **kwargs)

    def setAtomName(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setAtomName(self, *args, **kwargs)

    def setChemType(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_setChemType(self, *args, **kwargs)

    def atomPos(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_atomPos(self, *args, **kwargs)

    def atomVel(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_atomVel(self, *args, **kwargs)

    def atomMass(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_atomMass(self, *args, **kwargs)

    def atomFric(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_atomFric(self, *args, **kwargs)

    def atomCharge(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_atomCharge(self, *args, **kwargs)

    def segmentName(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_segmentName(self, *args, **kwargs)

    def residueName(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_residueName(self, *args, **kwargs)

    def residueNum(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_residueNum(self, *args, **kwargs)

    def atomName(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_atomName(self, *args, **kwargs)

    def chemType(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_chemType(self, *args, **kwargs)

    def addDependent(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_addDependent(self, *args, **kwargs)

    def removeDependent(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_removeDependent(self, *args, **kwargs)

    def markAsModified(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_markAsModified(self, *args, **kwargs)

    def modifiedID(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_modifiedID(self, *args, **kwargs)

    def select(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_select(self, *args, **kwargs)

    def calcKE(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_calcKE(self, *args, **kwargs)

    def getKE(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleMemberSimulation_getKE(self, *args, **kwargs)

class EnsembleMemberSimulationPtr(EnsembleMemberSimulation):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = EnsembleMemberSimulation

EnsembleMemberSimulation_swigregister = _ensembleSimulation.EnsembleMemberSimulation_swigregister
EnsembleMemberSimulation_swigregister(EnsembleMemberSimulation)

class EnsembleRMSD(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _ensembleSimulation.new_EnsembleRMSD(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def run(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleRMSD_run(self, *args, **kwargs)

    def rmsd(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleRMSD_rmsd(self, *args, **kwargs)

    def byResidue(self, *args, **kwargs):
        return _ensembleSimulation.EnsembleRMSD_byResidue(self, *args, **kwargs)
    __swig_destroy__ = _ensembleSimulation.delete_EnsembleRMSD
    __del__ = lambda self: None

class EnsembleRMSDPtr(EnsembleRMSD):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = EnsembleRMSD

EnsembleRMSD_swigregister = _ensembleSimulation.EnsembleRMSD_swigregister
EnsembleRMSD_swigregister(EnsembleRMSD)


threadSuffix=""
oldConstructor = EnsembleSimulation.__init__
def construct(self, *args, **kwargs):
    oldConstructor(self,*args,**kwargs)
    global threadSuffix
    sim = EnsembleSimulation_currentSimulation()
    if sim:
      threadSuffix=str(sim.member().memberIndex())
    else:
      threadSuffix=""
# disown - do the cleanup in C++
#    if self.size()>1:
    self.thisown=0
    return    
EnsembleSimulation.__init__ = construct

sizeOneSimulation     = EnsembleSimulation_sizeOneSimulation
getEnsembleSimulation = EnsembleSimulation_getEnsembleSimulation
currentSimulation     = EnsembleSimulation_currentSimulation



pyXplorHelp = help


def help(*args):
    return _ensembleSimulation.help(*args)
help = _ensembleSimulation.help

