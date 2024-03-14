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
        mname = '.'.join((pkg, '_symSimulation')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_symSimulation')
    _symSimulation = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_symSimulation', [dirname(__file__)])
        except ImportError:
            import _symSimulation
            return _symSimulation
        try:
            _mod = imp.load_module('_symSimulation', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _symSimulation = swig_import_helper()
    del swig_import_helper
else:
    import _symSimulation
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
    return _symSimulation.fromSimulation(*args, **kwargs)
fromSimulation = _symSimulation.fromSimulation
import simulation
class Modified(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    MOD_SELF = _symSimulation.Modified_MOD_SELF
    MOD_SIMULATION = _symSimulation.Modified_MOD_SIMULATION

    def __init__(self, *args, **kwargs):
        this = _symSimulation.new_Modified(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def set(self, *args, **kwargs):
        return _symSimulation.Modified_set(self, *args, **kwargs)

    def clear(self, *args, **kwargs):
        return _symSimulation.Modified_clear(self, *args, **kwargs)

    def update(self, *args, **kwargs):
        return _symSimulation.Modified_update(self, *args, **kwargs)

    def value(self, *args, **kwargs):
        return _symSimulation.Modified_value(self, *args, **kwargs)

    def __call__(self, *args, **kwargs):
        return _symSimulation.Modified___call__(self, *args, **kwargs)
    __swig_destroy__ = _symSimulation.delete_Modified
    __del__ = lambda self: None

class ModifiedPtr(Modified):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = Modified

Modified_swigregister = _symSimulation.Modified_swigregister
Modified_swigregister(Modified)

class ModifiedBase(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    modified = _swig_property(_symSimulation.ModifiedBase_modified_get, _symSimulation.ModifiedBase_modified_set)
    registeredSimulations = _swig_property(_symSimulation.ModifiedBase_registeredSimulations_get, _symSimulation.ModifiedBase_registeredSimulations_set)
    __swig_destroy__ = _symSimulation.delete_ModifiedBase
    __del__ = lambda self: None

    def registerTo(self, *args, **kwargs):
        return _symSimulation.ModifiedBase_registerTo(self, *args, **kwargs)

    def unRegister(self, *args, **kwargs):
        return _symSimulation.ModifiedBase_unRegister(self, *args, **kwargs)

    def updateValues(self, *args, **kwargs):
        return _symSimulation.ModifiedBase_updateValues(self, *args, **kwargs)

class ModifiedBasePtr(ModifiedBase):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = ModifiedBase

ModifiedBase_swigregister = _symSimulation.ModifiedBase_swigregister
ModifiedBase_swigregister(ModifiedBase)

class SymSimulationOp(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    rot = _swig_property(_symSimulation.SymSimulationOp_rot_get, _symSimulation.SymSimulationOp_rot_set)
    trans = _swig_property(_symSimulation.SymSimulationOp_trans_get, _symSimulation.SymSimulationOp_trans_set)
    centroidMult = _swig_property(_symSimulation.SymSimulationOp_centroidMult_get, _symSimulation.SymSimulationOp_centroidMult_set)
    centroidProj = _swig_property(_symSimulation.SymSimulationOp_centroidProj_get, _symSimulation.SymSimulationOp_centroidProj_set)
    transAtoms = _swig_property(_symSimulation.SymSimulationOp_transAtoms_get, _symSimulation.SymSimulationOp_transAtoms_set)
    transMult = _swig_property(_symSimulation.SymSimulationOp_transMult_get, _symSimulation.SymSimulationOp_transMult_set)
    rotAtoms = _swig_property(_symSimulation.SymSimulationOp_rotAtoms_get, _symSimulation.SymSimulationOp_rotAtoms_set)
    rotAmount = _swig_property(_symSimulation.SymSimulationOp_rotAmount_get, _symSimulation.SymSimulationOp_rotAmount_set)
    segidPrefix = _swig_property(_symSimulation.SymSimulationOp_segidPrefix_get, _symSimulation.SymSimulationOp_segidPrefix_set)
    segidSuffix = _swig_property(_symSimulation.SymSimulationOp_segidSuffix_get, _symSimulation.SymSimulationOp_segidSuffix_set)

    def __init__(self, *args, **kwargs):
        this = _symSimulation.new_SymSimulationOp(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _symSimulation.delete_SymSimulationOp
    __del__ = lambda self: None

class SymSimulationOpPtr(SymSimulationOp):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = SymSimulationOp

SymSimulationOp_swigregister = _symSimulation.SymSimulationOp_swigregister
SymSimulationOp_swigregister(SymSimulationOp)

class SymSimulation(simulation.Simulation, ModifiedBase):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _symSimulation.new_SymSimulation(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _symSimulation.delete_SymSimulation
    __del__ = lambda self: None

    def updateValues(self, *args, **kwargs):
        return _symSimulation.SymSimulation_updateValues(self, *args, **kwargs)

    def atomID(self, *args, **kwargs):
        return _symSimulation.SymSimulation_atomID(self, *args, **kwargs)

    def type(self, *args, **kwargs):
        return _symSimulation.SymSimulation_type(self, *args, **kwargs)

    def subSel(self, *args, **kwargs):
        return _symSimulation.SymSimulation_subSel(self, *args, **kwargs)

    def addCopy(self, *args, **kwargs):
        return _symSimulation.SymSimulation_addCopy(self, *args, **kwargs)

    def genCopy(self, *args, **kwargs):
        return _symSimulation.SymSimulation_genCopy(self, *args, **kwargs)

    def genCoords(self, *args, **kwargs):
        return _symSimulation.SymSimulation_genCoords(self, *args, **kwargs)

    def genVel(self, *args, **kwargs):
        return _symSimulation.SymSimulation_genVel(self, *args, **kwargs)

    def symOp(self, *args, **kwargs):
        return _symSimulation.SymSimulation_symOp(self, *args, **kwargs)

    def numCopies(self, *args, **kwargs):
        return _symSimulation.SymSimulation_numCopies(self, *args, **kwargs)

    def deleteAtoms(self, *args, **kwargs):
        return _symSimulation.SymSimulation_deleteAtoms(self, *args, **kwargs)

    def atomPosArr(self, *args, **kwargs):
        return _symSimulation.SymSimulation_atomPosArr(self, *args, **kwargs)

    def atomVelArr(self, *args, **kwargs):
        return _symSimulation.SymSimulation_atomVelArr(self, *args, **kwargs)

    def atomMassArr(self, *args, **kwargs):
        return _symSimulation.SymSimulation_atomMassArr(self, *args, **kwargs)

    def atomNameArr(self, *args, **kwargs):
        return _symSimulation.SymSimulation_atomNameArr(self, *args, **kwargs)

    def residueNameArr(self, *args, **kwargs):
        return _symSimulation.SymSimulation_residueNameArr(self, *args, **kwargs)

    def chemTypeArr(self, *args, **kwargs):
        return _symSimulation.SymSimulation_chemTypeArr(self, *args, **kwargs)

    def residueNumArr(self, *args, **kwargs):
        return _symSimulation.SymSimulation_residueNumArr(self, *args, **kwargs)

    def segmentNameArr(self, *args, **kwargs):
        return _symSimulation.SymSimulation_segmentNameArr(self, *args, **kwargs)

    def setAtomPosArr(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setAtomPosArr(self, *args, **kwargs)

    def setAtomVelArr(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setAtomVelArr(self, *args, **kwargs)

    def setAtomMassArr(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setAtomMassArr(self, *args, **kwargs)

    def setAtomPos(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setAtomPos(self, *args, **kwargs)

    def setAtomVel(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setAtomVel(self, *args, **kwargs)

    def setAtomMass(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setAtomMass(self, *args, **kwargs)

    def setAtomFric(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setAtomFric(self, *args, **kwargs)

    def setAtomCharge(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setAtomCharge(self, *args, **kwargs)

    def setSegmentName(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setSegmentName(self, *args, **kwargs)

    def setResidueName(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setResidueName(self, *args, **kwargs)

    def setResidueNum(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setResidueNum(self, *args, **kwargs)

    def setAtomName(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setAtomName(self, *args, **kwargs)

    def setChemType(self, *args, **kwargs):
        return _symSimulation.SymSimulation_setChemType(self, *args, **kwargs)

    def atomPos(self, *args, **kwargs):
        return _symSimulation.SymSimulation_atomPos(self, *args, **kwargs)

    def atomVel(self, *args, **kwargs):
        return _symSimulation.SymSimulation_atomVel(self, *args, **kwargs)

    def atomMass(self, *args, **kwargs):
        return _symSimulation.SymSimulation_atomMass(self, *args, **kwargs)

    def atomFric(self, *args, **kwargs):
        return _symSimulation.SymSimulation_atomFric(self, *args, **kwargs)

    def atomCharge(self, *args, **kwargs):
        return _symSimulation.SymSimulation_atomCharge(self, *args, **kwargs)

    def residueName(self, *args, **kwargs):
        return _symSimulation.SymSimulation_residueName(self, *args, **kwargs)

    def residueNum(self, *args, **kwargs):
        return _symSimulation.SymSimulation_residueNum(self, *args, **kwargs)

    def atomName(self, *args, **kwargs):
        return _symSimulation.SymSimulation_atomName(self, *args, **kwargs)

    def chemType(self, *args, **kwargs):
        return _symSimulation.SymSimulation_chemType(self, *args, **kwargs)

    def segmentName(self, *args, **kwargs):
        return _symSimulation.SymSimulation_segmentName(self, *args, **kwargs)

    def select(self, *args, **kwargs):
        return _symSimulation.SymSimulation_select(self, *args, **kwargs)

    def sync(self, *args, **kwargs):
        return _symSimulation.SymSimulation_sync(self, *args, **kwargs)

    def syncDerivs(self, *args, **kwargs):
        return _symSimulation.SymSimulation_syncDerivs(self, *args, **kwargs)


    #somewhat fragile...
    def __del__(self): #FIX: remove if rc_Simulation is implemented
      try:
        if self.thisown: 
          import simulation
          simulation.Simulation_deleteSimulation(self)
          self.thisown = False #don't delete a second time!
      except:
        pass


    def help(self, *args, **kwargs):
        return _symSimulation.SymSimulation_help(self, *args, **kwargs)

class SymSimulationPtr(SymSimulation):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = SymSimulation

SymSimulation_swigregister = _symSimulation.SymSimulation_swigregister
SymSimulation_swigregister(SymSimulation)


pyXplorHelp = help


def help(*args):
    return _symSimulation.help(*args)
help = _symSimulation.help


