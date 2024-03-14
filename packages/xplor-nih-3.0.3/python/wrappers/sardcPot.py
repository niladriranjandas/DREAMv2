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
        mname = '.'.join((pkg, '_sardcPot')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_sardcPot')
    _sardcPot = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_sardcPot', [dirname(__file__)])
        except ImportError:
            import _sardcPot
            return _sardcPot
        try:
            _mod = imp.load_module('_sardcPot', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _sardcPot = swig_import_helper()
    del swig_import_helper
else:
    import _sardcPot
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


class Modified(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    MOD_SELF = _sardcPot.Modified_MOD_SELF
    MOD_SIMULATION = _sardcPot.Modified_MOD_SIMULATION

    def __init__(self, *args, **kwargs):
        this = _sardcPot.new_Modified(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def set(self, *args, **kwargs):
        return _sardcPot.Modified_set(self, *args, **kwargs)

    def clear(self, *args, **kwargs):
        return _sardcPot.Modified_clear(self, *args, **kwargs)

    def update(self, *args, **kwargs):
        return _sardcPot.Modified_update(self, *args, **kwargs)

    def value(self, *args, **kwargs):
        return _sardcPot.Modified_value(self, *args, **kwargs)

    def __call__(self, *args, **kwargs):
        return _sardcPot.Modified___call__(self, *args, **kwargs)
    __swig_destroy__ = _sardcPot.delete_Modified
    __del__ = lambda self: None

class ModifiedPtr(Modified):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = Modified

Modified_swigregister = _sardcPot.Modified_swigregister
Modified_swigregister(Modified)

class ModifiedBase(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    modified = _swig_property(_sardcPot.ModifiedBase_modified_get, _sardcPot.ModifiedBase_modified_set)
    registeredSimulations = _swig_property(_sardcPot.ModifiedBase_registeredSimulations_get, _sardcPot.ModifiedBase_registeredSimulations_set)
    __swig_destroy__ = _sardcPot.delete_ModifiedBase
    __del__ = lambda self: None

    def registerTo(self, *args, **kwargs):
        return _sardcPot.ModifiedBase_registerTo(self, *args, **kwargs)

    def unRegister(self, *args, **kwargs):
        return _sardcPot.ModifiedBase_unRegister(self, *args, **kwargs)

    def updateValues(self, *args, **kwargs):
        return _sardcPot.ModifiedBase_updateValues(self, *args, **kwargs)

class ModifiedBasePtr(ModifiedBase):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = ModifiedBase

ModifiedBase_swigregister = _sardcPot.ModifiedBase_swigregister
ModifiedBase_swigregister(ModifiedBase)

class VarEnsWeights(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    ensWeights = _swig_property(_sardcPot.VarEnsWeights_ensWeights_get, _sardcPot.VarEnsWeights_ensWeights_set)
    mult = _swig_property(_sardcPot.VarEnsWeights_mult_get, _sardcPot.VarEnsWeights_mult_set)

    def __init__(self, *args, **kwargs):
        this = _sardcPot.new_VarEnsWeights(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _sardcPot.delete_VarEnsWeights
    __del__ = lambda self: None

class VarEnsWeightsPtr(VarEnsWeights):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = VarEnsWeights

VarEnsWeights_swigregister = _sardcPot.VarEnsWeights_swigregister
VarEnsWeights_swigregister(VarEnsWeights)

class EnsemblePot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _sardcPot.delete_EnsemblePot
    __del__ = lambda self: None

    def calcEnergy(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_calcEnergy(self, *args, **kwargs)

    def calcEnergyAndDerivs(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_calcEnergyAndDerivs(self, *args, **kwargs)

    def energyMaybeDerivs0(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_energyMaybeDerivs0(self, *args, **kwargs)

    def energyMaybeDerivs1(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_energyMaybeDerivs1(self, *args, **kwargs)

    def energyMaybeDerivs2(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_energyMaybeDerivs2(self, *args, **kwargs)

    def energyMaybeDerivs3(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_energyMaybeDerivs3(self, *args, **kwargs)

    def energyMaybeDerivs4(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_energyMaybeDerivs4(self, *args, **kwargs)

    def energyMaybeDerivsPre(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_energyMaybeDerivsPre(self, *args, **kwargs)

    def energyMaybeDerivsPost(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_energyMaybeDerivsPost(self, *args, **kwargs)

    def simulation(self, *args):
        return _sardcPot.EnsemblePot_simulation(self, *args)

    def ensWeight(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_ensWeight(self, *args, **kwargs)

    def ensWeights(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_ensWeights(self, *args, **kwargs)

    def setEnsWeights(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_setEnsWeights(self, *args, **kwargs)

    def addEnsWeights(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_addEnsWeights(self, *args, **kwargs)

    def getEnsWeights(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_getEnsWeights(self, *args, **kwargs)

    def clearEnsWeights(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_clearEnsWeights(self, *args, **kwargs)

    def updateEnsWeights(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_updateEnsWeights(self, *args, **kwargs)

    def useSimEnsWeights(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_useSimEnsWeights(self, *args, **kwargs)

    def setUseSimEnsWeights(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_setUseSimEnsWeights(self, *args, **kwargs)

    def calcWDerivs(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_calcWDerivs(self, *args, **kwargs)

    def setCalcWDerivs(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_setCalcWDerivs(self, *args, **kwargs)

    def ensWeightsInfo(self, *args, **kwargs):
        return _sardcPot.EnsemblePot_ensWeightsInfo(self, *args, **kwargs)

class EnsemblePotPtr(EnsemblePot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = EnsemblePot

EnsemblePot_swigregister = _sardcPot.EnsemblePot_swigregister
EnsemblePot_swigregister(EnsemblePot)

class rc_EnsemblePot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _sardcPot.new_rc_EnsemblePot(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _sardcPot.delete_rc_EnsemblePot
    __del__ = lambda self: None

class rc_EnsemblePotPtr(rc_EnsemblePot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = rc_EnsemblePot

rc_EnsemblePot_swigregister = _sardcPot.rc_EnsemblePot_swigregister
rc_EnsemblePot_swigregister(rc_EnsemblePot)

class SARDCPot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _sardcPot.new_SARDCPot(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def __deref__(self, *args, **kwargs):
        return _sardcPot.SARDCPot___deref__(self, *args, **kwargs)

    def __ref__(self, *args, **kwargs):
        return _sardcPot.SARDCPot___ref__(self, *args, **kwargs)

    def registerInstanceData(self, *args, **kwargs):
        return _sardcPot.SARDCPot_registerInstanceData(self, *args, **kwargs)

    def decrRefCnt(self, *args, **kwargs):
        return _sardcPot.SARDCPot_decrRefCnt(self, *args, **kwargs)

    def incrRefCnt(self, *args, **kwargs):
        return _sardcPot.SARDCPot_incrRefCnt(self, *args, **kwargs)

    def refCnt(self, *args, **kwargs):
        return _sardcPot.SARDCPot_refCnt(self, *args, **kwargs)

    def instanceData(self, *args, **kwargs):
        return _sardcPot.SARDCPot_instanceData(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _sardcPot.SARDCPot_help(self, *args, **kwargs)

    __oldinit__=__init__
    def __init__(self, *args):
        self.__oldinit__(*args)
        self.registerInstanceData(self)

    __swig_destroy__ = _sardcPot.delete_SARDCPot
    __del__ = lambda self: None
    tensor = _swig_property(_sardcPot.SARDCPot_tensor_get, _sardcPot.SARDCPot_tensor_set)

    def irredTensor(self, *args, **kwargs):
        return _sardcPot.SARDCPot_irredTensor(self, *args, **kwargs)

    def linkTo(self, *args, **kwargs):
        return _sardcPot.SARDCPot_linkTo(self, *args, **kwargs)

    def addDependent(self, *args, **kwargs):
        return _sardcPot.SARDCPot_addDependent(self, *args, **kwargs)

    def addRestraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_addRestraints(self, *args, **kwargs)

    def energyMaybeDerivs0(self, *args, **kwargs):
        return _sardcPot.SARDCPot_energyMaybeDerivs0(self, *args, **kwargs)

    def energyMaybeDerivs1(self, *args, **kwargs):
        return _sardcPot.SARDCPot_energyMaybeDerivs1(self, *args, **kwargs)

    def energyMaybeDerivs2(self, *args, **kwargs):
        return _sardcPot.SARDCPot_energyMaybeDerivs2(self, *args, **kwargs)

    def energyMaybeDerivs3(self, *args, **kwargs):
        return _sardcPot.SARDCPot_energyMaybeDerivs3(self, *args, **kwargs)

    def energyMaybeDerivs4(self, *args, **kwargs):
        return _sardcPot.SARDCPot_energyMaybeDerivs4(self, *args, **kwargs)

    def rms(self, *args, **kwargs):
        return _sardcPot.SARDCPot_rms(self, *args, **kwargs)

    def chisq(self, *args, **kwargs):
        return _sardcPot.SARDCPot_chisq(self, *args, **kwargs)

    def deviation(self, *args, **kwargs):
        return _sardcPot.SARDCPot_deviation(self, *args, **kwargs)

    def violations(self, *args, **kwargs):
        return _sardcPot.SARDCPot_violations(self, *args, **kwargs)

    def numRestraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_numRestraints(self, *args, **kwargs)

    def restraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_restraints(self, *args, **kwargs)

    def rawRestraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_rawRestraints(self, *args, **kwargs)

    def removeRestraint(self, *args, **kwargs):
        return _sardcPot.SARDCPot_removeRestraint(self, *args, **kwargs)

    def addRestraint(self, *args, **kwargs):
        return _sardcPot.SARDCPot_addRestraint(self, *args, **kwargs)

    def simulation(self, *args, **kwargs):
        return _sardcPot.SARDCPot_simulation(self, *args, **kwargs)

    def info(self, *args, **kwargs):
        return _sardcPot.SARDCPot_info(self, *args, **kwargs)

    def showRestraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_showRestraints(self, *args, **kwargs)

    def showViolations(self, *args, **kwargs):
        return _sardcPot.SARDCPot_showViolations(self, *args, **kwargs)
    avectorScaleUp = _swig_property(_sardcPot.SARDCPot_avectorScaleUp_get, _sardcPot.SARDCPot_avectorScaleUp_set)
    avectorScaleDown = _swig_property(_sardcPot.SARDCPot_avectorScaleDown_get, _sardcPot.SARDCPot_avectorScaleDown_set)

    def avectorScale(self, *args, **kwargs):
        return _sardcPot.SARDCPot_avectorScale(self, *args, **kwargs)

    def setAvectorScale(self, *args, **kwargs):
        return _sardcPot.SARDCPot_setAvectorScale(self, *args, **kwargs)

    def optimizeScale(self, *args, **kwargs):
        return _sardcPot.SARDCPot_optimizeScale(self, *args, **kwargs)

    def setOptimizeScale(self, *args, **kwargs):
        return _sardcPot.SARDCPot_setOptimizeScale(self, *args, **kwargs)

    def omitTensorGrad(self, *args, **kwargs):
        return _sardcPot.SARDCPot_omitTensorGrad(self, *args, **kwargs)

    def setOmitTensorGrad(self, *args, **kwargs):
        return _sardcPot.SARDCPot_setOmitTensorGrad(self, *args, **kwargs)

    def verbose(self, *args, **kwargs):
        return _sardcPot.SARDCPot_verbose(self, *args, **kwargs)

    def setVerbose(self, *args, **kwargs):
        return _sardcPot.SARDCPot_setVerbose(self, *args, **kwargs)

    def showAllRestraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_showAllRestraints(self, *args, **kwargs)

    def setShowAllRestraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_setShowAllRestraints(self, *args, **kwargs)

    def potType(self, *args, **kwargs):
        return _sardcPot.SARDCPot_potType(self, *args, **kwargs)

    def setPotType(self, *args, **kwargs):
        return _sardcPot.SARDCPot_setPotType(self, *args, **kwargs)

    def aveType(self, *args, **kwargs):
        return _sardcPot.SARDCPot_aveType(self, *args, **kwargs)

    def setAveType(self, *args, **kwargs):
        return _sardcPot.SARDCPot_setAveType(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _sardcPot.SARDCPot_help(self, *args, **kwargs)

    def calcEnergy(self, *args, **kwargs):
        return _sardcPot.SARDCPot_calcEnergy(self, *args, **kwargs)

    def calcEnergyAndDerivs(self, *args, **kwargs):
        return _sardcPot.SARDCPot_calcEnergyAndDerivs(self, *args, **kwargs)

    def energyMaybeDerivsPre(self, *args, **kwargs):
        return _sardcPot.SARDCPot_energyMaybeDerivsPre(self, *args, **kwargs)

    def energyMaybeDerivsPost(self, *args, **kwargs):
        return _sardcPot.SARDCPot_energyMaybeDerivsPost(self, *args, **kwargs)

    def ensWeight(self, *args, **kwargs):
        return _sardcPot.SARDCPot_ensWeight(self, *args, **kwargs)

    def ensWeights(self, *args, **kwargs):
        return _sardcPot.SARDCPot_ensWeights(self, *args, **kwargs)

    def setEnsWeights(self, *args, **kwargs):
        return _sardcPot.SARDCPot_setEnsWeights(self, *args, **kwargs)

    def addEnsWeights(self, *args, **kwargs):
        return _sardcPot.SARDCPot_addEnsWeights(self, *args, **kwargs)

    def getEnsWeights(self, *args, **kwargs):
        return _sardcPot.SARDCPot_getEnsWeights(self, *args, **kwargs)

    def clearEnsWeights(self, *args, **kwargs):
        return _sardcPot.SARDCPot_clearEnsWeights(self, *args, **kwargs)

    def updateEnsWeights(self, *args, **kwargs):
        return _sardcPot.SARDCPot_updateEnsWeights(self, *args, **kwargs)

    def useSimEnsWeights(self, *args, **kwargs):
        return _sardcPot.SARDCPot_useSimEnsWeights(self, *args, **kwargs)

    def setUseSimEnsWeights(self, *args, **kwargs):
        return _sardcPot.SARDCPot_setUseSimEnsWeights(self, *args, **kwargs)

    def calcWDerivs(self, *args, **kwargs):
        return _sardcPot.SARDCPot_calcWDerivs(self, *args, **kwargs)

    def setCalcWDerivs(self, *args, **kwargs):
        return _sardcPot.SARDCPot_setCalcWDerivs(self, *args, **kwargs)

    def ensWeightsInfo(self, *args, **kwargs):
        return _sardcPot.SARDCPot_ensWeightsInfo(self, *args, **kwargs)

    def potName(self, *args, **kwargs):
        return _sardcPot.SARDCPot_potName(self, *args, **kwargs)

    def instanceName(self, *args, **kwargs):
        return _sardcPot.SARDCPot_instanceName(self, *args, **kwargs)

    def resetPotName(self, *args, **kwargs):
        return _sardcPot.SARDCPot_resetPotName(self, *args, **kwargs)

    def resetInstanceName(self, *args, **kwargs):
        return _sardcPot.SARDCPot_resetInstanceName(self, *args, **kwargs)

    def scale(self, *args, **kwargs):
        return _sardcPot.SARDCPot_scale(self, *args, **kwargs)

    def setScale(self, *args, **kwargs):
        return _sardcPot.SARDCPot_setScale(self, *args, **kwargs)

    def threshold(self, *args, **kwargs):
        return _sardcPot.SARDCPot_threshold(self, *args, **kwargs)

    def setThreshold(self, *args, **kwargs):
        return _sardcPot.SARDCPot_setThreshold(self, *args, **kwargs)

    def updateValues(self, *args, **kwargs):
        return _sardcPot.SARDCPot_updateValues(self, *args, **kwargs)

    def updateDelta(self, *args, **kwargs):
        return _sardcPot.SARDCPot_updateDelta(self, *args, **kwargs)
    instanceData_ = _swig_property(_sardcPot.SARDCPot_instanceData__get, _sardcPot.SARDCPot_instanceData__set)
    instanceDataCreate = _swig_property(_sardcPot.SARDCPot_instanceDataCreate_get, _sardcPot.SARDCPot_instanceDataCreate_set)
    instanceDataCleanup = _swig_property(_sardcPot.SARDCPot_instanceDataCleanup_get, _sardcPot.SARDCPot_instanceDataCleanup_set)
    modified = _swig_property(_sardcPot.SARDCPot_modified_get, _sardcPot.SARDCPot_modified_set)
    registeredSimulations = _swig_property(_sardcPot.SARDCPot_registeredSimulations_get, _sardcPot.SARDCPot_registeredSimulations_set)

    def registerTo(self, *args, **kwargs):
        return _sardcPot.SARDCPot_registerTo(self, *args, **kwargs)

    def unRegister(self, *args, **kwargs):
        return _sardcPot.SARDCPot_unRegister(self, *args, **kwargs)

class SARDCPotPtr(SARDCPot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = SARDCPot

SARDCPot_swigregister = _sardcPot.SARDCPot_swigregister
SARDCPot_swigregister(SARDCPot)


realSARDCPot = SARDCPot
def SARDCPot(*args):
    from potProxy import PotProxy
    return PotProxy( realSARDCPot(*args) )

class Restraint_SARDCPot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def diff(self, *args, **kwargs):
        return _sardcPot.Restraint_SARDCPot_diff(self, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        this = _sardcPot.new_Restraint_SARDCPot(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def violated(self, *args, **kwargs):
        return _sardcPot.Restraint_SARDCPot_violated(self, *args, **kwargs)

    def name(self, *args, **kwargs):
        return _sardcPot.Restraint_SARDCPot_name(self, *args, **kwargs)

    def setName(self, *args, **kwargs):
        return _sardcPot.Restraint_SARDCPot_setName(self, *args, **kwargs)
    __swig_destroy__ = _sardcPot.delete_Restraint_SARDCPot
    __del__ = lambda self: None

class Restraint_SARDCPotPtr(Restraint_SARDCPot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = Restraint_SARDCPot

Restraint_SARDCPot_swigregister = _sardcPot.Restraint_SARDCPot_swigregister
Restraint_SARDCPot_swigregister(Restraint_SARDCPot)

class SARDCPot_LetterClass(EnsemblePot):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    HARMONIC = _sardcPot.SARDCPot_LetterClass_HARMONIC
    SQUARE = _sardcPot.SARDCPot_LetterClass_SQUARE
    LINEAR = _sardcPot.SARDCPot_LetterClass_LINEAR
    LINEARSQUARE = _sardcPot.SARDCPot_LetterClass_LINEARSQUARE
    AVE = _sardcPot.SARDCPot_LetterClass_AVE
    SUM = _sardcPot.SARDCPot_LetterClass_SUM
    PAIRWISE = _sardcPot.SARDCPot_LetterClass_PAIRWISE

    def __init__(self, *args, **kwargs):
        this = _sardcPot.new_SARDCPot_LetterClass(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _sardcPot.delete_SARDCPot_LetterClass
    __del__ = lambda self: None
    tensor = _swig_property(_sardcPot.SARDCPot_LetterClass_tensor_get, _sardcPot.SARDCPot_LetterClass_tensor_set)

    def irredTensor(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_irredTensor(self, *args, **kwargs)

    def linkTo(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_linkTo(self, *args, **kwargs)

    def addDependent(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_addDependent(self, *args, **kwargs)

    def addRestraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_addRestraints(self, *args, **kwargs)

    def energyMaybeDerivs0(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_energyMaybeDerivs0(self, *args, **kwargs)

    def energyMaybeDerivs1(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_energyMaybeDerivs1(self, *args, **kwargs)

    def energyMaybeDerivs2(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_energyMaybeDerivs2(self, *args, **kwargs)

    def energyMaybeDerivs3(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_energyMaybeDerivs3(self, *args, **kwargs)

    def energyMaybeDerivs4(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_energyMaybeDerivs4(self, *args, **kwargs)

    def rms(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_rms(self, *args, **kwargs)

    def chisq(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_chisq(self, *args, **kwargs)

    def deviation(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_deviation(self, *args, **kwargs)

    def violations(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_violations(self, *args, **kwargs)

    def numRestraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_numRestraints(self, *args, **kwargs)

    def restraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_restraints(self, *args, **kwargs)

    def rawRestraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_rawRestraints(self, *args, **kwargs)

    def removeRestraint(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_removeRestraint(self, *args, **kwargs)

    def addRestraint(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_addRestraint(self, *args, **kwargs)

    def simulation(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_simulation(self, *args, **kwargs)

    def info(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_info(self, *args, **kwargs)

    def showRestraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_showRestraints(self, *args, **kwargs)

    def showViolations(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_showViolations(self, *args, **kwargs)
    avectorScaleUp = _swig_property(_sardcPot.SARDCPot_LetterClass_avectorScaleUp_get, _sardcPot.SARDCPot_LetterClass_avectorScaleUp_set)
    avectorScaleDown = _swig_property(_sardcPot.SARDCPot_LetterClass_avectorScaleDown_get, _sardcPot.SARDCPot_LetterClass_avectorScaleDown_set)

    def avectorScale(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_avectorScale(self, *args, **kwargs)

    def setAvectorScale(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_setAvectorScale(self, *args, **kwargs)

    def optimizeScale(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_optimizeScale(self, *args, **kwargs)

    def setOptimizeScale(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_setOptimizeScale(self, *args, **kwargs)

    def omitTensorGrad(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_omitTensorGrad(self, *args, **kwargs)

    def setOmitTensorGrad(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_setOmitTensorGrad(self, *args, **kwargs)

    def verbose(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_verbose(self, *args, **kwargs)

    def setVerbose(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_setVerbose(self, *args, **kwargs)

    def showAllRestraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_showAllRestraints(self, *args, **kwargs)

    def setShowAllRestraints(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_setShowAllRestraints(self, *args, **kwargs)

    def potType(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_potType(self, *args, **kwargs)

    def setPotType(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_setPotType(self, *args, **kwargs)

    def aveType(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_aveType(self, *args, **kwargs)

    def setAveType(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_setAveType(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _sardcPot.SARDCPot_LetterClass_help(self, *args, **kwargs)

class SARDCPot_LetterClassPtr(SARDCPot_LetterClass):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = SARDCPot_LetterClass

SARDCPot_LetterClass_swigregister = _sardcPot.SARDCPot_LetterClass_swigregister
SARDCPot_LetterClass_swigregister(SARDCPot_LetterClass)

class SARDCPot_Restraint(Restraint_SARDCPot):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    verbose = _swig_property(_sardcPot.SARDCPot_Restraint_verbose_get, _sardcPot.SARDCPot_Restraint_verbose_set)

    def __init__(self, *args, **kwargs):
        this = _sardcPot.new_SARDCPot_Restraint(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _sardcPot.delete_SARDCPot_Restraint
    __del__ = lambda self: None

    def ok(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_ok(self, *args, **kwargs)

    def deriv(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_deriv(self, *args, **kwargs)

    def calcd(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_calcd(self, *args, **kwargs)

    def obs(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_obs(self, *args, **kwargs)

    def setObs(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_setObs(self, *args, **kwargs)

    def setErr(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_setErr(self, *args, **kwargs)

    def aveSize(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_aveSize(self, *args, **kwargs)

    def aSelection(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_aSelection(self, *args, **kwargs)

    def bSelection(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_bSelection(self, *args, **kwargs)

    def err(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_err(self, *args, **kwargs)

    def plusErr(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_plusErr(self, *args, **kwargs)

    def minusErr(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_minusErr(self, *args, **kwargs)

    def name(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_name(self, *args, **kwargs)

    def Dmax(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_Dmax(self, *args, **kwargs)

    def useSign(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_useSign(self, *args, **kwargs)

    def setUseSign(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_setUseSign(self, *args, **kwargs)

    def useDistance(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_useDistance(self, *args, **kwargs)

    def setUseDistance(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_setUseDistance(self, *args, **kwargs)

    def deviation(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_deviation(self, *args, **kwargs)

    def bondVectors(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_bondVectors(self, *args, **kwargs)

    def calcd_contrib(self, *args, **kwargs):
        return _sardcPot.SARDCPot_Restraint_calcd_contrib(self, *args, **kwargs)

class SARDCPot_RestraintPtr(SARDCPot_Restraint):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = SARDCPot_Restraint

SARDCPot_Restraint_swigregister = _sardcPot.SARDCPot_Restraint_swigregister
SARDCPot_Restraint_swigregister(SARDCPot_Restraint)

class rc_ptr_SARDCPot_Restraint(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _sardcPot.new_rc_ptr_SARDCPot_Restraint(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _sardcPot.delete_rc_ptr_SARDCPot_Restraint
    __del__ = lambda self: None

    def __deref__(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint___deref__(self, *args, **kwargs)

    def __ref__(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint___ref__(self, *args, **kwargs)

    def ptr(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_ptr(self, *args, **kwargs)

    def incr(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_incr(self, *args, **kwargs)

    def decr(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_decr(self, *args, **kwargs)

    def count(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_count(self, *args, **kwargs)

    def forceDelete(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_forceDelete(self, *args, **kwargs)

    def reset(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_reset(self, *args, **kwargs)

    def release(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_release(self, *args, **kwargs)
    verbose = _swig_property(_sardcPot.rc_ptr_SARDCPot_Restraint_verbose_get, _sardcPot.rc_ptr_SARDCPot_Restraint_verbose_set)

    def ok(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_ok(self, *args, **kwargs)

    def deriv(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_deriv(self, *args, **kwargs)

    def calcd(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_calcd(self, *args, **kwargs)

    def obs(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_obs(self, *args, **kwargs)

    def setObs(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_setObs(self, *args, **kwargs)

    def setErr(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_setErr(self, *args, **kwargs)

    def aveSize(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_aveSize(self, *args, **kwargs)

    def aSelection(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_aSelection(self, *args, **kwargs)

    def bSelection(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_bSelection(self, *args, **kwargs)

    def err(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_err(self, *args, **kwargs)

    def plusErr(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_plusErr(self, *args, **kwargs)

    def minusErr(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_minusErr(self, *args, **kwargs)

    def name(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_name(self, *args, **kwargs)

    def Dmax(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_Dmax(self, *args, **kwargs)

    def useSign(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_useSign(self, *args, **kwargs)

    def setUseSign(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_setUseSign(self, *args, **kwargs)

    def useDistance(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_useDistance(self, *args, **kwargs)

    def setUseDistance(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_setUseDistance(self, *args, **kwargs)

    def deviation(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_deviation(self, *args, **kwargs)

    def bondVectors(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_bondVectors(self, *args, **kwargs)

    def calcd_contrib(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_calcd_contrib(self, *args, **kwargs)

    def diff(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_diff(self, *args, **kwargs)

    def violated(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_violated(self, *args, **kwargs)

    def setName(self, *args, **kwargs):
        return _sardcPot.rc_ptr_SARDCPot_Restraint_setName(self, *args, **kwargs)

class rc_ptr_SARDCPot_RestraintPtr(rc_ptr_SARDCPot_Restraint):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = rc_ptr_SARDCPot_Restraint

rc_ptr_SARDCPot_Restraint_swigregister = _sardcPot.rc_ptr_SARDCPot_Restraint_swigregister
rc_ptr_SARDCPot_Restraint_swigregister(rc_ptr_SARDCPot_Restraint)

class CDSList_SARDCPot_Restraint(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __len__(self, *args, **kwargs):
        return _sardcPot.CDSList_SARDCPot_Restraint___len__(self, *args, **kwargs)

    def __init__(self, *args):
        this = _sardcPot.new_CDSList_SARDCPot_Restraint(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def __getitem__(self, *args):
        return _sardcPot.CDSList_SARDCPot_Restraint___getitem__(self, *args)

    def __delitem__(self, *args, **kwargs):
        return _sardcPot.CDSList_SARDCPot_Restraint___delitem__(self, *args, **kwargs)

    def append(self, *args, **kwargs):
        return _sardcPot.CDSList_SARDCPot_Restraint_append(self, *args, **kwargs)

    def remove(self, *args, **kwargs):
        return _sardcPot.CDSList_SARDCPot_Restraint_remove(self, *args, **kwargs)

    def removeAll(self, *args, **kwargs):
        return _sardcPot.CDSList_SARDCPot_Restraint_removeAll(self, *args, **kwargs)

    def __setitem__(self, *args, **kwargs):
        return _sardcPot.CDSList_SARDCPot_Restraint___setitem__(self, *args, **kwargs)

    def __getslice__(self, *args, **kwargs):
        return _sardcPot.CDSList_SARDCPot_Restraint___getslice__(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _sardcPot.CDSList_SARDCPot_Restraint_help(self, *args, **kwargs)
    __swig_destroy__ = _sardcPot.delete_CDSList_SARDCPot_Restraint
    __del__ = lambda self: None

class CDSList_SARDCPot_RestraintPtr(CDSList_SARDCPot_Restraint):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = CDSList_SARDCPot_Restraint

CDSList_SARDCPot_Restraint_swigregister = _sardcPot.CDSList_SARDCPot_Restraint_swigregister
CDSList_SARDCPot_Restraint_swigregister(CDSList_SARDCPot_Restraint)


pyXplorHelp = help


def help(*args):
    return _sardcPot.help(*args)
help = _sardcPot.help


