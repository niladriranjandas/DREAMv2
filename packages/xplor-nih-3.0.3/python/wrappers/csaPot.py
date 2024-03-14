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
        mname = '.'.join((pkg, '_csaPot')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_csaPot')
    _csaPot = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_csaPot', [dirname(__file__)])
        except ImportError:
            import _csaPot
            return _csaPot
        try:
            _mod = imp.load_module('_csaPot', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _csaPot = swig_import_helper()
    del swig_import_helper
else:
    import _csaPot
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
    MOD_SELF = _csaPot.Modified_MOD_SELF
    MOD_SIMULATION = _csaPot.Modified_MOD_SIMULATION

    def __init__(self, *args, **kwargs):
        this = _csaPot.new_Modified(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def set(self, *args, **kwargs):
        return _csaPot.Modified_set(self, *args, **kwargs)

    def clear(self, *args, **kwargs):
        return _csaPot.Modified_clear(self, *args, **kwargs)

    def update(self, *args, **kwargs):
        return _csaPot.Modified_update(self, *args, **kwargs)

    def value(self, *args, **kwargs):
        return _csaPot.Modified_value(self, *args, **kwargs)

    def __call__(self, *args, **kwargs):
        return _csaPot.Modified___call__(self, *args, **kwargs)
    __swig_destroy__ = _csaPot.delete_Modified
    __del__ = lambda self: None

class ModifiedPtr(Modified):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = Modified

Modified_swigregister = _csaPot.Modified_swigregister
Modified_swigregister(Modified)

class ModifiedBase(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    modified = _swig_property(_csaPot.ModifiedBase_modified_get, _csaPot.ModifiedBase_modified_set)
    registeredSimulations = _swig_property(_csaPot.ModifiedBase_registeredSimulations_get, _csaPot.ModifiedBase_registeredSimulations_set)
    __swig_destroy__ = _csaPot.delete_ModifiedBase
    __del__ = lambda self: None

    def registerTo(self, *args, **kwargs):
        return _csaPot.ModifiedBase_registerTo(self, *args, **kwargs)

    def unRegister(self, *args, **kwargs):
        return _csaPot.ModifiedBase_unRegister(self, *args, **kwargs)

    def updateValues(self, *args, **kwargs):
        return _csaPot.ModifiedBase_updateValues(self, *args, **kwargs)

class ModifiedBasePtr(ModifiedBase):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = ModifiedBase

ModifiedBase_swigregister = _csaPot.ModifiedBase_swigregister
ModifiedBase_swigregister(ModifiedBase)

class VarEnsWeights(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    ensWeights = _swig_property(_csaPot.VarEnsWeights_ensWeights_get, _csaPot.VarEnsWeights_ensWeights_set)
    mult = _swig_property(_csaPot.VarEnsWeights_mult_get, _csaPot.VarEnsWeights_mult_set)

    def __init__(self, *args, **kwargs):
        this = _csaPot.new_VarEnsWeights(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _csaPot.delete_VarEnsWeights
    __del__ = lambda self: None

class VarEnsWeightsPtr(VarEnsWeights):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = VarEnsWeights

VarEnsWeights_swigregister = _csaPot.VarEnsWeights_swigregister
VarEnsWeights_swigregister(VarEnsWeights)

class EnsemblePot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _csaPot.delete_EnsemblePot
    __del__ = lambda self: None

    def calcEnergy(self, *args, **kwargs):
        return _csaPot.EnsemblePot_calcEnergy(self, *args, **kwargs)

    def calcEnergyAndDerivs(self, *args, **kwargs):
        return _csaPot.EnsemblePot_calcEnergyAndDerivs(self, *args, **kwargs)

    def energyMaybeDerivs0(self, *args, **kwargs):
        return _csaPot.EnsemblePot_energyMaybeDerivs0(self, *args, **kwargs)

    def energyMaybeDerivs1(self, *args, **kwargs):
        return _csaPot.EnsemblePot_energyMaybeDerivs1(self, *args, **kwargs)

    def energyMaybeDerivs2(self, *args, **kwargs):
        return _csaPot.EnsemblePot_energyMaybeDerivs2(self, *args, **kwargs)

    def energyMaybeDerivs3(self, *args, **kwargs):
        return _csaPot.EnsemblePot_energyMaybeDerivs3(self, *args, **kwargs)

    def energyMaybeDerivs4(self, *args, **kwargs):
        return _csaPot.EnsemblePot_energyMaybeDerivs4(self, *args, **kwargs)

    def energyMaybeDerivsPre(self, *args, **kwargs):
        return _csaPot.EnsemblePot_energyMaybeDerivsPre(self, *args, **kwargs)

    def energyMaybeDerivsPost(self, *args, **kwargs):
        return _csaPot.EnsemblePot_energyMaybeDerivsPost(self, *args, **kwargs)

    def simulation(self, *args):
        return _csaPot.EnsemblePot_simulation(self, *args)

    def ensWeight(self, *args, **kwargs):
        return _csaPot.EnsemblePot_ensWeight(self, *args, **kwargs)

    def ensWeights(self, *args, **kwargs):
        return _csaPot.EnsemblePot_ensWeights(self, *args, **kwargs)

    def setEnsWeights(self, *args, **kwargs):
        return _csaPot.EnsemblePot_setEnsWeights(self, *args, **kwargs)

    def addEnsWeights(self, *args, **kwargs):
        return _csaPot.EnsemblePot_addEnsWeights(self, *args, **kwargs)

    def getEnsWeights(self, *args, **kwargs):
        return _csaPot.EnsemblePot_getEnsWeights(self, *args, **kwargs)

    def clearEnsWeights(self, *args, **kwargs):
        return _csaPot.EnsemblePot_clearEnsWeights(self, *args, **kwargs)

    def updateEnsWeights(self, *args, **kwargs):
        return _csaPot.EnsemblePot_updateEnsWeights(self, *args, **kwargs)

    def useSimEnsWeights(self, *args, **kwargs):
        return _csaPot.EnsemblePot_useSimEnsWeights(self, *args, **kwargs)

    def setUseSimEnsWeights(self, *args, **kwargs):
        return _csaPot.EnsemblePot_setUseSimEnsWeights(self, *args, **kwargs)

    def calcWDerivs(self, *args, **kwargs):
        return _csaPot.EnsemblePot_calcWDerivs(self, *args, **kwargs)

    def setCalcWDerivs(self, *args, **kwargs):
        return _csaPot.EnsemblePot_setCalcWDerivs(self, *args, **kwargs)

    def ensWeightsInfo(self, *args, **kwargs):
        return _csaPot.EnsemblePot_ensWeightsInfo(self, *args, **kwargs)

class EnsemblePotPtr(EnsemblePot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = EnsemblePot

EnsemblePot_swigregister = _csaPot.EnsemblePot_swigregister
EnsemblePot_swigregister(EnsemblePot)

class rc_EnsemblePot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _csaPot.new_rc_EnsemblePot(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _csaPot.delete_rc_EnsemblePot
    __del__ = lambda self: None

class rc_EnsemblePotPtr(rc_EnsemblePot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = rc_EnsemblePot

rc_EnsemblePot_swigregister = _csaPot.rc_EnsemblePot_swigregister
rc_EnsemblePot_swigregister(rc_EnsemblePot)

class CSAPot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _csaPot.new_CSAPot(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def __deref__(self, *args, **kwargs):
        return _csaPot.CSAPot___deref__(self, *args, **kwargs)

    def __ref__(self, *args, **kwargs):
        return _csaPot.CSAPot___ref__(self, *args, **kwargs)

    def registerInstanceData(self, *args, **kwargs):
        return _csaPot.CSAPot_registerInstanceData(self, *args, **kwargs)

    def decrRefCnt(self, *args, **kwargs):
        return _csaPot.CSAPot_decrRefCnt(self, *args, **kwargs)

    def incrRefCnt(self, *args, **kwargs):
        return _csaPot.CSAPot_incrRefCnt(self, *args, **kwargs)

    def refCnt(self, *args, **kwargs):
        return _csaPot.CSAPot_refCnt(self, *args, **kwargs)

    def instanceData(self, *args, **kwargs):
        return _csaPot.CSAPot_instanceData(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _csaPot.CSAPot_help(self, *args, **kwargs)

    __oldinit__=__init__
    def __init__(self, *args):
        self.__oldinit__(*args)
        self.registerInstanceData(self)

    __swig_destroy__ = _csaPot.delete_CSAPot
    __del__ = lambda self: None

    def addRestraints(self, *args, **kwargs):
        return _csaPot.CSAPot_addRestraints(self, *args, **kwargs)

    def energyMaybeDerivs0(self, *args, **kwargs):
        return _csaPot.CSAPot_energyMaybeDerivs0(self, *args, **kwargs)

    def energyMaybeDerivs1(self, *args, **kwargs):
        return _csaPot.CSAPot_energyMaybeDerivs1(self, *args, **kwargs)

    def energyMaybeDerivsPost(self, *args, **kwargs):
        return _csaPot.CSAPot_energyMaybeDerivsPost(self, *args, **kwargs)

    def rms(self, *args, **kwargs):
        return _csaPot.CSAPot_rms(self, *args, **kwargs)

    def deviation(self, *args, **kwargs):
        return _csaPot.CSAPot_deviation(self, *args, **kwargs)

    def violations(self, *args, **kwargs):
        return _csaPot.CSAPot_violations(self, *args, **kwargs)

    def numRestraints(self, *args, **kwargs):
        return _csaPot.CSAPot_numRestraints(self, *args, **kwargs)
    oTensor = _swig_property(_csaPot.CSAPot_oTensor_get, _csaPot.CSAPot_oTensor_set)

    def restraints(self, *args, **kwargs):
        return _csaPot.CSAPot_restraints(self, *args, **kwargs)

    def rawRestraints(self, *args, **kwargs):
        return _csaPot.CSAPot_rawRestraints(self, *args, **kwargs)

    def simulation(self, *args, **kwargs):
        return _csaPot.CSAPot_simulation(self, *args, **kwargs)

    def info(self, *args, **kwargs):
        return _csaPot.CSAPot_info(self, *args, **kwargs)

    def showRestraints(self, *args, **kwargs):
        return _csaPot.CSAPot_showRestraints(self, *args, **kwargs)

    def showViolations(self, *args, **kwargs):
        return _csaPot.CSAPot_showViolations(self, *args, **kwargs)

    def verbose(self, *args, **kwargs):
        return _csaPot.CSAPot_verbose(self, *args, **kwargs)

    def setVerbose(self, *args, **kwargs):
        return _csaPot.CSAPot_setVerbose(self, *args, **kwargs)

    def beta(self, *args, **kwargs):
        return _csaPot.CSAPot_beta(self, *args, **kwargs)

    def setBeta(self, *args, **kwargs):
        return _csaPot.CSAPot_setBeta(self, *args, **kwargs)

    def gamma(self, *args, **kwargs):
        return _csaPot.CSAPot_gamma(self, *args, **kwargs)

    def setGamma(self, *args, **kwargs):
        return _csaPot.CSAPot_setGamma(self, *args, **kwargs)

    def sigma(self, *args, **kwargs):
        return _csaPot.CSAPot_sigma(self, *args, **kwargs)

    def setSigma(self, *args, **kwargs):
        return _csaPot.CSAPot_setSigma(self, *args, **kwargs)

    def DaScale(self, *args, **kwargs):
        return _csaPot.CSAPot_DaScale(self, *args, **kwargs)

    def setDaScale(self, *args, **kwargs):
        return _csaPot.CSAPot_setDaScale(self, *args, **kwargs)

    def showAllRestraints(self, *args, **kwargs):
        return _csaPot.CSAPot_showAllRestraints(self, *args, **kwargs)

    def setShowAllRestraints(self, *args, **kwargs):
        return _csaPot.CSAPot_setShowAllRestraints(self, *args, **kwargs)

    def potType(self, *args, **kwargs):
        return _csaPot.CSAPot_potType(self, *args, **kwargs)

    def setPotType(self, *args, **kwargs):
        return _csaPot.CSAPot_setPotType(self, *args, **kwargs)

    def tensorClass(self, *args, **kwargs):
        return _csaPot.CSAPot_tensorClass(self, *args, **kwargs)

    def setTensorClass(self, *args, **kwargs):
        return _csaPot.CSAPot_setTensorClass(self, *args, **kwargs)

    def atomOrder(self, *args, **kwargs):
        return _csaPot.CSAPot_atomOrder(self, *args, **kwargs)

    def setAtomOrder(self, *args, **kwargs):
        return _csaPot.CSAPot_setAtomOrder(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _csaPot.CSAPot_help(self, *args, **kwargs)

    def calcEnergy(self, *args, **kwargs):
        return _csaPot.CSAPot_calcEnergy(self, *args, **kwargs)

    def calcEnergyAndDerivs(self, *args, **kwargs):
        return _csaPot.CSAPot_calcEnergyAndDerivs(self, *args, **kwargs)

    def energyMaybeDerivs2(self, *args, **kwargs):
        return _csaPot.CSAPot_energyMaybeDerivs2(self, *args, **kwargs)

    def energyMaybeDerivs3(self, *args, **kwargs):
        return _csaPot.CSAPot_energyMaybeDerivs3(self, *args, **kwargs)

    def energyMaybeDerivs4(self, *args, **kwargs):
        return _csaPot.CSAPot_energyMaybeDerivs4(self, *args, **kwargs)

    def energyMaybeDerivsPre(self, *args, **kwargs):
        return _csaPot.CSAPot_energyMaybeDerivsPre(self, *args, **kwargs)

    def ensWeight(self, *args, **kwargs):
        return _csaPot.CSAPot_ensWeight(self, *args, **kwargs)

    def ensWeights(self, *args, **kwargs):
        return _csaPot.CSAPot_ensWeights(self, *args, **kwargs)

    def setEnsWeights(self, *args, **kwargs):
        return _csaPot.CSAPot_setEnsWeights(self, *args, **kwargs)

    def addEnsWeights(self, *args, **kwargs):
        return _csaPot.CSAPot_addEnsWeights(self, *args, **kwargs)

    def getEnsWeights(self, *args, **kwargs):
        return _csaPot.CSAPot_getEnsWeights(self, *args, **kwargs)

    def clearEnsWeights(self, *args, **kwargs):
        return _csaPot.CSAPot_clearEnsWeights(self, *args, **kwargs)

    def updateEnsWeights(self, *args, **kwargs):
        return _csaPot.CSAPot_updateEnsWeights(self, *args, **kwargs)

    def useSimEnsWeights(self, *args, **kwargs):
        return _csaPot.CSAPot_useSimEnsWeights(self, *args, **kwargs)

    def setUseSimEnsWeights(self, *args, **kwargs):
        return _csaPot.CSAPot_setUseSimEnsWeights(self, *args, **kwargs)

    def calcWDerivs(self, *args, **kwargs):
        return _csaPot.CSAPot_calcWDerivs(self, *args, **kwargs)

    def setCalcWDerivs(self, *args, **kwargs):
        return _csaPot.CSAPot_setCalcWDerivs(self, *args, **kwargs)

    def ensWeightsInfo(self, *args, **kwargs):
        return _csaPot.CSAPot_ensWeightsInfo(self, *args, **kwargs)

    def potName(self, *args, **kwargs):
        return _csaPot.CSAPot_potName(self, *args, **kwargs)

    def instanceName(self, *args, **kwargs):
        return _csaPot.CSAPot_instanceName(self, *args, **kwargs)

    def resetPotName(self, *args, **kwargs):
        return _csaPot.CSAPot_resetPotName(self, *args, **kwargs)

    def resetInstanceName(self, *args, **kwargs):
        return _csaPot.CSAPot_resetInstanceName(self, *args, **kwargs)

    def scale(self, *args, **kwargs):
        return _csaPot.CSAPot_scale(self, *args, **kwargs)

    def setScale(self, *args, **kwargs):
        return _csaPot.CSAPot_setScale(self, *args, **kwargs)

    def threshold(self, *args, **kwargs):
        return _csaPot.CSAPot_threshold(self, *args, **kwargs)

    def setThreshold(self, *args, **kwargs):
        return _csaPot.CSAPot_setThreshold(self, *args, **kwargs)

    def updateValues(self, *args, **kwargs):
        return _csaPot.CSAPot_updateValues(self, *args, **kwargs)

    def updateDelta(self, *args, **kwargs):
        return _csaPot.CSAPot_updateDelta(self, *args, **kwargs)
    instanceData_ = _swig_property(_csaPot.CSAPot_instanceData__get, _csaPot.CSAPot_instanceData__set)
    instanceDataCreate = _swig_property(_csaPot.CSAPot_instanceDataCreate_get, _csaPot.CSAPot_instanceDataCreate_set)
    instanceDataCleanup = _swig_property(_csaPot.CSAPot_instanceDataCleanup_get, _csaPot.CSAPot_instanceDataCleanup_set)
    modified = _swig_property(_csaPot.CSAPot_modified_get, _csaPot.CSAPot_modified_set)
    registeredSimulations = _swig_property(_csaPot.CSAPot_registeredSimulations_get, _csaPot.CSAPot_registeredSimulations_set)

    def registerTo(self, *args, **kwargs):
        return _csaPot.CSAPot_registerTo(self, *args, **kwargs)

    def unRegister(self, *args, **kwargs):
        return _csaPot.CSAPot_unRegister(self, *args, **kwargs)

class CSAPotPtr(CSAPot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = CSAPot

CSAPot_swigregister = _csaPot.CSAPot_swigregister
CSAPot_swigregister(CSAPot)


realCSAPot = CSAPot
def CSAPot(*args):
    from potProxy import PotProxy
    return PotProxy( realCSAPot(*args) )

class Restraint_CSAPot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def diff(self, *args, **kwargs):
        return _csaPot.Restraint_CSAPot_diff(self, *args, **kwargs)

    def __init__(self, *args, **kwargs):
        this = _csaPot.new_Restraint_CSAPot(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def violated(self, *args, **kwargs):
        return _csaPot.Restraint_CSAPot_violated(self, *args, **kwargs)

    def name(self, *args, **kwargs):
        return _csaPot.Restraint_CSAPot_name(self, *args, **kwargs)

    def setName(self, *args, **kwargs):
        return _csaPot.Restraint_CSAPot_setName(self, *args, **kwargs)
    __swig_destroy__ = _csaPot.delete_Restraint_CSAPot
    __del__ = lambda self: None

class Restraint_CSAPotPtr(Restraint_CSAPot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = Restraint_CSAPot

Restraint_CSAPot_swigregister = _csaPot.Restraint_CSAPot_swigregister
Restraint_CSAPot_swigregister(Restraint_CSAPot)

class CSAPot_LetterClass(EnsemblePot):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    HARMONIC = _csaPot.CSAPot_LetterClass_HARMONIC
    SQUARE = _csaPot.CSAPot_LetterClass_SQUARE
    LINEAR = _csaPot.CSAPot_LetterClass_LINEAR
    LINEARSQUARE = _csaPot.CSAPot_LetterClass_LINEARSQUARE
    BOND = _csaPot.CSAPot_LetterClass_BOND
    BISECT = _csaPot.CSAPot_LetterClass_BISECT
    ORDER_123 = _csaPot.CSAPot_LetterClass_ORDER_123
    ORDER_132 = _csaPot.CSAPot_LetterClass_ORDER_132
    ORDER_231 = _csaPot.CSAPot_LetterClass_ORDER_231

    def __init__(self, *args, **kwargs):
        this = _csaPot.new_CSAPot_LetterClass(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _csaPot.delete_CSAPot_LetterClass
    __del__ = lambda self: None

    def addRestraints(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_addRestraints(self, *args, **kwargs)

    def energyMaybeDerivs0(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_energyMaybeDerivs0(self, *args, **kwargs)

    def energyMaybeDerivs1(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_energyMaybeDerivs1(self, *args, **kwargs)

    def energyMaybeDerivsPost(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_energyMaybeDerivsPost(self, *args, **kwargs)

    def rms(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_rms(self, *args, **kwargs)

    def deviation(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_deviation(self, *args, **kwargs)

    def violations(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_violations(self, *args, **kwargs)

    def numRestraints(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_numRestraints(self, *args, **kwargs)
    oTensor = _swig_property(_csaPot.CSAPot_LetterClass_oTensor_get, _csaPot.CSAPot_LetterClass_oTensor_set)

    def restraints(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_restraints(self, *args, **kwargs)

    def rawRestraints(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_rawRestraints(self, *args, **kwargs)

    def simulation(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_simulation(self, *args, **kwargs)

    def info(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_info(self, *args, **kwargs)

    def showRestraints(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_showRestraints(self, *args, **kwargs)

    def showViolations(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_showViolations(self, *args, **kwargs)

    def verbose(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_verbose(self, *args, **kwargs)

    def setVerbose(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_setVerbose(self, *args, **kwargs)

    def beta(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_beta(self, *args, **kwargs)

    def setBeta(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_setBeta(self, *args, **kwargs)

    def gamma(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_gamma(self, *args, **kwargs)

    def setGamma(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_setGamma(self, *args, **kwargs)

    def sigma(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_sigma(self, *args, **kwargs)

    def setSigma(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_setSigma(self, *args, **kwargs)

    def DaScale(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_DaScale(self, *args, **kwargs)

    def setDaScale(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_setDaScale(self, *args, **kwargs)

    def showAllRestraints(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_showAllRestraints(self, *args, **kwargs)

    def setShowAllRestraints(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_setShowAllRestraints(self, *args, **kwargs)

    def potType(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_potType(self, *args, **kwargs)

    def setPotType(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_setPotType(self, *args, **kwargs)

    def tensorClass(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_tensorClass(self, *args, **kwargs)

    def setTensorClass(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_setTensorClass(self, *args, **kwargs)

    def atomOrder(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_atomOrder(self, *args, **kwargs)

    def setAtomOrder(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_setAtomOrder(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _csaPot.CSAPot_LetterClass_help(self, *args, **kwargs)

class CSAPot_LetterClassPtr(CSAPot_LetterClass):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = CSAPot_LetterClass

CSAPot_LetterClass_swigregister = _csaPot.CSAPot_LetterClass_swigregister
CSAPot_LetterClass_swigregister(CSAPot_LetterClass)

class CSAPot_Restraint(Restraint_CSAPot):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    verbose = _swig_property(_csaPot.CSAPot_Restraint_verbose_get, _csaPot.CSAPot_Restraint_verbose_set)

    def __init__(self, *args, **kwargs):
        this = _csaPot.new_CSAPot_Restraint(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _csaPot.delete_CSAPot_Restraint
    __del__ = lambda self: None

    def ok(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_ok(self, *args, **kwargs)

    def deriv(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_deriv(self, *args, **kwargs)

    def calcd(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_calcd(self, *args, **kwargs)

    def obs(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_obs(self, *args, **kwargs)

    def setObs(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_setObs(self, *args, **kwargs)

    def Selection1(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_Selection1(self, *args, **kwargs)

    def Selection2(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_Selection2(self, *args, **kwargs)

    def Selection3(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_Selection3(self, *args, **kwargs)

    def plusErr(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_plusErr(self, *args, **kwargs)

    def minusErr(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_minusErr(self, *args, **kwargs)

    def name(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_name(self, *args, **kwargs)

    def tensor(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_tensor(self, *args, **kwargs)

    def eigenMatrix(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_eigenMatrix(self, *args, **kwargs)

    def deviation(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_deviation(self, *args, **kwargs)

    def calcd_contrib(self, *args, **kwargs):
        return _csaPot.CSAPot_Restraint_calcd_contrib(self, *args, **kwargs)

class CSAPot_RestraintPtr(CSAPot_Restraint):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = CSAPot_Restraint

CSAPot_Restraint_swigregister = _csaPot.CSAPot_Restraint_swigregister
CSAPot_Restraint_swigregister(CSAPot_Restraint)

class rc_ptr_CSAPot_Restraint(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _csaPot.new_rc_ptr_CSAPot_Restraint(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _csaPot.delete_rc_ptr_CSAPot_Restraint
    __del__ = lambda self: None

    def __deref__(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint___deref__(self, *args, **kwargs)

    def __ref__(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint___ref__(self, *args, **kwargs)

    def ptr(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_ptr(self, *args, **kwargs)

    def incr(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_incr(self, *args, **kwargs)

    def decr(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_decr(self, *args, **kwargs)

    def count(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_count(self, *args, **kwargs)

    def forceDelete(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_forceDelete(self, *args, **kwargs)

    def reset(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_reset(self, *args, **kwargs)

    def release(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_release(self, *args, **kwargs)
    verbose = _swig_property(_csaPot.rc_ptr_CSAPot_Restraint_verbose_get, _csaPot.rc_ptr_CSAPot_Restraint_verbose_set)

    def ok(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_ok(self, *args, **kwargs)

    def deriv(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_deriv(self, *args, **kwargs)

    def calcd(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_calcd(self, *args, **kwargs)

    def obs(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_obs(self, *args, **kwargs)

    def setObs(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_setObs(self, *args, **kwargs)

    def Selection1(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_Selection1(self, *args, **kwargs)

    def Selection2(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_Selection2(self, *args, **kwargs)

    def Selection3(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_Selection3(self, *args, **kwargs)

    def plusErr(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_plusErr(self, *args, **kwargs)

    def minusErr(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_minusErr(self, *args, **kwargs)

    def name(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_name(self, *args, **kwargs)

    def tensor(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_tensor(self, *args, **kwargs)

    def eigenMatrix(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_eigenMatrix(self, *args, **kwargs)

    def deviation(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_deviation(self, *args, **kwargs)

    def calcd_contrib(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_calcd_contrib(self, *args, **kwargs)

    def diff(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_diff(self, *args, **kwargs)

    def violated(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_violated(self, *args, **kwargs)

    def setName(self, *args, **kwargs):
        return _csaPot.rc_ptr_CSAPot_Restraint_setName(self, *args, **kwargs)

class rc_ptr_CSAPot_RestraintPtr(rc_ptr_CSAPot_Restraint):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = rc_ptr_CSAPot_Restraint

rc_ptr_CSAPot_Restraint_swigregister = _csaPot.rc_ptr_CSAPot_Restraint_swigregister
rc_ptr_CSAPot_Restraint_swigregister(rc_ptr_CSAPot_Restraint)

class CDSList_CSAPot_Restraint(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __len__(self, *args, **kwargs):
        return _csaPot.CDSList_CSAPot_Restraint___len__(self, *args, **kwargs)

    def __init__(self, *args):
        this = _csaPot.new_CDSList_CSAPot_Restraint(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def __getitem__(self, *args):
        return _csaPot.CDSList_CSAPot_Restraint___getitem__(self, *args)

    def __delitem__(self, *args, **kwargs):
        return _csaPot.CDSList_CSAPot_Restraint___delitem__(self, *args, **kwargs)

    def append(self, *args, **kwargs):
        return _csaPot.CDSList_CSAPot_Restraint_append(self, *args, **kwargs)

    def remove(self, *args, **kwargs):
        return _csaPot.CDSList_CSAPot_Restraint_remove(self, *args, **kwargs)

    def removeAll(self, *args, **kwargs):
        return _csaPot.CDSList_CSAPot_Restraint_removeAll(self, *args, **kwargs)

    def __setitem__(self, *args, **kwargs):
        return _csaPot.CDSList_CSAPot_Restraint___setitem__(self, *args, **kwargs)

    def __getslice__(self, *args, **kwargs):
        return _csaPot.CDSList_CSAPot_Restraint___getslice__(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _csaPot.CDSList_CSAPot_Restraint_help(self, *args, **kwargs)
    __swig_destroy__ = _csaPot.delete_CDSList_CSAPot_Restraint
    __del__ = lambda self: None

class CDSList_CSAPot_RestraintPtr(CDSList_CSAPot_Restraint):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = CDSList_CSAPot_Restraint

CDSList_CSAPot_Restraint_swigregister = _csaPot.CDSList_CSAPot_Restraint_swigregister
CDSList_CSAPot_Restraint_swigregister(CDSList_CSAPot_Restraint)


pyXplorHelp = help


def help(*args):
    return _csaPot.help(*args)
help = _csaPot.help
class CDSList_FloatPair(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __len__(self, *args, **kwargs):
        return _csaPot.CDSList_FloatPair___len__(self, *args, **kwargs)

    def __init__(self, *args):
        this = _csaPot.new_CDSList_FloatPair(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def __getitem__(self, *args):
        return _csaPot.CDSList_FloatPair___getitem__(self, *args)

    def __delitem__(self, *args, **kwargs):
        return _csaPot.CDSList_FloatPair___delitem__(self, *args, **kwargs)

    def append(self, *args, **kwargs):
        return _csaPot.CDSList_FloatPair_append(self, *args, **kwargs)

    def remove(self, *args, **kwargs):
        return _csaPot.CDSList_FloatPair_remove(self, *args, **kwargs)

    def removeAll(self, *args, **kwargs):
        return _csaPot.CDSList_FloatPair_removeAll(self, *args, **kwargs)

    def __setitem__(self, *args, **kwargs):
        return _csaPot.CDSList_FloatPair___setitem__(self, *args, **kwargs)

    def __getslice__(self, *args, **kwargs):
        return _csaPot.CDSList_FloatPair___getslice__(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _csaPot.CDSList_FloatPair_help(self, *args, **kwargs)
    __swig_destroy__ = _csaPot.delete_CDSList_FloatPair
    __del__ = lambda self: None

class CDSList_FloatPairPtr(CDSList_FloatPair):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = CDSList_FloatPair

CDSList_FloatPair_swigregister = _csaPot.CDSList_FloatPair_swigregister
CDSList_FloatPair_swigregister(CDSList_FloatPair)


def powderPattern(*args, **kwargs):
    return _csaPot.powderPattern(*args, **kwargs)
powderPattern = _csaPot.powderPattern

def convPowderPattern(*args, **kwargs):
    return _csaPot.convPowderPattern(*args, **kwargs)
convPowderPattern = _csaPot.convPowderPattern


