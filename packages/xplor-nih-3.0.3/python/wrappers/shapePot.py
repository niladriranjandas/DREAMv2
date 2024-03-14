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
        mname = '.'.join((pkg, '_shapePot')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_shapePot')
    _shapePot = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_shapePot', [dirname(__file__)])
        except ImportError:
            import _shapePot
            return _shapePot
        try:
            _mod = imp.load_module('_shapePot', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _shapePot = swig_import_helper()
    del swig_import_helper
else:
    import _shapePot
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
    MOD_SELF = _shapePot.Modified_MOD_SELF
    MOD_SIMULATION = _shapePot.Modified_MOD_SIMULATION

    def __init__(self, *args, **kwargs):
        this = _shapePot.new_Modified(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def set(self, *args, **kwargs):
        return _shapePot.Modified_set(self, *args, **kwargs)

    def clear(self, *args, **kwargs):
        return _shapePot.Modified_clear(self, *args, **kwargs)

    def update(self, *args, **kwargs):
        return _shapePot.Modified_update(self, *args, **kwargs)

    def value(self, *args, **kwargs):
        return _shapePot.Modified_value(self, *args, **kwargs)

    def __call__(self, *args, **kwargs):
        return _shapePot.Modified___call__(self, *args, **kwargs)
    __swig_destroy__ = _shapePot.delete_Modified
    __del__ = lambda self: None

class ModifiedPtr(Modified):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = Modified

Modified_swigregister = _shapePot.Modified_swigregister
Modified_swigregister(Modified)

class ModifiedBase(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    modified = _swig_property(_shapePot.ModifiedBase_modified_get, _shapePot.ModifiedBase_modified_set)
    registeredSimulations = _swig_property(_shapePot.ModifiedBase_registeredSimulations_get, _shapePot.ModifiedBase_registeredSimulations_set)
    __swig_destroy__ = _shapePot.delete_ModifiedBase
    __del__ = lambda self: None

    def registerTo(self, *args, **kwargs):
        return _shapePot.ModifiedBase_registerTo(self, *args, **kwargs)

    def unRegister(self, *args, **kwargs):
        return _shapePot.ModifiedBase_unRegister(self, *args, **kwargs)

    def updateValues(self, *args, **kwargs):
        return _shapePot.ModifiedBase_updateValues(self, *args, **kwargs)

class ModifiedBasePtr(ModifiedBase):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = ModifiedBase

ModifiedBase_swigregister = _shapePot.ModifiedBase_swigregister
ModifiedBase_swigregister(ModifiedBase)

class VarEnsWeights(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    ensWeights = _swig_property(_shapePot.VarEnsWeights_ensWeights_get, _shapePot.VarEnsWeights_ensWeights_set)
    mult = _swig_property(_shapePot.VarEnsWeights_mult_get, _shapePot.VarEnsWeights_mult_set)

    def __init__(self, *args, **kwargs):
        this = _shapePot.new_VarEnsWeights(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _shapePot.delete_VarEnsWeights
    __del__ = lambda self: None

class VarEnsWeightsPtr(VarEnsWeights):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = VarEnsWeights

VarEnsWeights_swigregister = _shapePot.VarEnsWeights_swigregister
VarEnsWeights_swigregister(VarEnsWeights)

class EnsemblePot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _shapePot.delete_EnsemblePot
    __del__ = lambda self: None

    def calcEnergy(self, *args, **kwargs):
        return _shapePot.EnsemblePot_calcEnergy(self, *args, **kwargs)

    def calcEnergyAndDerivs(self, *args, **kwargs):
        return _shapePot.EnsemblePot_calcEnergyAndDerivs(self, *args, **kwargs)

    def energyMaybeDerivs0(self, *args, **kwargs):
        return _shapePot.EnsemblePot_energyMaybeDerivs0(self, *args, **kwargs)

    def energyMaybeDerivs1(self, *args, **kwargs):
        return _shapePot.EnsemblePot_energyMaybeDerivs1(self, *args, **kwargs)

    def energyMaybeDerivs2(self, *args, **kwargs):
        return _shapePot.EnsemblePot_energyMaybeDerivs2(self, *args, **kwargs)

    def energyMaybeDerivs3(self, *args, **kwargs):
        return _shapePot.EnsemblePot_energyMaybeDerivs3(self, *args, **kwargs)

    def energyMaybeDerivs4(self, *args, **kwargs):
        return _shapePot.EnsemblePot_energyMaybeDerivs4(self, *args, **kwargs)

    def energyMaybeDerivsPre(self, *args, **kwargs):
        return _shapePot.EnsemblePot_energyMaybeDerivsPre(self, *args, **kwargs)

    def energyMaybeDerivsPost(self, *args, **kwargs):
        return _shapePot.EnsemblePot_energyMaybeDerivsPost(self, *args, **kwargs)

    def simulation(self, *args):
        return _shapePot.EnsemblePot_simulation(self, *args)

    def ensWeight(self, *args, **kwargs):
        return _shapePot.EnsemblePot_ensWeight(self, *args, **kwargs)

    def ensWeights(self, *args, **kwargs):
        return _shapePot.EnsemblePot_ensWeights(self, *args, **kwargs)

    def setEnsWeights(self, *args, **kwargs):
        return _shapePot.EnsemblePot_setEnsWeights(self, *args, **kwargs)

    def addEnsWeights(self, *args, **kwargs):
        return _shapePot.EnsemblePot_addEnsWeights(self, *args, **kwargs)

    def getEnsWeights(self, *args, **kwargs):
        return _shapePot.EnsemblePot_getEnsWeights(self, *args, **kwargs)

    def clearEnsWeights(self, *args, **kwargs):
        return _shapePot.EnsemblePot_clearEnsWeights(self, *args, **kwargs)

    def updateEnsWeights(self, *args, **kwargs):
        return _shapePot.EnsemblePot_updateEnsWeights(self, *args, **kwargs)

    def useSimEnsWeights(self, *args, **kwargs):
        return _shapePot.EnsemblePot_useSimEnsWeights(self, *args, **kwargs)

    def setUseSimEnsWeights(self, *args, **kwargs):
        return _shapePot.EnsemblePot_setUseSimEnsWeights(self, *args, **kwargs)

    def calcWDerivs(self, *args, **kwargs):
        return _shapePot.EnsemblePot_calcWDerivs(self, *args, **kwargs)

    def setCalcWDerivs(self, *args, **kwargs):
        return _shapePot.EnsemblePot_setCalcWDerivs(self, *args, **kwargs)

    def ensWeightsInfo(self, *args, **kwargs):
        return _shapePot.EnsemblePot_ensWeightsInfo(self, *args, **kwargs)

class EnsemblePotPtr(EnsemblePot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = EnsemblePot

EnsemblePot_swigregister = _shapePot.EnsemblePot_swigregister
EnsemblePot_swigregister(EnsemblePot)

class rc_EnsemblePot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _shapePot.new_rc_EnsemblePot(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _shapePot.delete_rc_EnsemblePot
    __del__ = lambda self: None

class rc_EnsemblePotPtr(rc_EnsemblePot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = rc_EnsemblePot

rc_EnsemblePot_swigregister = _shapePot.rc_EnsemblePot_swigregister
rc_EnsemblePot_swigregister(rc_EnsemblePot)

class ShapePot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _shapePot.new_ShapePot(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def __deref__(self, *args, **kwargs):
        return _shapePot.ShapePot___deref__(self, *args, **kwargs)

    def __ref__(self, *args, **kwargs):
        return _shapePot.ShapePot___ref__(self, *args, **kwargs)

    def registerInstanceData(self, *args, **kwargs):
        return _shapePot.ShapePot_registerInstanceData(self, *args, **kwargs)

    def decrRefCnt(self, *args, **kwargs):
        return _shapePot.ShapePot_decrRefCnt(self, *args, **kwargs)

    def incrRefCnt(self, *args, **kwargs):
        return _shapePot.ShapePot_incrRefCnt(self, *args, **kwargs)

    def refCnt(self, *args, **kwargs):
        return _shapePot.ShapePot_refCnt(self, *args, **kwargs)

    def instanceData(self, *args, **kwargs):
        return _shapePot.ShapePot_instanceData(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _shapePot.ShapePot_help(self, *args, **kwargs)

    __oldinit__=__init__
    def __init__(self, *args):
        self.__oldinit__(*args)
        self.registerInstanceData(self)

    __swig_destroy__ = _shapePot.delete_ShapePot
    __del__ = lambda self: None

    def energyMaybeDerivs0(self, *args, **kwargs):
        return _shapePot.ShapePot_energyMaybeDerivs0(self, *args, **kwargs)

    def energyMaybeDerivs1(self, *args, **kwargs):
        return _shapePot.ShapePot_energyMaybeDerivs1(self, *args, **kwargs)

    def energyMaybeDerivs2(self, *args, **kwargs):
        return _shapePot.ShapePot_energyMaybeDerivs2(self, *args, **kwargs)

    def energyMaybeDerivs3(self, *args, **kwargs):
        return _shapePot.ShapePot_energyMaybeDerivs3(self, *args, **kwargs)

    def rms(self, *args, **kwargs):
        return _shapePot.ShapePot_rms(self, *args, **kwargs)

    def violations(self, *args, **kwargs):
        return _shapePot.ShapePot_violations(self, *args, **kwargs)

    def numRestraints(self, *args, **kwargs):
        return _shapePot.ShapePot_numRestraints(self, *args, **kwargs)

    def simulation(self, *args, **kwargs):
        return _shapePot.ShapePot_simulation(self, *args, **kwargs)
    qc = _swig_property(_shapePot.ShapePot_qc_get, _shapePot.ShapePot_qc_set)
    qc_target = _swig_property(_shapePot.ShapePot_qc_target_get, _shapePot.ShapePot_qc_target_set)
    theta = _swig_property(_shapePot.ShapePot_theta_get, _shapePot.ShapePot_theta_set)
    cosTheta = _swig_property(_shapePot.ShapePot_cosTheta_get, _shapePot.ShapePot_cosTheta_set)
    contrib = _swig_property(_shapePot.ShapePot_contrib_get, _shapePot.ShapePot_contrib_set)
    valuesVectors = _swig_property(_shapePot.ShapePot_valuesVectors_get, _shapePot.ShapePot_valuesVectors_set)

    def getContrib(self, *args, **kwargs):
        return _shapePot.ShapePot_getContrib(self, *args, **kwargs)

    def getTarget(self, *args, **kwargs):
        return _shapePot.ShapePot_getTarget(self, *args, **kwargs)

    def qCenter(self, *args, **kwargs):
        return _shapePot.ShapePot_qCenter(self, *args, **kwargs)

    def qCenterTarget(self, *args, **kwargs):
        return _shapePot.ShapePot_qCenterTarget(self, *args, **kwargs)

    def rotation(self, *args, **kwargs):
        return _shapePot.ShapePot_rotation(self, *args, **kwargs)

    def atomSel(self, *args, **kwargs):
        return _shapePot.ShapePot_atomSel(self, *args, **kwargs)

    def info(self, *args, **kwargs):
        return _shapePot.ShapePot_info(self, *args, **kwargs)

    def showValues(self, *args, **kwargs):
        return _shapePot.ShapePot_showValues(self, *args, **kwargs)

    def showVectors(self, *args, **kwargs):
        return _shapePot.ShapePot_showVectors(self, *args, **kwargs)

    def sizeScale(self, *args, **kwargs):
        return _shapePot.ShapePot_sizeScale(self, *args, **kwargs)

    def setSizeScale(self, *args, **kwargs):
        return _shapePot.ShapePot_setSizeScale(self, *args, **kwargs)

    def orientScale(self, *args, **kwargs):
        return _shapePot.ShapePot_orientScale(self, *args, **kwargs)

    def setOrientScale(self, *args, **kwargs):
        return _shapePot.ShapePot_setOrientScale(self, *args, **kwargs)

    def targetSel(self, *args, **kwargs):
        return _shapePot.ShapePot_targetSel(self, *args, **kwargs)

    def setTargetSel(self, *args, **kwargs):
        return _shapePot.ShapePot_setTargetSel(self, *args, **kwargs)

    def sizeTol(self, *args, **kwargs):
        return _shapePot.ShapePot_sizeTol(self, *args, **kwargs)

    def setSizeTol(self, *args, **kwargs):
        return _shapePot.ShapePot_setSizeTol(self, *args, **kwargs)

    def orientTol(self, *args, **kwargs):
        return _shapePot.ShapePot_orientTol(self, *args, **kwargs)

    def setOrientTol(self, *args, **kwargs):
        return _shapePot.ShapePot_setOrientTol(self, *args, **kwargs)

    def degenerateTol(self, *args, **kwargs):
        return _shapePot.ShapePot_degenerateTol(self, *args, **kwargs)

    def setDegenerateTol(self, *args, **kwargs):
        return _shapePot.ShapePot_setDegenerateTol(self, *args, **kwargs)

    def targetType(self, *args, **kwargs):
        return _shapePot.ShapePot_targetType(self, *args, **kwargs)

    def setTargetType(self, *args, **kwargs):
        return _shapePot.ShapePot_setTargetType(self, *args, **kwargs)

    def sizePotType(self, *args, **kwargs):
        return _shapePot.ShapePot_sizePotType(self, *args, **kwargs)

    def orientPotType(self, *args, **kwargs):
        return _shapePot.ShapePot_orientPotType(self, *args, **kwargs)

    def setSizePotType(self, *args, **kwargs):
        return _shapePot.ShapePot_setSizePotType(self, *args, **kwargs)

    def setOrientPotType(self, *args, **kwargs):
        return _shapePot.ShapePot_setOrientPotType(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _shapePot.ShapePot_help(self, *args, **kwargs)

    def calcEnergy(self, *args, **kwargs):
        return _shapePot.ShapePot_calcEnergy(self, *args, **kwargs)

    def calcEnergyAndDerivs(self, *args, **kwargs):
        return _shapePot.ShapePot_calcEnergyAndDerivs(self, *args, **kwargs)

    def energyMaybeDerivs4(self, *args, **kwargs):
        return _shapePot.ShapePot_energyMaybeDerivs4(self, *args, **kwargs)

    def energyMaybeDerivsPre(self, *args, **kwargs):
        return _shapePot.ShapePot_energyMaybeDerivsPre(self, *args, **kwargs)

    def energyMaybeDerivsPost(self, *args, **kwargs):
        return _shapePot.ShapePot_energyMaybeDerivsPost(self, *args, **kwargs)

    def ensWeight(self, *args, **kwargs):
        return _shapePot.ShapePot_ensWeight(self, *args, **kwargs)

    def ensWeights(self, *args, **kwargs):
        return _shapePot.ShapePot_ensWeights(self, *args, **kwargs)

    def setEnsWeights(self, *args, **kwargs):
        return _shapePot.ShapePot_setEnsWeights(self, *args, **kwargs)

    def addEnsWeights(self, *args, **kwargs):
        return _shapePot.ShapePot_addEnsWeights(self, *args, **kwargs)

    def getEnsWeights(self, *args, **kwargs):
        return _shapePot.ShapePot_getEnsWeights(self, *args, **kwargs)

    def clearEnsWeights(self, *args, **kwargs):
        return _shapePot.ShapePot_clearEnsWeights(self, *args, **kwargs)

    def updateEnsWeights(self, *args, **kwargs):
        return _shapePot.ShapePot_updateEnsWeights(self, *args, **kwargs)

    def useSimEnsWeights(self, *args, **kwargs):
        return _shapePot.ShapePot_useSimEnsWeights(self, *args, **kwargs)

    def setUseSimEnsWeights(self, *args, **kwargs):
        return _shapePot.ShapePot_setUseSimEnsWeights(self, *args, **kwargs)

    def calcWDerivs(self, *args, **kwargs):
        return _shapePot.ShapePot_calcWDerivs(self, *args, **kwargs)

    def setCalcWDerivs(self, *args, **kwargs):
        return _shapePot.ShapePot_setCalcWDerivs(self, *args, **kwargs)

    def ensWeightsInfo(self, *args, **kwargs):
        return _shapePot.ShapePot_ensWeightsInfo(self, *args, **kwargs)

    def potName(self, *args, **kwargs):
        return _shapePot.ShapePot_potName(self, *args, **kwargs)

    def instanceName(self, *args, **kwargs):
        return _shapePot.ShapePot_instanceName(self, *args, **kwargs)

    def resetPotName(self, *args, **kwargs):
        return _shapePot.ShapePot_resetPotName(self, *args, **kwargs)

    def resetInstanceName(self, *args, **kwargs):
        return _shapePot.ShapePot_resetInstanceName(self, *args, **kwargs)

    def scale(self, *args, **kwargs):
        return _shapePot.ShapePot_scale(self, *args, **kwargs)

    def setScale(self, *args, **kwargs):
        return _shapePot.ShapePot_setScale(self, *args, **kwargs)

    def threshold(self, *args, **kwargs):
        return _shapePot.ShapePot_threshold(self, *args, **kwargs)

    def setThreshold(self, *args, **kwargs):
        return _shapePot.ShapePot_setThreshold(self, *args, **kwargs)

    def updateValues(self, *args, **kwargs):
        return _shapePot.ShapePot_updateValues(self, *args, **kwargs)

    def updateDelta(self, *args, **kwargs):
        return _shapePot.ShapePot_updateDelta(self, *args, **kwargs)
    instanceData_ = _swig_property(_shapePot.ShapePot_instanceData__get, _shapePot.ShapePot_instanceData__set)
    instanceDataCreate = _swig_property(_shapePot.ShapePot_instanceDataCreate_get, _shapePot.ShapePot_instanceDataCreate_set)
    instanceDataCleanup = _swig_property(_shapePot.ShapePot_instanceDataCleanup_get, _shapePot.ShapePot_instanceDataCleanup_set)
    modified = _swig_property(_shapePot.ShapePot_modified_get, _shapePot.ShapePot_modified_set)
    registeredSimulations = _swig_property(_shapePot.ShapePot_registeredSimulations_get, _shapePot.ShapePot_registeredSimulations_set)

    def registerTo(self, *args, **kwargs):
        return _shapePot.ShapePot_registerTo(self, *args, **kwargs)

    def unRegister(self, *args, **kwargs):
        return _shapePot.ShapePot_unRegister(self, *args, **kwargs)

class ShapePotPtr(ShapePot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = ShapePot

ShapePot_swigregister = _shapePot.ShapePot_swigregister
ShapePot_swigregister(ShapePot)


realShapePot = ShapePot
def ShapePot(*args):
    from potProxy import PotProxy
    return PotProxy( realShapePot(*args) )

class ShapePot_EigenPair(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    value = _swig_property(_shapePot.ShapePot_EigenPair_value_get, _shapePot.ShapePot_EigenPair_value_set)
    vector = _swig_property(_shapePot.ShapePot_EigenPair_vector_get, _shapePot.ShapePot_EigenPair_vector_set)

    def __init__(self, *args, **kwargs):
        this = _shapePot.new_ShapePot_EigenPair(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _shapePot.delete_ShapePot_EigenPair
    __del__ = lambda self: None

class ShapePot_EigenPairPtr(ShapePot_EigenPair):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = ShapePot_EigenPair

ShapePot_EigenPair_swigregister = _shapePot.ShapePot_EigenPair_swigregister
ShapePot_EigenPair_swigregister(ShapePot_EigenPair)

class ShapePot_LetterClass(EnsemblePot):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    HARMONIC = _shapePot.ShapePot_LetterClass_HARMONIC
    SQUARE = _shapePot.ShapePot_LetterClass_SQUARE
    PairwiseTarget = _shapePot.ShapePot_LetterClass_PairwiseTarget
    AverageTarget = _shapePot.ShapePot_LetterClass_AverageTarget
    MoleculeTarget = _shapePot.ShapePot_LetterClass_MoleculeTarget
    FixedTarget = _shapePot.ShapePot_LetterClass_FixedTarget

    def __init__(self, *args, **kwargs):
        this = _shapePot.new_ShapePot_LetterClass(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _shapePot.delete_ShapePot_LetterClass
    __del__ = lambda self: None

    def energyMaybeDerivs0(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_energyMaybeDerivs0(self, *args, **kwargs)

    def energyMaybeDerivs1(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_energyMaybeDerivs1(self, *args, **kwargs)

    def energyMaybeDerivs2(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_energyMaybeDerivs2(self, *args, **kwargs)

    def energyMaybeDerivs3(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_energyMaybeDerivs3(self, *args, **kwargs)

    def rms(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_rms(self, *args, **kwargs)

    def violations(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_violations(self, *args, **kwargs)

    def numRestraints(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_numRestraints(self, *args, **kwargs)

    def simulation(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_simulation(self, *args, **kwargs)
    qc = _swig_property(_shapePot.ShapePot_LetterClass_qc_get, _shapePot.ShapePot_LetterClass_qc_set)
    qc_target = _swig_property(_shapePot.ShapePot_LetterClass_qc_target_get, _shapePot.ShapePot_LetterClass_qc_target_set)
    theta = _swig_property(_shapePot.ShapePot_LetterClass_theta_get, _shapePot.ShapePot_LetterClass_theta_set)
    cosTheta = _swig_property(_shapePot.ShapePot_LetterClass_cosTheta_get, _shapePot.ShapePot_LetterClass_cosTheta_set)
    contrib = _swig_property(_shapePot.ShapePot_LetterClass_contrib_get, _shapePot.ShapePot_LetterClass_contrib_set)
    valuesVectors = _swig_property(_shapePot.ShapePot_LetterClass_valuesVectors_get, _shapePot.ShapePot_LetterClass_valuesVectors_set)

    def getContrib(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_getContrib(self, *args, **kwargs)

    def getTarget(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_getTarget(self, *args, **kwargs)

    def qCenter(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_qCenter(self, *args, **kwargs)

    def qCenterTarget(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_qCenterTarget(self, *args, **kwargs)

    def rotation(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_rotation(self, *args, **kwargs)

    def atomSel(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_atomSel(self, *args, **kwargs)

    def info(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_info(self, *args, **kwargs)

    def showValues(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_showValues(self, *args, **kwargs)

    def showVectors(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_showVectors(self, *args, **kwargs)

    def sizeScale(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_sizeScale(self, *args, **kwargs)

    def setSizeScale(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_setSizeScale(self, *args, **kwargs)

    def orientScale(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_orientScale(self, *args, **kwargs)

    def setOrientScale(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_setOrientScale(self, *args, **kwargs)

    def targetSel(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_targetSel(self, *args, **kwargs)

    def setTargetSel(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_setTargetSel(self, *args, **kwargs)

    def sizeTol(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_sizeTol(self, *args, **kwargs)

    def setSizeTol(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_setSizeTol(self, *args, **kwargs)

    def orientTol(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_orientTol(self, *args, **kwargs)

    def setOrientTol(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_setOrientTol(self, *args, **kwargs)

    def degenerateTol(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_degenerateTol(self, *args, **kwargs)

    def setDegenerateTol(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_setDegenerateTol(self, *args, **kwargs)

    def targetType(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_targetType(self, *args, **kwargs)

    def setTargetType(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_setTargetType(self, *args, **kwargs)

    def sizePotType(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_sizePotType(self, *args, **kwargs)

    def orientPotType(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_orientPotType(self, *args, **kwargs)

    def setSizePotType(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_setSizePotType(self, *args, **kwargs)

    def setOrientPotType(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_setOrientPotType(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _shapePot.ShapePot_LetterClass_help(self, *args, **kwargs)

class ShapePot_LetterClassPtr(ShapePot_LetterClass):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = ShapePot_LetterClass

ShapePot_LetterClass_swigregister = _shapePot.ShapePot_LetterClass_swigregister
ShapePot_LetterClass_swigregister(ShapePot_LetterClass)


pyXplorHelp = help


def help(*args):
    return _shapePot.help(*args)
help = _shapePot.help


