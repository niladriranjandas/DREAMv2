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
        mname = '.'.join((pkg, '_varTensor')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_varTensor')
    _varTensor = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_varTensor', [dirname(__file__)])
        except ImportError:
            import _varTensor
            return _varTensor
        try:
            _mod = imp.load_module('_varTensor', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _varTensor = swig_import_helper()
    del swig_import_helper
else:
    import _varTensor
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
    MOD_SELF = _varTensor.Modified_MOD_SELF
    MOD_SIMULATION = _varTensor.Modified_MOD_SIMULATION

    def __init__(self, *args, **kwargs):
        this = _varTensor.new_Modified(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def set(self, *args, **kwargs):
        return _varTensor.Modified_set(self, *args, **kwargs)

    def clear(self, *args, **kwargs):
        return _varTensor.Modified_clear(self, *args, **kwargs)

    def update(self, *args, **kwargs):
        return _varTensor.Modified_update(self, *args, **kwargs)

    def value(self, *args, **kwargs):
        return _varTensor.Modified_value(self, *args, **kwargs)

    def __call__(self, *args, **kwargs):
        return _varTensor.Modified___call__(self, *args, **kwargs)
    __swig_destroy__ = _varTensor.delete_Modified
    __del__ = lambda self: None

class ModifiedPtr(Modified):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = Modified

Modified_swigregister = _varTensor.Modified_swigregister
Modified_swigregister(Modified)

class ModifiedBase(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    modified = _swig_property(_varTensor.ModifiedBase_modified_get, _varTensor.ModifiedBase_modified_set)
    registeredSimulations = _swig_property(_varTensor.ModifiedBase_registeredSimulations_get, _varTensor.ModifiedBase_registeredSimulations_set)
    __swig_destroy__ = _varTensor.delete_ModifiedBase
    __del__ = lambda self: None

    def registerTo(self, *args, **kwargs):
        return _varTensor.ModifiedBase_registerTo(self, *args, **kwargs)

    def unRegister(self, *args, **kwargs):
        return _varTensor.ModifiedBase_unRegister(self, *args, **kwargs)

    def updateValues(self, *args, **kwargs):
        return _varTensor.ModifiedBase_updateValues(self, *args, **kwargs)

class ModifiedBasePtr(ModifiedBase):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = ModifiedBase

ModifiedBase_swigregister = _varTensor.ModifiedBase_swigregister
ModifiedBase_swigregister(ModifiedBase)

class VarEnsWeights(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr
    ensWeights = _swig_property(_varTensor.VarEnsWeights_ensWeights_get, _varTensor.VarEnsWeights_ensWeights_set)
    mult = _swig_property(_varTensor.VarEnsWeights_mult_get, _varTensor.VarEnsWeights_mult_set)

    def __init__(self, *args, **kwargs):
        this = _varTensor.new_VarEnsWeights(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _varTensor.delete_VarEnsWeights
    __del__ = lambda self: None

class VarEnsWeightsPtr(VarEnsWeights):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = VarEnsWeights

VarEnsWeights_swigregister = _varTensor.VarEnsWeights_swigregister
VarEnsWeights_swigregister(VarEnsWeights)

class EnsemblePot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')

    def __init__(self, *args, **kwargs):
        raise AttributeError("No constructor defined - class is abstract")
    __repr__ = _swig_repr
    __swig_destroy__ = _varTensor.delete_EnsemblePot
    __del__ = lambda self: None

    def calcEnergy(self, *args, **kwargs):
        return _varTensor.EnsemblePot_calcEnergy(self, *args, **kwargs)

    def calcEnergyAndDerivs(self, *args, **kwargs):
        return _varTensor.EnsemblePot_calcEnergyAndDerivs(self, *args, **kwargs)

    def energyMaybeDerivs0(self, *args, **kwargs):
        return _varTensor.EnsemblePot_energyMaybeDerivs0(self, *args, **kwargs)

    def energyMaybeDerivs1(self, *args, **kwargs):
        return _varTensor.EnsemblePot_energyMaybeDerivs1(self, *args, **kwargs)

    def energyMaybeDerivs2(self, *args, **kwargs):
        return _varTensor.EnsemblePot_energyMaybeDerivs2(self, *args, **kwargs)

    def energyMaybeDerivs3(self, *args, **kwargs):
        return _varTensor.EnsemblePot_energyMaybeDerivs3(self, *args, **kwargs)

    def energyMaybeDerivs4(self, *args, **kwargs):
        return _varTensor.EnsemblePot_energyMaybeDerivs4(self, *args, **kwargs)

    def energyMaybeDerivsPre(self, *args, **kwargs):
        return _varTensor.EnsemblePot_energyMaybeDerivsPre(self, *args, **kwargs)

    def energyMaybeDerivsPost(self, *args, **kwargs):
        return _varTensor.EnsemblePot_energyMaybeDerivsPost(self, *args, **kwargs)

    def simulation(self, *args):
        return _varTensor.EnsemblePot_simulation(self, *args)

    def ensWeight(self, *args, **kwargs):
        return _varTensor.EnsemblePot_ensWeight(self, *args, **kwargs)

    def ensWeights(self, *args, **kwargs):
        return _varTensor.EnsemblePot_ensWeights(self, *args, **kwargs)

    def setEnsWeights(self, *args, **kwargs):
        return _varTensor.EnsemblePot_setEnsWeights(self, *args, **kwargs)

    def addEnsWeights(self, *args, **kwargs):
        return _varTensor.EnsemblePot_addEnsWeights(self, *args, **kwargs)

    def getEnsWeights(self, *args, **kwargs):
        return _varTensor.EnsemblePot_getEnsWeights(self, *args, **kwargs)

    def clearEnsWeights(self, *args, **kwargs):
        return _varTensor.EnsemblePot_clearEnsWeights(self, *args, **kwargs)

    def updateEnsWeights(self, *args, **kwargs):
        return _varTensor.EnsemblePot_updateEnsWeights(self, *args, **kwargs)

    def useSimEnsWeights(self, *args, **kwargs):
        return _varTensor.EnsemblePot_useSimEnsWeights(self, *args, **kwargs)

    def setUseSimEnsWeights(self, *args, **kwargs):
        return _varTensor.EnsemblePot_setUseSimEnsWeights(self, *args, **kwargs)

    def calcWDerivs(self, *args, **kwargs):
        return _varTensor.EnsemblePot_calcWDerivs(self, *args, **kwargs)

    def setCalcWDerivs(self, *args, **kwargs):
        return _varTensor.EnsemblePot_setCalcWDerivs(self, *args, **kwargs)

    def ensWeightsInfo(self, *args, **kwargs):
        return _varTensor.EnsemblePot_ensWeightsInfo(self, *args, **kwargs)

class EnsemblePotPtr(EnsemblePot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = EnsemblePot

EnsemblePot_swigregister = _varTensor.EnsemblePot_swigregister
EnsemblePot_swigregister(EnsemblePot)

class rc_EnsemblePot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _varTensor.new_rc_EnsemblePot(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _varTensor.delete_rc_EnsemblePot
    __del__ = lambda self: None

class rc_EnsemblePotPtr(rc_EnsemblePot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = rc_EnsemblePot

rc_EnsemblePot_swigregister = _varTensor.rc_EnsemblePot_swigregister
rc_EnsemblePot_swigregister(rc_EnsemblePot)

class CDSList_Pot(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __len__(self, *args, **kwargs):
        return _varTensor.CDSList_Pot___len__(self, *args, **kwargs)

    def __init__(self, *args):
        this = _varTensor.new_CDSList_Pot(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def __getitem__(self, *args):
        return _varTensor.CDSList_Pot___getitem__(self, *args)

    def __delitem__(self, *args, **kwargs):
        return _varTensor.CDSList_Pot___delitem__(self, *args, **kwargs)

    def append(self, *args, **kwargs):
        return _varTensor.CDSList_Pot_append(self, *args, **kwargs)

    def remove(self, *args, **kwargs):
        return _varTensor.CDSList_Pot_remove(self, *args, **kwargs)

    def removeAll(self, *args, **kwargs):
        return _varTensor.CDSList_Pot_removeAll(self, *args, **kwargs)

    def __setitem__(self, *args, **kwargs):
        return _varTensor.CDSList_Pot___setitem__(self, *args, **kwargs)

    def __getslice__(self, *args, **kwargs):
        return _varTensor.CDSList_Pot___getslice__(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _varTensor.CDSList_Pot_help(self, *args, **kwargs)
    __swig_destroy__ = _varTensor.delete_CDSList_Pot
    __del__ = lambda self: None

class CDSList_PotPtr(CDSList_Pot):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = CDSList_Pot

CDSList_Pot_swigregister = _varTensor.CDSList_Pot_swigregister
CDSList_Pot_swigregister(CDSList_Pot)

class VarTensor(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _varTensor.new_VarTensor(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def __deref__(self, *args, **kwargs):
        return _varTensor.VarTensor___deref__(self, *args, **kwargs)

    def __ref__(self, *args, **kwargs):
        return _varTensor.VarTensor___ref__(self, *args, **kwargs)

    def registerInstanceData(self, *args, **kwargs):
        return _varTensor.VarTensor_registerInstanceData(self, *args, **kwargs)

    def decrRefCnt(self, *args, **kwargs):
        return _varTensor.VarTensor_decrRefCnt(self, *args, **kwargs)

    def incrRefCnt(self, *args, **kwargs):
        return _varTensor.VarTensor_incrRefCnt(self, *args, **kwargs)

    def refCnt(self, *args, **kwargs):
        return _varTensor.VarTensor_refCnt(self, *args, **kwargs)

    def instanceData(self, *args, **kwargs):
        return _varTensor.VarTensor_instanceData(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _varTensor.VarTensor_help(self, *args, **kwargs)

    __oldinit__=__init__
    def __init__(self, *args):
        self.__oldinit__(*args)
        self.registerInstanceData(self)

    __swig_destroy__ = _varTensor.delete_VarTensor
    __del__ = lambda self: None

    def Da(self, *args, **kwargs):
        return _varTensor.VarTensor_Da(self, *args, **kwargs)

    def aveDa(self, *args, **kwargs):
        return _varTensor.VarTensor_aveDa(self, *args, **kwargs)

    def Rh(self, *args, **kwargs):
        return _varTensor.VarTensor_Rh(self, *args, **kwargs)

    def aveRh(self, *args, **kwargs):
        return _varTensor.VarTensor_aveRh(self, *args, **kwargs)

    def energyMaybeDerivs0(self, *args, **kwargs):
        return _varTensor.VarTensor_energyMaybeDerivs0(self, *args, **kwargs)

    def energyMaybeDerivs1(self, *args, **kwargs):
        return _varTensor.VarTensor_energyMaybeDerivs1(self, *args, **kwargs)

    def rms(self, *args, **kwargs):
        return _varTensor.VarTensor_rms(self, *args, **kwargs)

    def violations(self, *args, **kwargs):
        return _varTensor.VarTensor_violations(self, *args, **kwargs)

    def numRestraints(self, *args, **kwargs):
        return _varTensor.VarTensor_numRestraints(self, *args, **kwargs)

    def setDa(self, *args, **kwargs):
        return _varTensor.VarTensor_setDa(self, *args, **kwargs)

    def setRh(self, *args, **kwargs):
        return _varTensor.VarTensor_setRh(self, *args, **kwargs)
    DaMax_ = _swig_property(_varTensor.VarTensor_DaMax__get, _varTensor.VarTensor_DaMax__set)

    def DaMax(self, *args, **kwargs):
        return _varTensor.VarTensor_DaMax(self, *args, **kwargs)

    def setDaMax(self, *args, **kwargs):
        return _varTensor.VarTensor_setDaMax(self, *args, **kwargs)

    def simulation(self, *args, **kwargs):
        return _varTensor.VarTensor_simulation(self, *args, **kwargs)

    def info(self, *args, **kwargs):
        return _varTensor.VarTensor_info(self, *args, **kwargs)

    def ensembleAxis(self, *args, **kwargs):
        return _varTensor.VarTensor_ensembleAxis(self, *args, **kwargs)

    def ensembleDa(self, *args, **kwargs):
        return _varTensor.VarTensor_ensembleDa(self, *args, **kwargs)

    def ensembleRh(self, *args, **kwargs):
        return _varTensor.VarTensor_ensembleRh(self, *args, **kwargs)

    def setEnsembleAxis(self, *args, **kwargs):
        return _varTensor.VarTensor_setEnsembleAxis(self, *args, **kwargs)

    def setEnsembleDa(self, *args, **kwargs):
        return _varTensor.VarTensor_setEnsembleDa(self, *args, **kwargs)

    def setEnsembleRh(self, *args, **kwargs):
        return _varTensor.VarTensor_setEnsembleRh(self, *args, **kwargs)

    def configTensorAtoms(self, *args, **kwargs):
        return _varTensor.VarTensor_configTensorAtoms(self, *args, **kwargs)

    def scaleDa(self, *args, **kwargs):
        return _varTensor.VarTensor_scaleDa(self, *args, **kwargs)

    def setScaleDa(self, *args, **kwargs):
        return _varTensor.VarTensor_setScaleDa(self, *args, **kwargs)

    def scaleRh(self, *args, **kwargs):
        return _varTensor.VarTensor_scaleRh(self, *args, **kwargs)

    def setScaleRh(self, *args, **kwargs):
        return _varTensor.VarTensor_setScaleRh(self, *args, **kwargs)

    def spreadDa(self, *args, **kwargs):
        return _varTensor.VarTensor_spreadDa(self, *args, **kwargs)

    def setSpreadDa(self, *args, **kwargs):
        return _varTensor.VarTensor_setSpreadDa(self, *args, **kwargs)

    def spreadRh(self, *args, **kwargs):
        return _varTensor.VarTensor_spreadRh(self, *args, **kwargs)

    def setSpreadRh(self, *args, **kwargs):
        return _varTensor.VarTensor_setSpreadRh(self, *args, **kwargs)

    def oAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_oAtom(self, *args, **kwargs)

    def setOAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_setOAtom(self, *args, **kwargs)

    def xAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_xAtom(self, *args, **kwargs)

    def setXAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_setXAtom(self, *args, **kwargs)

    def yAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_yAtom(self, *args, **kwargs)

    def setYAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_setYAtom(self, *args, **kwargs)

    def zAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_zAtom(self, *args, **kwargs)

    def setZAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_setZAtom(self, *args, **kwargs)

    def o2Atom(self, *args, **kwargs):
        return _varTensor.VarTensor_o2Atom(self, *args, **kwargs)

    def setO2Atom(self, *args, **kwargs):
        return _varTensor.VarTensor_setO2Atom(self, *args, **kwargs)

    def p1Atom(self, *args, **kwargs):
        return _varTensor.VarTensor_p1Atom(self, *args, **kwargs)

    def setP1Atom(self, *args, **kwargs):
        return _varTensor.VarTensor_setP1Atom(self, *args, **kwargs)

    def p2Atom(self, *args, **kwargs):
        return _varTensor.VarTensor_p2Atom(self, *args, **kwargs)

    def setP2Atom(self, *args, **kwargs):
        return _varTensor.VarTensor_setP2Atom(self, *args, **kwargs)

    def oAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_oAtomSel(self, *args, **kwargs)

    def setOAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_setOAtomSel(self, *args, **kwargs)

    def xAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_xAtomSel(self, *args, **kwargs)

    def setXAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_setXAtomSel(self, *args, **kwargs)

    def yAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_yAtomSel(self, *args, **kwargs)

    def setYAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_setYAtomSel(self, *args, **kwargs)

    def zAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_zAtomSel(self, *args, **kwargs)

    def setZAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_setZAtomSel(self, *args, **kwargs)

    def o2AtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_o2AtomSel(self, *args, **kwargs)

    def setO2AtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_setO2AtomSel(self, *args, **kwargs)

    def p1AtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_p1AtomSel(self, *args, **kwargs)

    def setP1AtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_setP1AtomSel(self, *args, **kwargs)

    def p2AtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_p2AtomSel(self, *args, **kwargs)

    def setP2AtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_setP2AtomSel(self, *args, **kwargs)

    def freedom(self, *args, **kwargs):
        return _varTensor.VarTensor_freedom(self, *args, **kwargs)

    def setFreedom(self, *args, **kwargs):
        return _varTensor.VarTensor_setFreedom(self, *args, **kwargs)

    def expts(self, *args, **kwargs):
        return _varTensor.VarTensor_expts(self, *args, **kwargs)

    def setExpts(self, *args, **kwargs):
        return _varTensor.VarTensor_setExpts(self, *args, **kwargs)
    dO = _swig_property(_varTensor.VarTensor_dO_get, _varTensor.VarTensor_dO_set)
    dX = _swig_property(_varTensor.VarTensor_dX_get, _varTensor.VarTensor_dX_set)
    dY = _swig_property(_varTensor.VarTensor_dY_get, _varTensor.VarTensor_dY_set)
    dZ = _swig_property(_varTensor.VarTensor_dZ_get, _varTensor.VarTensor_dZ_set)
    dP1 = _swig_property(_varTensor.VarTensor_dP1_get, _varTensor.VarTensor_dP1_set)
    dP2 = _swig_property(_varTensor.VarTensor_dP2_get, _varTensor.VarTensor_dP2_set)

    def initDerivs(self, *args, **kwargs):
        return _varTensor.VarTensor_initDerivs(self, *args, **kwargs)

    def accumDerivs(self, *args, **kwargs):
        return _varTensor.VarTensor_accumDerivs(self, *args, **kwargs)
    U = _swig_property(_varTensor.VarTensor_U_get, _varTensor.VarTensor_U_set)
    normQ = _swig_property(_varTensor.VarTensor_normQ_get, _varTensor.VarTensor_normQ_set)
    xp1hat = _swig_property(_varTensor.VarTensor_xp1hat_get, _varTensor.VarTensor_xp1hat_set)
    zp1hat = _swig_property(_varTensor.VarTensor_zp1hat_get, _varTensor.VarTensor_zp1hat_set)
    zp2hat = _swig_property(_varTensor.VarTensor_zp2hat_get, _varTensor.VarTensor_zp2hat_set)
    qp1 = _swig_property(_varTensor.VarTensor_qp1_get, _varTensor.VarTensor_qp1_set)
    qp2 = _swig_property(_varTensor.VarTensor_qp2_get, _varTensor.VarTensor_qp2_set)
    normQp1 = _swig_property(_varTensor.VarTensor_normQp1_get, _varTensor.VarTensor_normQp1_set)
    normQp2 = _swig_property(_varTensor.VarTensor_normQp2_get, _varTensor.VarTensor_normQp2_set)
    sqrtOneMinusZP1Squared = _swig_property(_varTensor.VarTensor_sqrtOneMinusZP1Squared_get, _varTensor.VarTensor_sqrtOneMinusZP1Squared_set)

    def help(self, *args, **kwargs):
        return _varTensor.VarTensor_help(self, *args, **kwargs)

    def calcEnergy(self, *args, **kwargs):
        return _varTensor.VarTensor_calcEnergy(self, *args, **kwargs)

    def calcEnergyAndDerivs(self, *args, **kwargs):
        return _varTensor.VarTensor_calcEnergyAndDerivs(self, *args, **kwargs)

    def energyMaybeDerivs2(self, *args, **kwargs):
        return _varTensor.VarTensor_energyMaybeDerivs2(self, *args, **kwargs)

    def energyMaybeDerivs3(self, *args, **kwargs):
        return _varTensor.VarTensor_energyMaybeDerivs3(self, *args, **kwargs)

    def energyMaybeDerivs4(self, *args, **kwargs):
        return _varTensor.VarTensor_energyMaybeDerivs4(self, *args, **kwargs)

    def energyMaybeDerivsPre(self, *args, **kwargs):
        return _varTensor.VarTensor_energyMaybeDerivsPre(self, *args, **kwargs)

    def energyMaybeDerivsPost(self, *args, **kwargs):
        return _varTensor.VarTensor_energyMaybeDerivsPost(self, *args, **kwargs)

    def ensWeight(self, *args, **kwargs):
        return _varTensor.VarTensor_ensWeight(self, *args, **kwargs)

    def ensWeights(self, *args, **kwargs):
        return _varTensor.VarTensor_ensWeights(self, *args, **kwargs)

    def setEnsWeights(self, *args, **kwargs):
        return _varTensor.VarTensor_setEnsWeights(self, *args, **kwargs)

    def addEnsWeights(self, *args, **kwargs):
        return _varTensor.VarTensor_addEnsWeights(self, *args, **kwargs)

    def getEnsWeights(self, *args, **kwargs):
        return _varTensor.VarTensor_getEnsWeights(self, *args, **kwargs)

    def clearEnsWeights(self, *args, **kwargs):
        return _varTensor.VarTensor_clearEnsWeights(self, *args, **kwargs)

    def updateEnsWeights(self, *args, **kwargs):
        return _varTensor.VarTensor_updateEnsWeights(self, *args, **kwargs)

    def useSimEnsWeights(self, *args, **kwargs):
        return _varTensor.VarTensor_useSimEnsWeights(self, *args, **kwargs)

    def setUseSimEnsWeights(self, *args, **kwargs):
        return _varTensor.VarTensor_setUseSimEnsWeights(self, *args, **kwargs)

    def calcWDerivs(self, *args, **kwargs):
        return _varTensor.VarTensor_calcWDerivs(self, *args, **kwargs)

    def setCalcWDerivs(self, *args, **kwargs):
        return _varTensor.VarTensor_setCalcWDerivs(self, *args, **kwargs)

    def ensWeightsInfo(self, *args, **kwargs):
        return _varTensor.VarTensor_ensWeightsInfo(self, *args, **kwargs)

    def potName(self, *args, **kwargs):
        return _varTensor.VarTensor_potName(self, *args, **kwargs)

    def instanceName(self, *args, **kwargs):
        return _varTensor.VarTensor_instanceName(self, *args, **kwargs)

    def resetPotName(self, *args, **kwargs):
        return _varTensor.VarTensor_resetPotName(self, *args, **kwargs)

    def resetInstanceName(self, *args, **kwargs):
        return _varTensor.VarTensor_resetInstanceName(self, *args, **kwargs)

    def scale(self, *args, **kwargs):
        return _varTensor.VarTensor_scale(self, *args, **kwargs)

    def setScale(self, *args, **kwargs):
        return _varTensor.VarTensor_setScale(self, *args, **kwargs)

    def threshold(self, *args, **kwargs):
        return _varTensor.VarTensor_threshold(self, *args, **kwargs)

    def setThreshold(self, *args, **kwargs):
        return _varTensor.VarTensor_setThreshold(self, *args, **kwargs)

    def updateValues(self, *args, **kwargs):
        return _varTensor.VarTensor_updateValues(self, *args, **kwargs)

    def updateDelta(self, *args, **kwargs):
        return _varTensor.VarTensor_updateDelta(self, *args, **kwargs)
    instanceData_ = _swig_property(_varTensor.VarTensor_instanceData__get, _varTensor.VarTensor_instanceData__set)
    instanceDataCreate = _swig_property(_varTensor.VarTensor_instanceDataCreate_get, _varTensor.VarTensor_instanceDataCreate_set)
    instanceDataCleanup = _swig_property(_varTensor.VarTensor_instanceDataCleanup_get, _varTensor.VarTensor_instanceDataCleanup_set)
    modified = _swig_property(_varTensor.VarTensor_modified_get, _varTensor.VarTensor_modified_set)
    registeredSimulations = _swig_property(_varTensor.VarTensor_registeredSimulations_get, _varTensor.VarTensor_registeredSimulations_set)

    def registerTo(self, *args, **kwargs):
        return _varTensor.VarTensor_registerTo(self, *args, **kwargs)

    def unRegister(self, *args, **kwargs):
        return _varTensor.VarTensor_unRegister(self, *args, **kwargs)

class VarTensorPtr(VarTensor):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = VarTensor

VarTensor_swigregister = _varTensor.VarTensor_swigregister
VarTensor_swigregister(VarTensor)


realVarTensor = VarTensor
def VarTensor(*args):
    from potProxy import PotProxy
    return PotProxy( realVarTensor(*args) )

class VarTensor_LetterClass(EnsemblePot):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _varTensor.new_VarTensor_LetterClass(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def Da(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_Da(self, *args, **kwargs)

    def aveDa(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_aveDa(self, *args, **kwargs)

    def Rh(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_Rh(self, *args, **kwargs)

    def aveRh(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_aveRh(self, *args, **kwargs)

    def energyMaybeDerivs0(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_energyMaybeDerivs0(self, *args, **kwargs)

    def energyMaybeDerivs1(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_energyMaybeDerivs1(self, *args, **kwargs)

    def rms(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_rms(self, *args, **kwargs)

    def violations(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_violations(self, *args, **kwargs)

    def numRestraints(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_numRestraints(self, *args, **kwargs)

    def setDa(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setDa(self, *args, **kwargs)

    def setRh(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setRh(self, *args, **kwargs)
    DaMax_ = _swig_property(_varTensor.VarTensor_LetterClass_DaMax__get, _varTensor.VarTensor_LetterClass_DaMax__set)

    def DaMax(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_DaMax(self, *args, **kwargs)

    def setDaMax(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setDaMax(self, *args, **kwargs)

    def simulation(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_simulation(self, *args, **kwargs)

    def info(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_info(self, *args, **kwargs)

    def ensembleAxis(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_ensembleAxis(self, *args, **kwargs)

    def ensembleDa(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_ensembleDa(self, *args, **kwargs)

    def ensembleRh(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_ensembleRh(self, *args, **kwargs)

    def setEnsembleAxis(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setEnsembleAxis(self, *args, **kwargs)

    def setEnsembleDa(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setEnsembleDa(self, *args, **kwargs)

    def setEnsembleRh(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setEnsembleRh(self, *args, **kwargs)

    def configTensorAtoms(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_configTensorAtoms(self, *args, **kwargs)

    def scaleDa(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_scaleDa(self, *args, **kwargs)

    def setScaleDa(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setScaleDa(self, *args, **kwargs)

    def scaleRh(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_scaleRh(self, *args, **kwargs)

    def setScaleRh(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setScaleRh(self, *args, **kwargs)

    def spreadDa(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_spreadDa(self, *args, **kwargs)

    def setSpreadDa(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setSpreadDa(self, *args, **kwargs)

    def spreadRh(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_spreadRh(self, *args, **kwargs)

    def setSpreadRh(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setSpreadRh(self, *args, **kwargs)

    def oAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_oAtom(self, *args, **kwargs)

    def setOAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setOAtom(self, *args, **kwargs)

    def xAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_xAtom(self, *args, **kwargs)

    def setXAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setXAtom(self, *args, **kwargs)

    def yAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_yAtom(self, *args, **kwargs)

    def setYAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setYAtom(self, *args, **kwargs)

    def zAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_zAtom(self, *args, **kwargs)

    def setZAtom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setZAtom(self, *args, **kwargs)

    def o2Atom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_o2Atom(self, *args, **kwargs)

    def setO2Atom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setO2Atom(self, *args, **kwargs)

    def p1Atom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_p1Atom(self, *args, **kwargs)

    def setP1Atom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setP1Atom(self, *args, **kwargs)

    def p2Atom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_p2Atom(self, *args, **kwargs)

    def setP2Atom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setP2Atom(self, *args, **kwargs)

    def oAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_oAtomSel(self, *args, **kwargs)

    def setOAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setOAtomSel(self, *args, **kwargs)

    def xAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_xAtomSel(self, *args, **kwargs)

    def setXAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setXAtomSel(self, *args, **kwargs)

    def yAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_yAtomSel(self, *args, **kwargs)

    def setYAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setYAtomSel(self, *args, **kwargs)

    def zAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_zAtomSel(self, *args, **kwargs)

    def setZAtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setZAtomSel(self, *args, **kwargs)

    def o2AtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_o2AtomSel(self, *args, **kwargs)

    def setO2AtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setO2AtomSel(self, *args, **kwargs)

    def p1AtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_p1AtomSel(self, *args, **kwargs)

    def setP1AtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setP1AtomSel(self, *args, **kwargs)

    def p2AtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_p2AtomSel(self, *args, **kwargs)

    def setP2AtomSel(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setP2AtomSel(self, *args, **kwargs)

    def freedom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_freedom(self, *args, **kwargs)

    def setFreedom(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setFreedom(self, *args, **kwargs)

    def expts(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_expts(self, *args, **kwargs)

    def setExpts(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_setExpts(self, *args, **kwargs)
    dO = _swig_property(_varTensor.VarTensor_LetterClass_dO_get, _varTensor.VarTensor_LetterClass_dO_set)
    dX = _swig_property(_varTensor.VarTensor_LetterClass_dX_get, _varTensor.VarTensor_LetterClass_dX_set)
    dY = _swig_property(_varTensor.VarTensor_LetterClass_dY_get, _varTensor.VarTensor_LetterClass_dY_set)
    dZ = _swig_property(_varTensor.VarTensor_LetterClass_dZ_get, _varTensor.VarTensor_LetterClass_dZ_set)
    dP1 = _swig_property(_varTensor.VarTensor_LetterClass_dP1_get, _varTensor.VarTensor_LetterClass_dP1_set)
    dP2 = _swig_property(_varTensor.VarTensor_LetterClass_dP2_get, _varTensor.VarTensor_LetterClass_dP2_set)

    def initDerivs(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_initDerivs(self, *args, **kwargs)

    def accumDerivs(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_accumDerivs(self, *args, **kwargs)
    U = _swig_property(_varTensor.VarTensor_LetterClass_U_get, _varTensor.VarTensor_LetterClass_U_set)
    normQ = _swig_property(_varTensor.VarTensor_LetterClass_normQ_get, _varTensor.VarTensor_LetterClass_normQ_set)
    xp1hat = _swig_property(_varTensor.VarTensor_LetterClass_xp1hat_get, _varTensor.VarTensor_LetterClass_xp1hat_set)
    zp1hat = _swig_property(_varTensor.VarTensor_LetterClass_zp1hat_get, _varTensor.VarTensor_LetterClass_zp1hat_set)
    zp2hat = _swig_property(_varTensor.VarTensor_LetterClass_zp2hat_get, _varTensor.VarTensor_LetterClass_zp2hat_set)
    qp1 = _swig_property(_varTensor.VarTensor_LetterClass_qp1_get, _varTensor.VarTensor_LetterClass_qp1_set)
    qp2 = _swig_property(_varTensor.VarTensor_LetterClass_qp2_get, _varTensor.VarTensor_LetterClass_qp2_set)
    normQp1 = _swig_property(_varTensor.VarTensor_LetterClass_normQp1_get, _varTensor.VarTensor_LetterClass_normQp1_set)
    normQp2 = _swig_property(_varTensor.VarTensor_LetterClass_normQp2_get, _varTensor.VarTensor_LetterClass_normQp2_set)
    sqrtOneMinusZP1Squared = _swig_property(_varTensor.VarTensor_LetterClass_sqrtOneMinusZP1Squared_get, _varTensor.VarTensor_LetterClass_sqrtOneMinusZP1Squared_set)

    def help(self, *args, **kwargs):
        return _varTensor.VarTensor_LetterClass_help(self, *args, **kwargs)
    __swig_destroy__ = _varTensor.delete_VarTensor_LetterClass
    __del__ = lambda self: None

class VarTensor_LetterClassPtr(VarTensor_LetterClass):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = VarTensor_LetterClass

VarTensor_LetterClass_swigregister = _varTensor.VarTensor_LetterClass_swigregister
VarTensor_LetterClass_swigregister(VarTensor_LetterClass)


pyXplorHelp = help


def help(*args):
    return _varTensor.help(*args)
help = _varTensor.help


