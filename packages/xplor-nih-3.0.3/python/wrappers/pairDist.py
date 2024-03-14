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
        mname = '.'.join((pkg, '_pairDist')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_pairDist')
    _pairDist = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_pairDist', [dirname(__file__)])
        except ImportError:
            import _pairDist
            return _pairDist
        try:
            _mod = imp.load_module('_pairDist', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _pairDist = swig_import_helper()
    del swig_import_helper
else:
    import _pairDist
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


import vec3
class PairDist(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _pairDist.new_PairDist(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def probWidth(self, *args, **kwargs):
        return _pairDist.PairDist_probWidth(self, *args, **kwargs)

    def setProbWidth(self, *args, **kwargs):
        return _pairDist.PairDist_setProbWidth(self, *args, **kwargs)

    def probType(self, *args, **kwargs):
        return _pairDist.PairDist_probType(self, *args, **kwargs)

    def setProbType(self, *args, **kwargs):
        return _pairDist.PairDist_setProbType(self, *args, **kwargs)

    def calcGradient(self, *args, **kwargs):
        return _pairDist.PairDist_calcGradient(self, *args, **kwargs)

    def setCalcGradient(self, *args, **kwargs):
        return _pairDist.PairDist_setCalcGradient(self, *args, **kwargs)

    def updateValues(self, *args, **kwargs):
        return _pairDist.PairDist_updateValues(self, *args, **kwargs)

    def Pr(self, *args, **kwargs):
        return _pairDist.PairDist_Pr(self, *args, **kwargs)

    def getGradient(self, *args, **kwargs):
        return _pairDist.PairDist_getGradient(self, *args, **kwargs)
    __swig_destroy__ = _pairDist.delete_PairDist
    __del__ = lambda self: None

class PairDistPtr(PairDist):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = PairDist

PairDist_swigregister = _pairDist.PairDist_swigregister
PairDist_swigregister(PairDist)

class CDSList_IntVec3Pair(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __len__(self, *args, **kwargs):
        return _pairDist.CDSList_IntVec3Pair___len__(self, *args, **kwargs)

    def __init__(self, *args):
        this = _pairDist.new_CDSList_IntVec3Pair(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def __getitem__(self, *args):
        return _pairDist.CDSList_IntVec3Pair___getitem__(self, *args)

    def __delitem__(self, *args, **kwargs):
        return _pairDist.CDSList_IntVec3Pair___delitem__(self, *args, **kwargs)

    def append(self, *args, **kwargs):
        return _pairDist.CDSList_IntVec3Pair_append(self, *args, **kwargs)

    def remove(self, *args, **kwargs):
        return _pairDist.CDSList_IntVec3Pair_remove(self, *args, **kwargs)

    def removeAll(self, *args, **kwargs):
        return _pairDist.CDSList_IntVec3Pair_removeAll(self, *args, **kwargs)

    def __setitem__(self, *args, **kwargs):
        return _pairDist.CDSList_IntVec3Pair___setitem__(self, *args, **kwargs)

    def __getslice__(self, *args, **kwargs):
        return _pairDist.CDSList_IntVec3Pair___getslice__(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _pairDist.CDSList_IntVec3Pair_help(self, *args, **kwargs)
    __swig_destroy__ = _pairDist.delete_CDSList_IntVec3Pair
    __del__ = lambda self: None

class CDSList_IntVec3PairPtr(CDSList_IntVec3Pair):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = CDSList_IntVec3Pair

CDSList_IntVec3Pair_swigregister = _pairDist.CDSList_IntVec3Pair_swigregister
CDSList_IntVec3Pair_swigregister(CDSList_IntVec3Pair)


pyXplorHelp = help



