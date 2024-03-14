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
        mname = '.'.join((pkg, '_xplorWrap')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_xplorWrap')
    _xplorWrap = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_xplorWrap', [dirname(__file__)])
        except ImportError:
            import _xplorWrap
            return _xplorWrap
        try:
            _mod = imp.load_module('_xplorWrap', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _xplorWrap = swig_import_helper()
    del swig_import_helper
else:
    import _xplorWrap
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


class XplorOutputState(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _xplorWrap.new_XplorOutputState(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    echo = _swig_property(_xplorWrap.XplorOutputState_echo_get, _xplorWrap.XplorOutputState_echo_set)
    wrnlev = _swig_property(_xplorWrap.XplorOutputState_wrnlev_get, _xplorWrap.XplorOutputState_wrnlev_set)
    __swig_destroy__ = _xplorWrap.delete_XplorOutputState
    __del__ = lambda self: None

class XplorOutputStatePtr(XplorOutputState):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = XplorOutputState

XplorOutputState_swigregister = _xplorWrap.XplorOutputState_swigregister
XplorOutputState_swigregister(XplorOutputState)

class XplorWrap(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _xplorWrap.new_XplorWrap(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
    __swig_destroy__ = _xplorWrap.delete_XplorWrap
    __del__ = lambda self: None

    def simulation(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_simulation(self, *args, **kwargs)

    def command(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_command(self, *args, **kwargs)

    def fastCommand(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_fastCommand(self, *args, **kwargs)

    def shell(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_shell(self, *args, **kwargs)

    def deleteAtoms(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_deleteAtoms(self, *args, **kwargs)

    def openFile(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_openFile(self, *args, **kwargs)

    def closeFile(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_closeFile(self, *args, **kwargs)

    def energy(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_energy(self, *args, **kwargs)

    def derivs(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_derivs(self, *args, **kwargs)

    def select(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_select(self, *args, **kwargs)

    def disableOutput(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_disableOutput(self, *args, **kwargs)

    def enableOutput(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_enableOutput(self, *args, **kwargs)

    def setRandomSeed(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_setRandomSeed(self, *args, **kwargs)

    def randomSeed(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_randomSeed(self, *args, **kwargs)

    def uniformRandom(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_uniformRandom(self, *args, **kwargs)
    xplorVars = staticmethod(_xplorWrap.XplorWrap_xplorVars)
    resetXplorVars = staticmethod(_xplorWrap.XplorWrap_resetXplorVars)
    local = _swig_property(_xplorWrap.XplorWrap_local_get, _xplorWrap.XplorWrap_local_set)

    def help(self, *args, **kwargs):
        return _xplorWrap.XplorWrap_help(self, *args, **kwargs)

class XplorWrapPtr(XplorWrap):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = XplorWrap

XplorWrap_swigregister = _xplorWrap.XplorWrap_swigregister
XplorWrap_swigregister(XplorWrap)

def XplorWrap_xplorVars(*args):
    return _xplorWrap.XplorWrap_xplorVars(*args)
XplorWrap_xplorVars = _xplorWrap.XplorWrap_xplorVars

def XplorWrap_resetXplorVars(*args, **kwargs):
    return _xplorWrap.XplorWrap_resetXplorVars(*args, **kwargs)
XplorWrap_resetXplorVars = _xplorWrap.XplorWrap_resetXplorVars


pyXplorHelp = help


def help(*args):
    return _xplorWrap.help(*args)
help = _xplorWrap.help

