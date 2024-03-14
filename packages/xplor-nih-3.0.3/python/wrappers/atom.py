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
        mname = '.'.join((pkg, '_atom')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_atom')
    _atom = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_atom', [dirname(__file__)])
        except ImportError:
            import _atom
            return _atom
        try:
            _mod = imp.load_module('_atom', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _atom = swig_import_helper()
    del swig_import_helper
else:
    import _atom
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



from vec3 import Vec3

class Atom(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args):
        this = _atom.new_Atom(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def segmentName(self, *args, **kwargs):
        return _atom.Atom_segmentName(self, *args, **kwargs)

    def residueName(self, *args, **kwargs):
        return _atom.Atom_residueName(self, *args, **kwargs)

    def atomName(self, *args, **kwargs):
        return _atom.Atom_atomName(self, *args, **kwargs)

    def chemType(self, *args, **kwargs):
        return _atom.Atom_chemType(self, *args, **kwargs)

    def residueNum(self, *args, **kwargs):
        return _atom.Atom_residueNum(self, *args, **kwargs)

    def pos_ref(self, *args, **kwargs):
        return _atom.Atom_pos_ref(self, *args, **kwargs)

    def vel_ref(self, *args, **kwargs):
        return _atom.Atom_vel_ref(self, *args, **kwargs)

    def mass(self, *args, **kwargs):
        return _atom.Atom_mass(self, *args, **kwargs)

    def fric(self, *args, **kwargs):
        return _atom.Atom_fric(self, *args, **kwargs)

    def charge(self, *args, **kwargs):
        return _atom.Atom_charge(self, *args, **kwargs)

    def setSegmentName(self, *args, **kwargs):
        return _atom.Atom_setSegmentName(self, *args, **kwargs)

    def setResidueName(self, *args, **kwargs):
        return _atom.Atom_setResidueName(self, *args, **kwargs)

    def setAtomName(self, *args, **kwargs):
        return _atom.Atom_setAtomName(self, *args, **kwargs)

    def setChemType(self, *args, **kwargs):
        return _atom.Atom_setChemType(self, *args, **kwargs)

    def setResidueNum(self, *args, **kwargs):
        return _atom.Atom_setResidueNum(self, *args, **kwargs)

    def setPos(self, *args, **kwargs):
        return _atom.Atom_setPos(self, *args, **kwargs)

    def setVel(self, *args, **kwargs):
        return _atom.Atom_setVel(self, *args, **kwargs)

    def setMass(self, *args, **kwargs):
        return _atom.Atom_setMass(self, *args, **kwargs)

    def setFric(self, *args, **kwargs):
        return _atom.Atom_setFric(self, *args, **kwargs)

    def setCharge(self, *args, **kwargs):
        return _atom.Atom_setCharge(self, *args, **kwargs)

    def string(self, *args, **kwargs):
        return _atom.Atom_string(self, *args, **kwargs)

    def simulation(self, *args):
        return _atom.Atom_simulation(self, *args)

    def index(self, *args, **kwargs):
        return _atom.Atom_index(self, *args, **kwargs)

    def isValid(self, *args, **kwargs):
        return _atom.Atom_isValid(self, *args, **kwargs)

    def bondedTo(self, *args, **kwargs):
        return _atom.Atom_bondedTo(self, *args, **kwargs)

    def __eq__(self,other):
      if type(self) != type(other):
          return False

      if (self.simulation().atomID()==other.simulation().atomID() and
          self.index()==other.index()                                 ):
        return True
      else:
        return False
    def __ne__(self,other):
      return not self.__eq__(other)
    def pos(self):
      return Vec3(self.pos_ref())
    def vel(self):
      return Vec3(self.vel_ref())


    def help(self, *args, **kwargs):
        return _atom.Atom_help(self, *args, **kwargs)
    __swig_destroy__ = _atom.delete_Atom
    __del__ = lambda self: None

class AtomPtr(Atom):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = Atom

Atom_swigregister = _atom.Atom_swigregister
Atom_swigregister(Atom)
cvar = _atom.cvar
Atom.INVALID_COORD = _atom.cvar.Atom_INVALID_COORD


def cdsMapConvertToInt(*args, **kwargs):
    return _atom.cdsMapConvertToInt(*args, **kwargs)
cdsMapConvertToInt = _atom.cdsMapConvertToInt

pyXplorHelp = help


def help(*args):
    return _atom.help(*args)
help = _atom.help

