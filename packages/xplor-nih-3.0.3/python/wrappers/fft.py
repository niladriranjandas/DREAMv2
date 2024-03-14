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
        mname = '.'.join((pkg, '_fft')).lstrip('.')
        try:
            return importlib.import_module(mname)
        except ImportError:
            return importlib.import_module('_fft')
    _fft = swig_import_helper()
    del swig_import_helper
elif _swig_python_version_info >= (2, 6, 0):
    def swig_import_helper():
        from os.path import dirname
        import imp
        fp = None
        try:
            fp, pathname, description = imp.find_module('_fft', [dirname(__file__)])
        except ImportError:
            import _fft
            return _fft
        try:
            _mod = imp.load_module('_fft', fp, pathname, description)
        finally:
            if fp is not None:
                fp.close()
        return _mod
    _fft = swig_import_helper()
    del swig_import_helper
else:
    import _fft
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


__sthead_hh__ = _fft.__sthead_hh__

def omp_get_thread_num():
    return _fft.omp_get_thread_num()
omp_get_thread_num = _fft.omp_get_thread_num

def omp_get_max_threads():
    return _fft.omp_get_max_threads()
omp_get_max_threads = _fft.omp_get_max_threads
FALSE = _fft.FALSE
TRUE = _fft.TRUE
PI = _fft.PI
class FFT_FLOATTYPE(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _fft.new_FFT_FLOATTYPE(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def fft(self, *args, **kwargs):
        return _fft.FFT_FLOATTYPE_fft(self, *args, **kwargs)

    def fft_inv(self, *args, **kwargs):
        return _fft.FFT_FLOATTYPE_fft_inv(self, *args, **kwargs)

    def resize(self, *args, **kwargs):
        return _fft.FFT_FLOATTYPE_resize(self, *args, **kwargs)
    __swig_destroy__ = _fft.delete_FFT_FLOATTYPE
    __del__ = lambda self: None

class FFT_FLOATTYPEPtr(FFT_FLOATTYPE):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = FFT_FLOATTYPE

FFT_FLOATTYPE_swigregister = _fft.FFT_FLOATTYPE_swigregister
FFT_FLOATTYPE_swigregister(FFT_FLOATTYPE)

class RFFT_FLOATTYPE(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __init__(self, *args, **kwargs):
        this = _fft.new_RFFT_FLOATTYPE(*args, **kwargs)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def fft(self, *args, **kwargs):
        return _fft.RFFT_FLOATTYPE_fft(self, *args, **kwargs)

    def fft_inv(self, *args, **kwargs):
        return _fft.RFFT_FLOATTYPE_fft_inv(self, *args, **kwargs)

    def resize(self, *args, **kwargs):
        return _fft.RFFT_FLOATTYPE_resize(self, *args, **kwargs)
    __swig_destroy__ = _fft.delete_RFFT_FLOATTYPE
    __del__ = lambda self: None

class RFFT_FLOATTYPEPtr(RFFT_FLOATTYPE):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = RFFT_FLOATTYPE

RFFT_FLOATTYPE_swigregister = _fft.RFFT_FLOATTYPE_swigregister
RFFT_FLOATTYPE_swigregister(RFFT_FLOATTYPE)

class CDSVector_DComplex(object):
    thisown = _swig_property(lambda x: x.this.own(), lambda x, v: x.this.own(v), doc='The membership flag')
    __repr__ = _swig_repr

    def __len__(self, *args, **kwargs):
        return _fft.CDSVector_DComplex___len__(self, *args, **kwargs)

    def resize(self, *args, **kwargs):
        return _fft.CDSVector_DComplex_resize(self, *args, **kwargs)

    def set(self, *args, **kwargs):
        return _fft.CDSVector_DComplex_set(self, *args, **kwargs)

    def __init__(self, *args):
        this = _fft.new_CDSVector_DComplex(*args)
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this

    def get(self, *args, **kwargs):
        return _fft.CDSVector_DComplex_get(self, *args, **kwargs)

    def fromList(s,l):
        s.resize(len(l)) 
        for i in range( len(s) ):
            s[i] = l[i]
        return s


    def __setitem__(self, *args, **kwargs):
        return _fft.CDSVector_DComplex___setitem__(self, *args, **kwargs)

    def __getslice__(self, *args, **kwargs):
        return _fft.CDSVector_DComplex___getslice__(self, *args, **kwargs)

    def __getitem__(self, *args, **kwargs):
        arg = args[0]
        if type(arg) is slice:
            if not (arg.step==None or arg.step==1):
                raise Exception("slice step!=1 not supported: " +
                                str(arg.step))
            start=arg.start if arg.start!=None else 0
            stop=arg.stop if arg.stop!=None else len(self)
            return self.__getslice__(start,stop)
        else:
            return self.get(arg)


    def __pow__(self, *args, **kwargs):
        return _fft.CDSVector_DComplex___pow__(self, *args, **kwargs)

    def help(self, *args, **kwargs):
        return _fft.CDSVector_DComplex_help(self, *args, **kwargs)

    def scale(self, *args, **kwargs):
        return _fft.CDSVector_DComplex_scale(self, *args, **kwargs)

    def __add__(self, *args):
        return _fft.CDSVector_DComplex___add__(self, *args)

    def __radd__(self, *args, **kwargs):
        return _fft.CDSVector_DComplex___radd__(self, *args, **kwargs)

    def __sub__(self, *args):
        return _fft.CDSVector_DComplex___sub__(self, *args)

    def __rsub__(self, *args, **kwargs):
        return _fft.CDSVector_DComplex___rsub__(self, *args, **kwargs)

    def __mul__(self, *args):
        return _fft.CDSVector_DComplex___mul__(self, *args)

    def __rmul__(self, *args, **kwargs):
        return _fft.CDSVector_DComplex___rmul__(self, *args, **kwargs)

    def __truediv__(self, *args):
        return _fft.CDSVector_DComplex___truediv__(self, *args)

    def __rtruediv__(self, *args, **kwargs):
        return _fft.CDSVector_DComplex___rtruediv__(self, *args, **kwargs)

    def __imul__(self, *args, **kwargs):
        return _fft.CDSVector_DComplex___imul__(self, *args, **kwargs)

    def __itruediv__(self, *args, **kwargs):
        return _fft.CDSVector_DComplex___itruediv__(self, *args, **kwargs)

    def __iadd__(self, *args):
        return _fft.CDSVector_DComplex___iadd__(self, *args)

    def __isub__(self, *args):
        return _fft.CDSVector_DComplex___isub__(self, *args)

    def __neg__(self, *args, **kwargs):
        return _fft.CDSVector_DComplex___neg__(self, *args, **kwargs)

    def __str__(self, *args, **kwargs):
        return _fft.CDSVector_DComplex___str__(self, *args, **kwargs)
    __swig_destroy__ = _fft.delete_CDSVector_DComplex
    __del__ = lambda self: None

class CDSVector_DComplexPtr(CDSVector_DComplex):
    def __init__(self, this):
        try:
            self.this.append(this)
        except __builtin__.Exception:
            self.this = this
        self.this.own(0)
        self.__class__ = CDSVector_DComplex

CDSVector_DComplex_swigregister = _fft.CDSVector_DComplex_swigregister
CDSVector_DComplex_swigregister(CDSVector_DComplex)


def min(*args, **kwargs):
    return _fft.min(*args, **kwargs)
min = _fft.min

def max(*args, **kwargs):
    return _fft.max(*args, **kwargs)
max = _fft.max

def sum(*args, **kwargs):
    return _fft.sum(*args, **kwargs)
sum = _fft.sum

def cat(*args, **kwargs):
    return _fft.cat(*args, **kwargs)
cat = _fft.cat

def cat3(*args, **kwargs):
    return _fft.cat3(*args, **kwargs)
cat3 = _fft.cat3

def cat4(*args, **kwargs):
    return _fft.cat4(*args, **kwargs)
cat4 = _fft.cat4

import cdsVector

FFT  =  FFT_FLOATTYPE
RFFT = RFFT_FLOATTYPE

default_fft =   FFT()
default_rfft = RFFT()

fft     = default_fft.fft
fft_inv = default_fft.fft_inv

rfft     = default_rfft.fft
rfft_inv = default_rfft.fft_inv



pyXplorHelp = help


def help(*args):
    return _fft.help(*args)
help = _fft.help


