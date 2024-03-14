"""tools to perform logical operations on atom selections

NOTE: not tested! These routines need to be fixed.
"""

#
# FIX: should use AtomSel objects
# select() should be renamed.
# should this be placed in selectTools.py?
#


from simulation import currentSimulation

class BoolArray:
    """Array whose values are zero or one. The array's size is determined by
    the size argument to the constructor. All members of the array are zero
    except those corresponding to indices specified by the second argument to
    the constructor.

    Simple boolean operations can be performed on BoolArrays.
    An array of the nonzero elements can be obtained using the simpleArray
    method."""
    def __init__(s,size,arr):
        s.arr = []
        for i in range(0,size):
            s.arr.append(0)
            pass
        for i in arr:
            s.arr[i] = 1
            pass
        return
    def simpleArray(s):
        ret = []
        for i in range(0,len(s.arr)):
            if s.arr[i]==1:
                ret.append(i)
                pass
            pass
        return ret
    def __and__(s,x):
        ret = BoolArray(len(s.arr),[])
        for i in range(0,len(s.arr)):
            if (s.arr[i] and x.arr[i]):
                ret.arr[i]=1
                pass
            pass
        return ret
    def __or__(s,x):
        ret = BoolArray(len(s.arr),[])
        for i in range(0,len(s.arr)):
            if (s.arr[i] or x.arr[i]):
                ret.arr[i]=1
                pass
            pass
        return ret
    def __neg__(s):
        ret = BoolArray(len(s.arr),[])
        for i in range(0,len(s.arr)):
            if not s.arr[i]:
                ret.arr[i] = 1
            else:
                ret.arr[i] = 0
                pass
            pass
        return ret
            
    pass
    
            

def select(str,strNS=None,size=None):
    """perform logical operations on select arrays.
    For instance, if
     sel1 = select('name CA')
     sel2 = select('name N')
    then
     selectLogic('sel1 or sel2')
    returns the indices of all CA and N atoms.
    """
    import re
    if strNS is None: strNS = globals()
    if size is None: size = currentSimulation().numAtoms()
    ns = {}
    ns['__name__'] = "select namespace"
    names = re.findall(r"\b[a-zA-Z0-9_]+\b",str)
    for name in names:
        if not name in ("and","or","not"):
            ns[name] = BoolArray(size,strNS[name])
            pass
        pass

    # replace \band\b w/ &
    str = re.sub(r"\band\b","&",str)
    # replace \bor\b w/ |
    str = re.sub(r"\bor\b","|",str)
    # replace \bnot\ w/ -
    str = re.sub(r"\bnot\b","-",str)
    ret = eval(str,ns)
    
    return ret.simpleArray()

#selectLogic(str,

