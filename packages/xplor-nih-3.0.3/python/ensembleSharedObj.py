"""create a Python object shared across an EnsembleSimulation

Method of sharing python objects across between processes associated with
<m ensembleSimulation>.EnsembleSimulations.
"""

import pickle
from os import unlink, getpid
import collections

class SharedObj:
    """
    a Python object shared across ensemble simulation members.
    The shared object must be picklable- meaning it should be a Python type.
    Strange memory errors will occur if non-picklable objects are shared.
    SWIG objects do not seem to be pickable.
    """
    def __init__(s,ens=0,obj=0):
        """
        arguments:
          ens - optional EnsembleSimulation object.
          obj - optional initial value for variable.
        """
        from tempfile import NamedTemporaryFile

        if not ens:
            from ensembleSimulation import EnsembleSimulation_currentSimulation
            ens = EnsembleSimulation_currentSimulation()
            
            if not ens:
                from ensembleSimulation import EnsembleSimulation_sizeOneSimulation
                ens=EnsembleSimulation_sizeOneSimulation()
            pass

        s.ens = ens
        

        sharedString = ens.sharedString( 1024 )

        if ens.singleThread():
            file = NamedTemporaryFile(prefix="%ld_"%getpid() )
            s.filename = file.name
            sharedString.set(s.filename)
            file.write( pickle.dumps(obj) )
            # so that destructor always works ok.
            del file
            open(s.filename,'w').write('something\n') 
            pass
        s.ens.multiThread()
        
        s.filename = sharedString.get()

        #FIX: this shouldn't be necessary
        sharedString.thisown=1
        return
    def __del__(s):
        #can't use s.ens, because it may already have been deleted.
        from ensembleSimulation import singleThread, multiThread
        if singleThread():
            unlink(s.filename)
        multiThread()
        return
    def set(s,obj):
        " set the value of the shared object"
        open(s.filename,'wb').write( pickle.dumps(obj) )
        return
    def get(s):
        """ return the value of the shared object. """
        return pickle.loads( open(s.filename,"rb").read() )
    def barrierGet(s):
        """ return the non-shared value of the shared object w/ barrier.
        Must be called when all threads are active.
        """
        s.ens.barrier()
        ret=s.get()
        s.ens.barrier()
        return ret
    pass

#FIX: need to add regression test
from ensembleSimulation import singleThread, multiThread
class Random:
    """an ensemble-safe Random object- the same result from random is
    distributed throughout the <m ensembleSimulation>.EnsembleSimulation
    """
    def __init__(s):
        s.shared = SharedObj()
        import random
        methods=[x for x in dir(random) if not x.startswith('_')]
        for m in methods:
            if isinstance(eval("random.%s"%m), collections.Callable):
                exec("""
def %s(*args):
    from ensembleSharedObj import SharedObj
    shared = SharedObj()
    from ensembleSimulation import singleThread, multiThread
    if singleThread():
      import random
      ret =  random.%s(*args)
      shared.set(ret)
    multiThread()
    ret = shared.barrierGet()
    return ret"""% (m,m),s.__dict__)
        return
    pass

        


# member 0 generates a random filename
# write 0 to it and close
# distribute the filename to the other members

# set: write pickled object to file, close it

# get: read pickled object from file

def collect(esim,
            val):
    """
    Return a list of size esim.size(), containing each value val from
    each process of an <m ensembleSimulation>.EnsembleSimulation.

    val must be picklable.
    """

    sharedVals=[]
    for i in range(esim.size()):
        sharedVals.append( SharedObj(esim) )
        pass
    sharedVals[esim.member().memberIndex()].set(val)
    
    ret=[]
    for sharedVal in sharedVals:
        ret.append( sharedVal.barrierGet() )
        pass

    return ret
