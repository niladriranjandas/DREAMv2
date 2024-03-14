"""Xplor-NIH interface to the VMD-XPLOR visualization package

constructor:

  VMDInter(host,port) -- see __init__ below.


methods:       

makeObj   - defined below
getObj    - defined below
loadFiles - defined below

tclCommand(cmd) -
  execute a tcl command within VMD
  
deleteAll() -
  delete all VMD-XPLOR objects.

deleteObj(name) -
  delete named VMD-XPLOR object.

makeTop(name) -
  make the named VMD-XPLOR object `top,' or active.

setCenter(vec3) -
  set the center point for VMD rotations.

setColor(name,color) -
  set the color of the named VMD-XPLOR object.

vmdLocal() -
   boolean - returns true if VMD-XPLOR is running on the local machine.


the following methods should usually not be called directly, rather through
the VMDObject interface:

  appendFrame()
  labels()
  loadFrame()

"""
from publicVMDInter import PublicVMDInter
from simulation import currentSimulation
from atomSel import AtomSel
import atom
import simulation
import vec3
import os
from pdbTool import PDBTool

class VMDInter(PublicVMDInter):
    def __init__(s,host="",port=-1):
        """the arguments are optional.
        host: where VMD-XPLOR is running. defaults to localhost or
        the host specified using -host on the command line.
        port: the port to which VMD-XPLOR is listening. Defaults to the
        value given using -port on the command line.
        """
        if host=="":
            try:
                host = os.environ["DP_HOST"]
            except:
                host = "localhost"
                pass
            pass
        if port==-1:
            try:
                port = int(os.environ["DP_PORT"])
            except:
                raise "VMDInter: no port specified"
                pass
            pass
        PublicVMDInter.__init__(s,host,port)
        s.labelSegID    = 0
        s.labelResName  = 0
        s.labelResID    = 1
        s.labelAtomName = 1
        s.labelSize     = 1
        return
    def makeObj(s,name):
        """create a graphical object in VMD with the given name"""
        ok=0
        cnt=0
        name =name.replace(' ','_')   #spaces invalid
        while not ok:
            gname=name
            try:
                if cnt>0: gname += "_%d" % cnt
                s.addObject(gname)
                ret= VMDObject(s,gname)
                ok=1
            except SystemError as excpt:
                print(excpt)
                pass
            cnt += 1
            pass
        return ret
    def getObj(s,name):
        """make the named object the top object and return it's Xplor-NIH
        representation.
        """
        s.makeTop(name) # will throw exception if name doesn't exist
        return VMDObject(s,name)
    def loadFiles(s,
                  filename,
                  structRange=[0],
                  command="",
                  color="",
                  displaySel="known",
                  loadSel="known"):
        """ loadFiles(filename[,optional keyword arguments])
               load the given files into an VMD representation. Filename can
               be a single file, or a template containing the string
               STRUCTURE, which is replaced by an integer from the structure
               range.
               keywords:
                  structRange: a list of integers used to generate filename.
                  command:     a string or function exec'ed after coordinates
                               are read, and before structures are loaded.
                  color:       color to set loaded structures.
                  displaySel:  string selecting which atoms are show in
                               the VMD representation.
                  loadSel:     selection of which atoms are read from the
                               coordinate file.
        """
        template = filename
        objs = []
        import os.path
        for cnt in structRange:
            filename = template.replace("STRUCTURE","%d" %cnt)
            PDBTool(filename,AtomSel(loadSel)).read()
            if command and type(command)==type("string"):
                global_dict = getouterframes( currentframe() )[1][0].f_globals
                local_dict = getouterframes( currentframe() )[1][0].f_locals
                exec( command, global_dict, local_dict )
            elif command:
                command()
                pass
            obj = s.makeObj(os.path.split(filename)[-1])
            obj.bonds(displaySel)
            if color: obj.color(color)
            objs.append(obj)
            pass
        return objs
    def help(s):
        return open(os.environ["XPLOR_HELPDIR"]+"nih-py-vmdInter").read()
    pass
        



#
class VMDObject:
    def __init__(s,conn,name):
        """these object are created using the makeObj method of VMDInter.
        """
        s.conn = conn
        s.name = name
        s.labelSegID    = conn.labelSegID    
        s.labelResName  = conn.labelResName  
        s.labelResID    = conn.labelResID    
        s.labelAtomName = conn.labelAtomName 
        s.labelSize     = conn.labelSize     
        return
    def color(s,colorSpec):
        """set color (allow rgb or color name)
        """
        s.conn.setColor(s.name,colorSpec)
        return
    def delete(s):
        """delete the associated vmd object
        """
        s.conn.deleteObj(s.name)
        return
    def labels(s,sel):
        """
         make vmd text objects
           the following boolean members are consulted:
             labelSegID      
             labelResName    
             labelResID
             labelAtomName
        """
        if type(sel) == type("string"):
            sel = AtomSel(sel)
            pass

        # might waht to break style out to
        # separate call
        #construct label list
        sim = sel.simulation()
        labelList = []
        for i in sel.indices():
            atom = sim.atomByID(i)
            if not atom.isValid():
                continue
            text = ""
            if s.labelSegID:    text += atom.segmentName() + ":"
            if s.labelResName:  text += atom.residueName()
            if s.labelResID:    text += repr(atom.residueNum()) + ":"
            if s.labelAtomName: text += atom.atomName()
            if text[-1]==":": text = text[:-1]
            labelList.append( (atom.pos(),text,s.labelSize) )
            pass
        s.conn.labels(s.name,labelList)
        return
    def bonds(s,aSel,style="lines",colorBy="Molecule"):
        """create a bonded representation between atoms specified in
        the <m atomSel>.AtomSel aSel.

        style   is the name of the VMD drawing method.
        colorBy is the VMD coloring method.
        
        """
        s.aSel = aSel
        if type(aSel).__name__ == type("string").__name__:
            s.aSel = AtomSel(aSel)
            pass
        # loadFrame checks that atoms are not INVALID 
        s.conn.loadFrame(s.name,s.aSel,style,colorBy)
        return
    def append(s):
        """append bond frames to this object."""
        s.conn.appendFrame(s.name,s.aSel)
        return
    def lines(s,sel):
        """not yet implemented.
        """
        # do bonds (hbonds ps>) command
        print("not implemented")
        return

    def help(s):
        return open(os.environ["XPLOR_HELPDIR"]+"nih-vmdObject").read()
    pass

class VMDTraj:
    """a trajectory object to be used with an <m ivm>.IVM object, used to
    display coordinates in a VMD window during a trajectory run.
    Use the IVM setTrajectory method to to view a given ivm trajectory.
    """
    def __init__(s,vmdObj,saveInterval=1):
        """specify the vmdObject to append frames. The frequency of update is
        controlled via the saveInterval argument (defaults to 1 [every frame]).
        """
        s.vmdObj = vmdObj
        s.saveInterval = saveInterval
        return
    def write(s,iter,stepsize):
        if iter%s.saveInterval != 0: return

        if iter==0:
            s.vmdObj.bonds(s.vmdObj.aSel)
        else:
            s.vmdObj.append()
            pass
        return
    pass

        
            


def help():
    return open(os.environ["XPLOR_HELPDIR"]+"nih-py-vmdInter").read()
pyXplorHelp=help
