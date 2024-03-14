"""
 Tools for manipulating torsion angles
"""
# Definition of the DICT in the script
"""
torStatCnt  :dict of the form X[UID][State]=count where X[UID] is a dictionary
             and X[UID][State] another dictionary
RotAngleName:dict of the from X[UID]=angleName which is a List of tuple
             of the form [('CHI1',('N', 'CA', 'CB', 'OG')),
             ('CHI2',('CA','CB', 'OG', 'HG'))]
torStateIndex:dict of the form X[UID][State]=Index which is a list of integers
              that stores the structure index value in the atomCoordinates
              List.
scTorsionAtoms: a dict of the form
                X['SER'] = [('CHI1',('N', 'CA', 'CB', 'OG')),
                            ('CHI2',('CA','CB', 'OG', 'HG'))]              
"""
from sys import argv
from selectTools import minResid, maxResid
from dihedral import Dihedral
from atomSel import AtomSel

def getAngles(resType,resid,segid):
    segSel = 'and segid "%s"' % segid
    angles = [("phi",  ("name C  and resid %d %s" % (resid-1,segSel),
                        "name N  and resid %d %s" % (resid,segSel),
                        "name CA and resid %d %s" % (resid,segSel),
                        "name C  and resid %d %s" % (resid,segSel) )),
              ("psi",  ("name N  and resid %d %s" % (resid,segSel),
                        "name CA and resid %d %s" % (resid,segSel),
                        "name C  and resid %d %s" % (resid,segSel),
                        "name N  and resid %d %s" % (resid+1,segSel))),
              ("omega",("name CA and resid %d %s" % (resid,segSel),
                        "name C  and resid %d %s" % (resid,segSel),
                        "name N  and resid %d %s" % (resid+1,segSel),
                        "name CA and resid %d %s" % (resid+1,segSel)))]
    try:
        for (name,atoms) in scTorsionAtoms[resType]:
            angles.append((name,
                           ("name %s and resid %d %s" % (atoms[0],
                                                         resid,segSel),
                            "name %s and resid %d %s" % (atoms[1],
                                                         resid,segSel),
                            "name %s and resid %d %s" % (atoms[2],
                                                         resid,segSel),
                            "name %s and resid %d %s" % (atoms[3],
                                                             resid,segSel))))
            pass
        pass
    except KeyError:
        print("unsupported residue type:", resType)
        pass
        
    return angles

#Edited by Robin to remove torsion angle defined by protons
scTorsionAtoms = {}
scTorsionAtoms['GLY'] = []
scTorsionAtoms['ALA'] = []
scTorsionAtoms['SER'] = [('CHI1',('N', 'CA', 'CB', 'OG'))]
                        # [('CHI2',('CA','CB', 'OG', 'HG'))]
scTorsionAtoms['THR'] = [('CHI1' , ('N' , 'CA', 'CB', 'OG1'))]
                        #[('CHI21', ('CA', 'CB', 'OG1', 'HG1'))]
scTorsionAtoms['LYS'] = [('CHI1', ('N',  'CA', 'CB', 'CG' )),
                         ('CHI2', ('CA', 'CB', 'CG', 'CD' )),
                         ('CHI3', ('CB', 'CG', 'CD', 'CE' )),
                         ('CHI4', ('CG', 'CD', 'CE', 'NZ' ))]
scTorsionAtoms['CYS'] =[('CHI1', ('N',  'CA', 'CB', 'SG'))]
                       #[('CHI2', ('CA', 'CB', 'SG', 'HG'))]
scTorsionAtoms['MET'] = [('CHI1', ('N',  'CA', 'CB', 'CG' )),
                         ('CHI2', ('CA', 'CB', 'CG', 'SD' )),
                         ('CHI3', ('CB', 'CG', 'SD', 'CE' ))]
scTorsionAtoms['VAL'] =[('CHI1',  ('N',  'CA', 'CB',  'CG1' ))]
scTorsionAtoms['ILE'] = [('CHI1',  ('N',  'CA',  'CB',  'CG1')),
                         ('CHI21', ('CA', 'CB',  'CG1', 'CD1'))]
scTorsionAtoms['LEU'] = [('CHI1',  ('N',  'CA',  'CB',  'CG'   )),
                         ('CHI2',  ('CA', 'CB',  'CG',  'CD1'  ))]
scTorsionAtoms['ASP'] = [('CHI1',  ('N',  'CA', 'CB',  'CG'  )),
                         ('CHI2',  ('CA', 'CB', 'CG',  'OD1' ))]
scTorsionAtoms['ASN'] = [('CHI1',  ('N',  'CA', 'CB',  'CG'  )),
                         ('CHI2',  ('CA', 'CB', 'CG',  'OD1' ))]
scTorsionAtoms['GLU'] = [('CHI1',  ('N',  'CA', 'CB',  'CG'  )),
                         ('CHI2',  ('CA', 'CB', 'CG',  'CD'  )),
                         ('CHI3',  ('CB', 'CG', 'CD',  'OE1' ))]
scTorsionAtoms['GLN'] = [('CHI1',  ('N',  'CA', 'CB',  'CG'  )),
                         ('CHI2',  ('CA', 'CB', 'CG',  'CD'  )),
                         ('CHI3',  ('CB', 'CG', 'CD',  'OE1' ))]
scTorsionAtoms['ARG'] = [('CHI1',  ('N',  'CA', 'CB',  'CG'  )),
                         ('CHI2',  ('CA', 'CB', 'CG',  'CD'  )),
                         ('CHI3',  ('CB', 'CG', 'CD',  'NE'  )),
                         ('CHI4',  ('CG', 'CD', 'NE',  'CZ'  ))]
scTorsionAtoms['PRO'] = []
scTorsionAtoms['HIS'] = [('CHI1',  ('N',  'CA', 'CB', 'CG'  )),
                         ('CHI2',  ('CA', 'CB', 'CG', 'ND1' ))]
scTorsionAtoms['PHE'] = [('CHI1',  ('N',  'CA', 'CB', 'CG' )),
                         ('CHI2',  ('CA', 'CB', 'CG', 'CD1'))]
scTorsionAtoms['TYR'] = [('CHI1',   ('N',   'CA', 'CB', 'CG' )),
                         ('CHI2',   ('CA',  'CB', 'CG', 'CD1'))]
                         #[('CHI6',   ('CE1', 'CZ', 'OH', 'HH' ))]
scTorsionAtoms['TRP'] = [('CHI1',  ('N',  'CA', 'CB', 'CG' )),
                         ('CHI2',  ('CA', 'CB', 'CG', 'CD1'))]

def getscTorsionAtoms(resname):
    return scTorsionAtoms[resname]


binData={}

binData['PHE_CHI2'] = (('+90',90), ('-90',-90))
binData['TYR_CHI2'] = (('+90',90), ('-90',-90))
binData['TRP_CHI2'] = (('+90',90), ('-90',-90))
binData['HIS_CHI2'] = (('+90',90), ('-90',-90))
binData['ASP_CHI2'] = (('0',0), ('180',180))
binData['ASN_CHI2'] = (('0',0), ('180',180))
binData['GLU_CHI3'] = (('0',0), ('180',180))
binData['GLN_CHI3'] = (('0',0), ('180',180))

def getBranch(theta1,theta2):
    """return the branch number n such that n*2*PI+theta1-theta2 is
    minimal"""
    return round((theta2-theta1)/360)

def getBin(resType,taname,val):
    closest = ('single',0,val)
    bins=[]
    if taname.startswith('CHI'):
        bins = (('180',180), ('+60',60), ('-60',-60))
        pass
    try:
        bins = binData[resType+'_'+taname]
    except KeyError:
        pass
    
    if bins:
        closest = ('',361,361)
        for (name,mean) in bins:
            val += 360*getBranch(val,mean)
            dist=abs(val-mean)
            if dist<closest[1]:
                closest = (name,dist,val)
                pass
            #print 'branch',mean,val,getBranch(val,mean), closest
            #print closest
            pass
        pass
        
    return closest
        
def getangleBins():
    return angleBins
    
angleBins={}
def binAngle(uid,resType,taName,val):
    global angleBins
    #print "taName:",taName
    if uid not in angleBins:
        angleBins[uid] = {}
        pass
    #print "Arguments for getBin()",
    #print resType,taName,val
    (binName,binMean,branchVal) = getBin(resType,taName,val)
    #print "getBin Return",
    #print binName,binMean,branchVal
    while 1:
        try:
            tabin = angleBins[uid][taName]
        except KeyError:
            angleBins[uid][taName] = {}
            tabin = angleBins[uid][taName]
            pass
        try:
            tabin[binName].append(branchVal)
            break
        except KeyError:
            tabin[binName] = []
            pass
        pass
    StateList.append(str(binName)+" ")
    return

from math import pi
global torStateCnt
torStateCnt={}
global torStateIndx
torStateIndx={}

def binAngles(filename):
    from pdbTool import PDBTool
    pdbFile = PDBTool(filename)
    pdbFile.read()
    afile = open(pdbFile.filename()+".angles",'w')
    #FIX: this assumes only a single segid
    for atom in AtomSel("tag"):
        resid = atom.residueNum()
        segid = atom.segmentName()
        resType = atom.residueName()
        uid = "%4s %4d %4s" % (segid,resid,resType)
        #print >> afile, uid,
        afile.write(uid)
        #print "getAngles()",getAngles(resType,resid,segid);
        angleName=[]
        global StateList
        StateList=[]
        if uid not in torStateCnt:
           torStateCnt[uid] = {}
           pass
        for (aname,atomSel) in getAngles(resType,resid,segid):
            angleName.append(aname)
            val=0
            try:
                val=Dihedral(*atomSel).value() * 180./pi
                afile.write(' %7.2f '%val)
                #print >> afile,  "%7.2f" % val,
                binAngle(uid,resType,aname,val)
            except IndexError:
                #print >> afile, "[%5s]"%aname,
                afile.write('[%5s] '%aname)
                pass
            
            #print >> afile
            pass
        afile.write('\n')
        StateStr=''.join(StateList)
        if StateStr not in torStateCnt[uid]:
           torStateCnt[uid][StateStr] = 0
           pass
        torStateCnt[uid][StateStr]=torStateCnt[uid][StateStr]+1
    afile.close() 
    return


def getIndex():
    return torStateIndx

def getRotomericStates():
    return torStateCnt

global RotAngleNames
RotAngleName={}

def binAngleState(index): 
    #FIX: this assumes only a single segid
    from atomSel import AtomSel
    for atom in AtomSel("tag"):
        resid = atom.residueNum()
        segid = atom.segmentName()
        resType = atom.residueName()
        uid = "%4s %4d %4s" % (segid,resid,resType)
        #print "uid->",uid
        #print "getAngles()",getAngles(resType,resid,segid);
        angleName=[]
        global StateList
        StateList=[]
        if uid not in torStateCnt:
           torStateCnt[uid] = {}
           pass
        if uid not in torStateIndx:
           torStateIndx[uid] = {}
           pass
        for (aname,atomSel) in getAngles(resType,resid,segid):
            if (aname == 'phi' or aname == 'psi' or aname == 'omega'):
                continue
                pass
            angleName.append(aname)
            val=0
            try:
                val=Dihedral(*atomSel).value() * 180./pi
                binAngle(uid,resType,aname,val)
            except IndexError:
                #print "IndexError[%5s]"%aname
                pass
            pass
        if uid not in RotAngleName:
           RotAngleName[uid]=angleName
           pass
        
        StateStr=''.join(StateList)
        #print StateStr
        if StateStr not in torStateCnt[uid]:
           torStateCnt[uid][StateStr] = 0
           torStateIndx[uid][StateStr]=[]
           pass
        torStateCnt[uid][StateStr]=torStateCnt[uid][StateStr]+1
        torStateIndx[uid][StateStr].append(index)
    return

    
def getRotomericStateNames():
    return RotAngleName

def RotomericStates(atomCoordList):
    import xplor
    sim = xplor.simulation
    index=0
    for atmCoords in atomCoordList:
        sim.setAtomPosArr(atmCoords)
        binAngleState(index)
        index=index+1
        pass
    return 

def average(l):
    ret=0.
    if not l: return ret
    
    for el in l:
        ret += el
        pass
    ret /= len(l)
    return ret
from math import sin,cos,atan,pi,sqrt,log,pow

def CStat(l):
   """
   Added by Robin A Thottungal
   Calculates the average of Circular Data(in this case angles)
   Reference Book: Statistical analysis of circular data by
                   N. I. Fisher [ISBN:0521350182] 
   """
   ret=0.
   sum_sin=0.
   sum_cos=0.
   radius=0.
   if not l: return ret
   for elm in l:
      rad=elm*pi/180.
      sum_sin+=sin(rad)
      sum_cos+=cos(rad)
      pass
   radius=sqrt(pow(sum_sin,2)+pow(sum_cos,2))/len(l)
   V=1-radius
   sd=sqrt(-2*log(1-V))
   tol=1e-10
   if sum_sin>tol and sum_cos>tol:     #1st Quad
      ave=atan(sum_sin/sum_cos) 
   elif sum_sin>tol and sum_cos<-tol:   #2nd Quad
      ave=atan(sum_sin/sum_cos)+pi 
   elif sum_sin<-tol and sum_cos<-tol:   #3rd Quad
      ave=atan(sum_sin/sum_cos)-pi  
   elif sum_sin>tol and sum_cos<-tol:   #4th Quad 
      ave=atan(sum_sin/sum_cos)
   else:
       ave=0
       print("Warning: Standard Deviation is not defined")
       if abs(sum_cos)<tol:
           print("Warning:Sum of cosine is zero. Average not defined")
   return (ave*180/pi,sd*180/pi)
       

def deviation(l):
    from math import sqrt
    ret=0.
    if not l: return ret

    ave = average(l)

    for el in l:
        ret += (el-ave)**2
        pass
    ret /= len(l)
    return sqrt(ret)

#### Written by CDS #####
# error is < 1

def setTorsions(angles=[],
                selection="all",
                ivm=None,
                verbose=False,
                simulation=0):
    """
    Set selected torsion angles to specified values.
    
    angles is a list of one or more tuples of the form (sel, val), where val is
    the desired dihedral value (in degrees) and sel is a four-membered tuple,
    each member a selection string identifying one of the four atoms involved
    in the torsion angle.  

    The optional ivm argument specifies an <m ivm>.IVM object referred to for
    torsion angle definitions; if omitted, an <m ivm>.IVM object is created to
    allow all standard torsion angles.  One would create an IVM object
    externally primarily for optimization purposes.

    This function will set all torsion angles specified by the angles list that
    are additionally included in the argument selection (a selection string).

    The simulation argument can be used to specify a <m simulation>.Simulation
    other than the current one.
    """

    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection,simulation)
    simulation = selection.simulation()

    from atomSel import notSelection, intersection
    if ivm==None:
        from ivm import IVM
        import protocol
        ivm=IVM(simulation)
        ivm.group(notSelection(selection))
        protocol.torsionTopology(ivm,flexRiboseRing="not all")
    else:
        ivm.init()
        pass

    from math import pi
    deg2rad = pi/180.
    rad2deg = 180./pi

    from dihedral import Dihedral

    fullCount=0
    for (sel,val) in angles:
        count=0
        val *= deg2rad
##        if type(sel)==type("string"):
##            sel = map(lambda name: "name %s" % name, angleLookup[sel])
##            pass
        atoms1 = intersection(AtomSel(sel[1],simulation),
                              selection)

                      #added "name "+sel[] to fix the problem by RAT
        for node in ivm.nodeList():
            if node.type()=='torsion':
                
                id1=-1
                id2=-1
                if atoms1.containsIndex( node.atoms()[0] ):
                    id1 = node.atoms()[0]
                    id2 = node.parentAtom()
                elif atoms1.containsIndex( node.parentAtom() ):
                    id2 = node.atoms()[0]
                    id1 = node.parentAtom()
                    pass

                  
                if id1<0: continue
                
                try:                                                
                    d = Dihedral(AtomSel("%s and bondedto id %d" %
                                         (sel[0],id1+1), simulation),
                                 AtomSel("id %d"                 %
                                         (id1+1),        simulation),
                                 AtomSel("id %d"                 %
                                         (id2+1),        simulation),
                                 AtomSel("%s and bondedto id %d" %
                                         (sel[3],id2+1), simulation))
                    
                    internalCoords=ivm.pos()
                    index=node.startIndex()

                    internalCoords[index] += val-d.value()
                    count += 1
                    
                    #print 'old val: ', d.value()*rad2deg
                    #print 'target:  ', val*rad2deg
                    target=val*rad2deg
                    ivm.setPos(internalCoords)           #for testing
#                    print 'new val: ', d.value()*rad2deg
                except IndexError as obj:
                    #print "IndexError", obj
                    pass
                pass
            pass
        if count==0:
            print('setTorsions: WARNING: no angle found for')
            print('    ',sel)
            pass
        fullCount+= count
        pass
    if verbose:
        print('set %d torsion angle(s)' % fullCount)
        pass
    pass
   
##angleLookup={}
##angleLookup['psi'] = "N", "CA", "C", "N"
##angleLookup['omega'] = "CA", "C", "N", "CA"
    
    
def setTorsionsFromTable(table):
    """Set dihedral angles to the nominal values specified by table.

    The input table specifies the dihedral angles to be modified and the target
    values (the function <m torsionTools>.setTorsion() sets the new angles).
    table contains XPLOR <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/cdih_syntax.html
    RESTraints DIHEdral statements> of the type associated with the
    <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/cdih.html CDIH term>.
    table can be:
       1) a string with such statements
       2) a string with the name of a file containing the statements
       3) a sequence with any combination of the above two

    """

    from parseTools import findNested


    import os

    string = ''  # string to contain all dihedral restraints
    
    if type(table) == str: table = [table]
    for group in table:
        if os.access(group, os.F_OK):  # group is a filename
            infile = open(group, 'r')
            string += infile.read()
            infile.close()
        else:  # assume group is a str
            string += group
    

    lines = string.split('\n')

    #strip ! comments
    cnt=0
    import re
    while cnt<len(lines):
        lines[cnt] = re.sub(r"!.*$","",lines[cnt])
        cnt+=1
        pass
    
    buf = '\n'.join(lines)

    #strip {} comments
    import potUtils
    buf = potUtils.stripBracketedComments(buf)
    
    list=[]

    #process assi statements
    while True:
        m=re.search(r"^[ \t]*assi",buf,re.IGNORECASE|re.MULTILINE)

        if m==None:
            break

        sels=[]
        indx = m.start(0)
        buf = buf[indx:]
        for i in range(4):
            start=buf.find('(')
            end = findNested('(',')',start,buf,0)
            sels.append( buf[start+1:end] )
            buf = buf[end:]
            pass

        val = float(buf[1:].split()[1])


        list.append( (sels,val) )
        pass

    setTorsions(list)


    if len(list)==0:
        raise Exception("no restraints found in restraint table: %s" % table)
    return

class RotomericStats:
    
    def __init__(s):
        
        """
        A class to analyze the rotomeric states of a given ensemble
        -----dictionarie used in the class----------
        torStatCnt  :dict of the form X[UID][State]=count
                     where X[UID] is a dictionary
                     and X[UID][State] another dictionary
        RotAngleName:dict of the from X[UID]=angleName which is a List of tuple
                     of the form [('CHI1',('N', 'CA', 'CB', 'OG')),
                                  ('CHI2',('CA','CB', 'OG', 'HG'))]
        torStateIndex:dict of the form X[UID][State]=Index which is a
                      list of integers that stores the structure index
                      value in the atomCoordinates List.
        scTorsionAtoms: a dict of the form
                       X['SER'] = [('CHI1',('N', 'CA', 'CB', 'OG')),
                                   ('CHI2',('CA','CB', 'OG', 'HG'))]           
       """
        s.torStateCnt={}
        s.RotAngleName={}
        s.torStateIndex={}
        s.angleBins={}
        s.StateList=[]
        s.torStateIndx={}
        
       
    def getAngles(s,resType,resid,segid):
        """
        Given a residue, the method returns the rotomerics state
        definition for the residue. 
        """
        segSel = 'and segid "%s"' % segid
        angles = [("Phi",  ("name C  and resid %d %s" % (resid-1,segSel),
                            "name N  and resid %d %s" % (resid,segSel),
                            "name CA and resid %d %s" % (resid,segSel),
                            "name C  and resid %d %s" % (resid,segSel) )),
                  ("Psi",  ("name N  and resid %d %s" % (resid,segSel),
                            "name CA and resid %d %s" % (resid,segSel),
                            "name C  and resid %d %s" % (resid,segSel),
                            "name N  and resid %d %s" % (resid+1,segSel))),
                  ("omega",("name CA and resid %d %s" % (resid,segSel),
                            "name C  and resid %d %s" % (resid,segSel),
                            "name N  and resid %d %s" % (resid+1,segSel),
                            "name CA and resid %d %s" % (resid+1,segSel)))]
        try:
            for (name,atoms) in scTorsionAtoms[resType]:
                angles.append((name,
                               ("name %s and resid %d %s" % (atoms[0],
                                                             resid,segSel),
                                "name %s and resid %d %s" % (atoms[1],
                                                             resid,segSel),
                                "name %s and resid %d %s" % (atoms[2],
                                                             resid,segSel),
                                "name %s and resid %d %s" % (atoms[3],
                                                             resid,segSel))))
                pass
            pass
        except KeyError:
            print("unsupported residue type:", resType)
            pass
        
        return angles

    def getBranch(s,theta1,theta2):
        """return the branch number n such that n*2*PI+theta1-theta2 is
        minimal"""
        return round((theta2-theta1)/360)

    def getBin(s,resType,taname,val):
        """
        Given the residue type and torsion angle name, the method finds out
        the states in which the torsion angle can exit and return the angle
        that is closest to the val.
        """
        closest = ('single',0,val)
        bins=[]
        if taname.startswith('CHI'):
            bins = (('180',180), ('+60',60), ('-60',-60))
            pass
        try:
            bins = binData[resType+'_'+taname]
        except KeyError:
            pass
    
        if bins:
            closest = ('',361,361)
            for (name,mean) in bins:
                val += 360*s.getBranch(val,mean)
                dist=abs(val-mean)
                if dist<closest[1]:
                    closest = (name,dist,val)
                    pass
                #print 'branch',mean,val,getBranch(val,mean), closest
                #print closest
                pass
            pass
        
        return closest
        
    def binAngle(s,uid,resType,taName,val):
        """
        creates a dict named angleBins
        angleBins[uid][torsionangleName] : contains the angle value 
                 corresponding to the torsionangleName
                                           
        """
        #print "taName:",taName
        if uid not in s.angleBins:
            s.angleBins[uid] = {}
            pass
        #print "Arguments for getBin()",
        #print resType,taName,val
        (binName,binMean,branchVal) = s.getBin(resType,taName,val)
        #print "getBin Return",
        #print binName,binMean,branchVal
        while 1:
            try:
                tabin = s.angleBins[uid][taName]
            except KeyError:
                s.angleBins[uid][taName] = {}
                tabin = s.angleBins[uid][taName]
                pass
            try:
                tabin[binName].append(branchVal)
                break
            except KeyError:
                tabin[binName] = []
                pass
            pass
        s.StateList.append(str(binName)+" ")
        return


    def binAngles(s,angleFilename=None):
        """
        Takes in a set of pdb filename and analyze the
        rotomerics state of each pdb file. The result is
        written in the file name [filename].angles

        FIX: currently angleFilename must be specified.
        """
        afile = open(angleFilename,'w') if angleFilename else None
        #FIX: this assumes only a single segid
        for atom in AtomSel("tag"):
            resid = atom.residueNum()
            segid = atom.segmentName()
            resType = atom.residueName()
            uid = "%4s %4d %4s" % (segid,resid,resType)
            #print >> afile, uid,
            afile.write(uid)
            #print "getAngles()",getAngles(resType,resid,segid);
            s.angleName=[]
            s.StateList=[]
            if uid not in s.torStateCnt:
                s.torStateCnt[uid] = {}
                pass
            for (aname,atomSel) in s.getAngles(resType,resid,segid):
                s.angleName.append(aname)
                val=0
                try:
                    val=Dihedral(*atomSel).value() * 180./pi
                    afile.write(' %7.2f '%val)
                    #print >> afile,  "%7.2f" % val,
                    s.binAngle(uid,resType,aname,val)
                except IndexError:
                    #print >> afile, "[%5s]"%aname,
                    afile.write('[%5s] '%aname)
                    pass
            
                 #print >> afile
                pass
            afile.write('\n')
            StateStr=''.join(s.StateList)
            if StateStr not in s.torStateCnt[uid]:
               s.torStateCnt[uid][StateStr] = 0
               pass
            s.torStateCnt[uid][StateStr]=s.torStateCnt[uid][StateStr]+1
        afile.close() 
        return

    def binAngleState(s,index):
        """
        creates two dict:
           torStateCnt: dict of the form X[UID][State]=count
                        where X[UID] is a dictionary
                        and X[UID][State] another dictionary
           torStateIndx:dict of the form X[UID][State]=Index which is a
                        list of integers that stores the structure index
                        value in the atomCoordinates List. 
        
        """
        #FIX: this assumes only a single segid
        from atomSel import AtomSel
        for atom in AtomSel("tag"):
            resid = atom.residueNum()
            segid = atom.segmentName()
            resType = atom.residueName()
            uid = "%4s %4d %4s" % (segid,resid,resType)
            #print "uid->",uid
            #print "getAngles()",getAngles(resType,resid,segid);
            s.angleName=[]
            s.StateList=[]
            if uid not in s.torStateCnt:
                s.torStateCnt[uid] = {}
                pass
            if uid not in s.torStateIndx:
                s.torStateIndx[uid] = {}
                pass
            for (aname,atomSel) in s.getAngles(resType,resid,segid):
                if (aname == 'Phi' or aname == 'Psi' or aname == 'omega'):
                    continue
                    pass
                s.angleName.append(aname)
                val=0
                try:
                    val=Dihedral(*atomSel).value() * 180./pi
                    s.binAngle(uid,resType,aname,val)
                except IndexError:
                    #print "IndexError[%5s]"%aname
                    pass
                pass
            if uid not in s.RotAngleName:
                s.RotAngleName[uid]=s.angleName
                pass
        
            StateStr=''.join(s.StateList)
            #print StateStr
            if StateStr not in s.torStateCnt[uid]:
                s.torStateCnt[uid][StateStr] = 0
                s.torStateIndx[uid][StateStr]=[]
                pass
            s.torStateCnt[uid][StateStr]=s.torStateCnt[uid][StateStr]+1
            s.torStateIndx[uid][StateStr].append(index)
        return

    def RotomericStates(s,atomCoordList):
        import xplor
        sim = xplor.simulation
        index=0
        for atmCoords in atomCoordList:
            sim.setAtomPosArr(atmCoords)
            s.binAngleState(index)
            index=index+1
            pass
        return 
  


    def circularStat(s,l):
        """
        Added by Robin A Thottungal
        Calculates the average of Circular Data(in this case angles)
        Reference Book: Statistical analysis of circular data by
        N. I. Fisher [ISBN:0521350182] 
        """
        from math import sin,cos,atan,pi,sqrt,log,pow
        ret=0.
        sum_sin=0.
        sum_cos=0.
        radius=0.
        if not l: return ret
        for elm in l:
            rad=elm*pi/180.
            sum_sin+=sin(rad)
            sum_cos+=cos(rad)
            pass
        #print sum_sin,sum_cos
        radius=sqrt(pow(sum_sin,2)+pow(sum_cos,2))/len(l)
        V=1-radius
        try:
            sd=sqrt(-2*log(1-V))
        except ValueError:
            sd=0
            pass
            
        tol=1e-10
        if sum_sin>tol and sum_cos>tol:     #1st Quad
            ave=atan(sum_sin/sum_cos) 
        elif sum_sin>tol and sum_cos<-tol:   #2nd Quad
            ave=atan(sum_sin/sum_cos)+pi 
        elif sum_sin<-tol and sum_cos<-tol:   #3rd Quad
            ave=atan(sum_sin/sum_cos)-pi  
        elif sum_sin<-tol and sum_cos>tol:   #4th Quad 
            ave=atan(sum_sin/sum_cos)
        else:
            ave=0
            #print "Warning: Standard Deviation is not defined"
            if abs(sum_cos)<tol:
                print("Warning:Sum of cosine is zero. Average not defined")
        return (ave*180/pi,sd*180/pi)

    # linearAverage / linearDeviation does not give correct answer
    # for circular data always. Use circularStat() instead..
    def linearAverage(s,l):
        ret=0.
        if not l: return ret
    
        for el in l:
            ret += el
            pass
        ret /= len(l)
        return ret   
    def linearDeviation(s,l):
        from math import sqrt
        ret=0.
        if not l: return ret

        ave = average(l)

        for el in l:
            ret += (el-ave)**2
            pass
        ret /= len(l)
        return sqrt(ret)

def genXplorRestraint(name,residueNum,value,width,
                      segmentName=None,
                      checkAtoms=False):
    """Generate one XPLOR torsion angle CDIH restraint, given its name,
    a residue number, the value and width. The segment name may also optionally
    be specified.

    Supported names are currently 'phi' and 'psi'. Capitalization is ignored.

    The return value is the restraint string, or an empy string, if one or more
    of the specified atoms is missing. If the checkAtoms argument is set to True,
    an exception will be thrown in the case of missing atoms.
    """

    residM=residueNum-1
    residP=residueNum+1

    name = name.lower()
    if name=='phi':
        atomsResids=list(zip("C N CA C".split(),
                        (residM,residueNum,residueNum,residueNum)))
    elif name=='psi':
        atomsResids=list(zip("N CA C N".split(),
                        (residueNum,residueNum,residueNum,residP)))
    else:
        raise Exception("unknown torsion angle name")
        pass
    
    ret='assign\n'

    segStr = 'segid "'+segmentName+'" and' if segmentName else ""
    for name,resid in atomsResids:
        atomSelStr = '%s resid %3d and name %4s' % (segStr,resid,name)
        if checkAtoms:
            if len(AtomSel(atomSelStr))!=1:
                if checkAtoms!="NoException":
                    raise Exception("no such atom: " + atomSelStr)
                return ""
            pass
        ret += '  (%s)\n' % atomSelStr
        pass
    ret += '  1.0 %.2f %.2f 2\n' % (value,width)
    return ret

    
