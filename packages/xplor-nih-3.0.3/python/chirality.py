"""
Use the  Cahn-Ingold-Prelog priority rules to determine absolute
R or S chiral configuration.
"""

#
# bond length thresholds in decreasing order
#
# bonds longer than the first number are single bonds, those between
# the first and second, double, etc.
#
bondOrderByLength={ "cc" : (1.41,1.25),
                    "co" : (1.34,),
                    }



def getBondParams(simulation):
    """Return full matrix mapping chemical type pairs to bond parameters
    """
    from xplorSimulation import getXplorSimulation
    from simulationTools import mktemp
    xSim = getXplorSimulation(simulation)

    outputState=xSim.disableOutput()

    parFilename=mktemp('param-output')

    xSim.fastCommand("write param output=%s end" % parFilename)
    xSim.enableOutput(outputState)

    paramOutput = [line.strip() for line in open(parFilename).readlines()]
    import os
    os.unlink(parFilename)
    bondLines=[line for line in paramOutput if line.startswith('BOND')]

#    print 'read %d bond params' % len(bondLines)


    from potUtils import stripBracketedComments
    bondParams = [stripBracketedComments(line).split()[1:5]
                  for line in bondLines]

    ret={}
    for type1,type2,fc,l0 in bondParams:
        fc = float(fc)
        l0=float(l0)
        
        ret[type1+' '+type2] = (fc,l0)
        ret[type2+' '+type1] = (fc,l0)
        pass

    return ret

def readBondOrders(filename):
    """
    Read entries from a bond order file. Entries have the form:

    double name1 name2
    triple name1 name2

    where name1 and name2 are atom names. By default all bonds are
    single bonds, so they do not need to be specified.
    """
    lines=open(filename).readlines()
    ret={}
    for line in [line.strip() for line in lines]:
        if line.startswith("#"):
            continue
        elif line.lower().startswith("double"):
            atom1,atom2 = (s.lower() for s in line.split()[1:3])
            ret[atom1+' '+atom2] = 2
            ret[atom2+' '+atom1] = 2
        elif line.lower().startswith("triple"):
            atom1,atom2 = (s.lower() for s in line.split()[1:3])
            ret[atom1+' '+atom2] = 3
            ret[atom2+' '+atom1] = 3
        else:
            raise Exception("Error in bond order entry: " + line)
        pass
#    print "readBondOrders: ", ret
    return ret



def sign(v):
    if v>=0.:
        return 1
    return -1

def absoluteChirality(atomName,bondOrderFile=None,
                      verbose=False
                      ):
    """Determine S or R absolute configuration. Return value is a tuple of
    center name, tuple of ordered bond atom names, and string representing
    absolute chirality.
    """

    centerSel="name %s" % atomName
    from atomSel import AtomSel
    center=AtomSel(centerSel)
    centerAtom=center[0]

    bondParams={}
    if bondOrderFile:
        bondOrders = readBondOrders(bondOrderFile)
    else:
        bondParams = getBondParams(centerAtom.simulation())
        pass

    bonded=AtomSel("bondedto %s" % centerSel)

    def priority(atom):
        """atom priority goes with atomic weight
        """
        ret="error"

        if atom.atomName().startswith('H'):
            ret=1
        elif atom.atomName().startswith('C'):
            ret=10
        elif atom.atomName().startswith('N'):
            ret=100
        elif atom.atomName().startswith('O'):
            ret=1000
        else:
            ret=10000 #assume a heavier element
            pass
        return ret
    #    for neighbor in AtomSel("bondedto id %d" %atom.id()):
    #        if neighbor==prevAtom:
    #            continue
    def bondOrderByDistance(atom1,atom2,bondParams):
        """Deduce bond order by atom-atom distance and atom types.
        """
        from vec3 import norm
        key = atom1.chemType() + ' ' + atom2.chemType()
        if key in list(bondParams.keys()):
            blen = bondParams[key][1]
        else:
            blen = norm(atom1.pos()-atom2.pos())
            pass
        
        ret=1
        key=''.join(sorted(atom1.atomName()[0].lower()+
                          atom2.atomName()[0].lower()))
        if True:
            if key in list(bondOrderByLength.keys()):
                lengths=bondOrderByLength[key]
                for index,length in enumerate(lengths):
                    if blen>lengths[index]:
                        break
                    ret += 1
                    pass
                pass
            pass
        else:
            #old code with hard-coded distance thresholds: remove
            types = [atom1.atomName()[0].lower(),
                     atom2.atomName()[0].lower()]
    #        print 'bond length', blen, key
    
            if types[0] == "c" and types[1] == "c":
                # deduced from
                # http://chemistry.osu.edu/~woodward/ch121/ch8_bondorder.htm
                if blen<1.25:
                    ret=3
                elif blen<1.41:
                    ret=2
                    pass
                pass
            elif 'c' in types and 'o' in types:
                #deduced from
                # http://www.wiredchemist.com/chemistry/data/bond_energies_lengths.html
                if blen<1.34:
                    ret=2
                    pass
                pass
            pass
        
        if ret>1 and verbose:
            print('bondOrderByDistance: ', \
                  atom1.atomName(),atom2.atomName(), ret)

        return ret

    def determineBondOrder(atom1,atom2):
        if bondOrderFile:
            key = (atom1.atomName() + ' ' + atom2.atomName()).lower()
#            print key, bondOrders.keys()
            if key in list(bondOrders.keys()):
                if verbose:
                    print('determineBondOrder:', key, bondOrders[key])
                    pass
                return bondOrders[key]
            else:
                return 1
            pass
        else: # use bond length from parameters
            return bondOrderByDistance(atom1,atom2,bondParams)
        return
        

    if len(bonded)!=4:
        raise Exception("atom %s has %d bonds. Cannot determine chirality" % (
            center.string(), len(bonded)))



    class AtomNode:
        def __init__(s,atom,parentNode=None):
            s.atom=atom
            s.parentNode=parentNode
            s.depth=parentNode.depth+1 if parentNode else 0
            #print 'ctor',atom.string(),s.depth
            s.children=[]
            return
        def addChildren(s,depth):
            " Add children at specified depth."
            #print 'addChildren', s.depth,depth, len(s.children)
            if s.depth>depth:
                return
            if not s.children :
                bonded = AtomSel("bondedto index %d" % s.atom.index())
                for atom in bonded:
                    if atom!=s.parentNode.atom:
                        s.children.append( AtomNode(atom,s) )
                        pass
                    pass
                pass
            else:
                for child in s.children:
                    child.addChildren(depth)
                    pass
                pass
            return
        def getAtomsAtDepth(s,depth):
            if s.depth==depth:
                return s.children
            ret = []
            if s.depth>depth:
                return ret
            for child in s.children:
                ret += child.getAtomsAtDepth(depth)
                pass
            return ret
        pass
    
    #
    # this algorithm should be simple enough to grok quickly.
    # A tree is constructed deep enough so that all 4 branches can be
    # distinguished.
    #


    centerNode=AtomNode(centerAtom)
    trees = [AtomNode(atom,centerNode) for atom in bonded]
    priorities = [[priority(atom)] for atom in bonded]

    level=0
    maxLevel=100
    while True:
        done = True
        computeNextLevel=[] # atoms which need to go to next level
        #print 'level:',level
        for i,atomi in enumerate(bonded):
            for j,atomi in enumerate(bonded):
                if j<=i:
                    continue
                if (len(priorities[i])>level and
                    len(priorities[j])>level and
                    priorities[i][level]==priorities[j][level]):
                    #catch the case where we're at the end of branches
                    done = False
                    if not i in computeNextLevel:
                        computeNextLevel.append(i)
                    if not j in computeNextLevel:
                        computeNextLevel.append(j)
                        pass
                    pass
                pass
            pass

        if done: break

        childrenAdded=0
        level += 1
        if level>maxLevel:
            raise Exception("Too many recursion levels. " +
                            "Achiral center with loops?")
        for i in computeNextLevel:
            trees[i].addChildren(level)

            curScore=0
            for child in trees[i].getAtomsAtDepth(level):
                childrenAdded += 1
                bo = determineBondOrder(child.atom,child.parentNode.atom)
                curScore += bo * priority(child.atom)
                pass
            priorities[i].append( curScore )
            pass
        if childrenAdded==0:
            raise Exception("Degeneracy detected. Is this a chiral center?")
        
        pass
    #print priorities

        
            
    bondedPriorities = list(zip(bonded,priorities))
    def cmp(a, b):
        return (a > b) - (a < b)
    def sorter(a,b):
        aAtom,aPriorities=a
        bAtom,bPriorities=b
        level = 0 
        while aPriorities[level]==bPriorities[level]:
            level += 1
            pass
        return cmp(aPriorities[level],bPriorities[level])
    from functools import cmp_to_key
    bondedPriorities.sort( key=cmp_to_key(sorter) )
        
    sortedAtoms = [atom for atom,p in bondedPriorities]
    
    ret = [centerAtom.atomName()]
    ret.append( [atom.atomName() for atom in sortedAtoms] )

    #now can determine chirality
    from vec3 import unitVec, cross, dot
    relativeVec = [unitVec(atom.pos()-centerAtom.pos())
                   for atom in sortedAtoms]

    c1 = dot(cross(relativeVec[3],relativeVec[2]),relativeVec[0])
    c2 = dot(cross(relativeVec[2],relativeVec[1]),relativeVec[0])

    if sign(c1)!=sign(c2):
        print("absoluteChirality: signs don't match: c1=%f, c2=%f" %(c1,c2))
        ret.append('T')
    elif sign(c1)>0.:
        ret.append("R")
    else:
        ret.append("S")
        pass
    return ret

def reportChiralities(centers,bondOderFile=None):
    """ Return a nicely formatted string reporting the chirality of each chiral
    center specified in centers.
    """
    ret=""
    brief="chirality in brief: "
    lines=[]
    for center in centers:
        cname,sortedBoundNames,chirality = absoluteChirality(center,
                                                             bondOderFile)
        ret += "center %-4s" % center
        ret += " [sorted bound atoms: "
        sortedBoundNames = ["%-4s" %n for n in sortedBoundNames]
        ret += ",".join(sortedBoundNames)
        ret += "] " + chirality + '\n'
        brief += chirality
        pass
        
    ret += brief
    return ret

    
    
