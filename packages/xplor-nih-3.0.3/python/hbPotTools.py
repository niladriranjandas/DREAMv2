"""
Helper functions for creating and analyzing <m hbPot>.HBPot energy terms.
"""

# Location of database within Xplor-NIH directory.
# (Added HBPOT environment variable in ../bin/xplor.in)
from os import environ
dbDir = environ.get('HBPOT','.')

from trace import notrace_decorate
@notrace_decorate
def create_HBPot(name,
                 selection="not PSEUDO",
                 configFilename="da.dat",
                 threshold='99.95',
                 verbose=1):
    """
    create and return a <m hbPot>.HBPot energy term with instanceName name, and
    covering atoms specified in selection.

    Arguments: name
               selection="not PSEUDO",
               configFilename="da.dat",
               threshold='99.95',
               verbose=1

       selection - the term can be applied to a subset of atoms by specifying
                   an alternate value.
       configFilename - file which specifies donors and acceptors for selected
                        residues and a mapping of all pairs to a
                        <m hbSurf>.HBSurf energy surface.
       threshold      - percentage threshold used to determine whether a
                        hydrogen bond is violated. If the energy for a
                        particular geometry is larger than that corresponding
                        to threshold percent of the input datapoints, it is
                        considered to be violated. Valid values are currently
                        '98', '99' and '99.95'.
     
    """

    (aList,dList,mappings) = readConfigFile(configFilename)
    
    if verbose:
        print("found %d acceptor definitions" % len(aList))
        print("found %d donor definitions" % len(dList))
        print("found %d mappings" % len(mappings))
        pass

    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)

    from hbPot import HBPot
    ret = HBPot(name,aList,dList,mappings,
                selection,threshold,verbose>1)

    return ret

def readConfigFile(configFilename):
    """
    Read file defining hydrogen bond acceptors, donors and mappings of
    acceptor/donor pairs to database surface. Acceptor lines are of the format

        acceptor aName sel1 [sel2]

    where aName is an arbitrary label for the acceptor, sel1 is an
    atom selection string (<m atomSelLang>), and sel2 is an optional second
    selection string which defines the bound-to atom. The second selection string
    is necessary if the acceptor is bound to more than one atom. If sel2 is
    specified, sel1 must be enclosed in parentheses.

    Donor lines are of the form

        donor dName sel

    where dName is a label for the donor and sel is the selection string.

    Mappings from acceptor/donor pairs to database source are of the form

        mapping aName:dName deltaResid filename

    where aName and dName are one of the acceptor or donor labels, respectively,
    or the "*" character to match any. The deltaResid argument specifies a
    difference in primary sequence of the donor and acceptor residues, or
    the "*" character to match any. The filename argument specifies a file for
    the hydrogen bonding potential of mean force relative to the FIX...

      
    """


    global dbDir
    
    altDir="."
    import os
    if os.path.exists(configFilename):
        altDir = os.path.dirname(configFilename)
    else:
        configFilename=os.path.join(dbDir,configFilename)
        if not os.path.exists(configFilename):
            raise Exception("Could not find config file in . or " +
                            dbDir)
            pass
        pass

    #read acceptor/donor file
    lines = open(configFilename).readlines()
    aLines = [line.strip() for line in lines
              if line.strip().startswith("acceptor")]
    aList=[]
    from parseTools import findNested
    for line in aLines:
        a, typeStr, sels = line.split(None,2)
        if sels.startswith('('):
            end = findNested('(',')',1,sels)
            sel = sels[1:end]
            sels = sels[end+1:].strip()
            if sels.startswith('('):
                end = findNested('(',')',1,sels)
                nextSel = sels[1:end]
                aList.append((typeStr,sel,nextSel))
                continue
                pass
        else:
            sel = sels
            pass
        aList.append((typeStr,sel))
        pass
    dLines = [line.strip() for line in lines
              if line.strip().startswith("donor")]
    dList=[]
    for line in dLines:
        d, typeStr, sel = line.split(None,2)
        dList.append( (typeStr,sel) )
        pass
    mLines = [line.strip() for line in lines
              if line.strip().startswith("mapping")]
    
    mappings=[]
    #fix: move this logic to C++ layer. change type of deltaResid to String
    for line in mLines:
        fields=line.split()
        m, typePair, deltaResid, filename = fields[:4]
        aType,dType=typePair.split(":")
        deltaResid = int(deltaResid) if deltaResid!='*' else 7
        import os
        if os.path.exists(os.path.join(altDir,filename)):
            filename = os.path.join(altDir,filename)
        else:
            filename = os.path.join(dbDir,filename)
            pass
        if not os.path.exists(filename):
            raise Exception("Could not find HBSurf data for typePair " +
                            typePair + " in " + altDir + " or " + dbDir)
        scale=float(fields[4]) if len(fields)>4 else 1.
        mappings.append( (typePair,filename,deltaResid,scale) )
        pass

    return aList,dList,mappings

def hBondInfo(term):
    """
    Return a pretty-printed string with information on each detected hydrogen
    bond.
    """

    ret = \
  "({:^19}) <-- ({:^19}) {:7} {:^5} {:<6} {:^6} {:^4} {}\n".format("Acceptor",
                                                                   "Donor",
                                                                   "Energy",
                                                                   "r","theta",
                                                                   "phi",
                                                                   "viol",
                                                                   "db filename")
    for pair in term.hbonds():
        a = term.acceptors[pair.a]
        d = term.donors[pair.b]
        aAtom=a.atom
        dAtom=d.atom
        energy = term.energy(a,d)
        from os import path
        mapStr = path.split(term.getMapString(a,d))[-1]
        geom = term.getGeom(a,d)
        violated = term.violated(a,d)
        
        ret += "(%s) <-- (%s) %7.3f %5.2f %6.2f %6.2f %4s %s\n" %(
            aAtom.string(),
            dAtom.string(),
            energy,
            geom.r,
            geom.theta,
            geom.phi,
            '*' if violated else '',
            mapStr)
#                                                              info.r,
#                                                              info.theta,
#                                                              info.phi,
#                                                              info.dbName)
        pass
    return ret

def numHBonds(term):
    return len(term.hbonds())
               
    
def analyze(potlist):
    """Perform analysis of HBPot terms.

    """
    from simulationTools import getPotTerms
    terms = getPotTerms(potlist, 'HBPot')

    if not terms: return ''  # this goes to pdb file header

    instanceNames = [term.instanceName() for term in terms]

    ret = "%10s  %8s\n" % ("","# hbonds")
    summary = ""
    for term in terms:
        summary += "\n  HBPot term: %s\n\n" % term.instanceName()
        summary += hBondInfo(term) + "\n";
        summary += term.info()
        summary += "database directory: " + dbDir + '\n'
        ret += "{:<10}  {:^8d}\n".format(term.instanceName(), numHBonds(term))
#        ret += "%10s  %8d\n" % (term.instanceName(), numHBonds(term))
        pass
    
    print(summary)      # this goes to .viols file
    return ret  # this goes to pdb file header


from simulationTools import registerTerm, registerExtraStats
registerTerm(analyze, 'HBPot terms', 'HBPot',
r"""
For each <m hbPot>.HBPot term, print the number of detected hydrogen bonds.
""")
    
registerExtraStats("HBPot","# hbonds",numHBonds)


                    


            
            

        
    

                     
                    
