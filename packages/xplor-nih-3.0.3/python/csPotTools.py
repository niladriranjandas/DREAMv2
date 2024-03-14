
def create_CSPot(name,
                 restraints="",
                 filename=None,
                 format="plain",
                 saveSet=None,
                 version="1",
                 segid=None,
                 ambiguousGlyHA=False,
                 defaultError=1,
                 useTableErrors=False,
                 selection='not PSEUDO',
                 verbose=False
                 ):

    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)

    from csPot import csPot

    if restraints==None: restraints=""


    if type(restraints)==type('string'):
        if restraints: restraints += '\n'
        if filename:   restraints += open(filename).read()

        from chemShiftTools import convertRestraints
        restraints = convertRestraints(restraints,
                                       format=format,
                                       saveSet=saveSet,
                                       segid=segid,
                                       ambiguousGlyHA=ambiguousGlyHA,
                                       verbose=verbose)
        pass

    if not restraints:
        raise Exception("no chemical shifts read")

    from os import environ

    #FIX: should probbaly not be hardcoded
    camshiftData=environ["XPLOR_DIR"]+"/databases/camshift/"
    params=[
        "1CAorgpairs.par",
        "1HAorgpairs.par",
        "1Norgpairs.par",
        "1Horgpairs.par",
        "1Corgpairs.par",
        "1CBorgpairs.par",
        ]

    shifts=[
        "CAshifts.dat",
        "HAshifts.dat",
        "Nshifts.dat",
        "Hshifts.dat",
        "Cshifts.dat",
        "CBshifts.dat",
        ]
    return csPot(name,
                 restraints,
                 version,
                 selection)
    
                

def analyze(potList):
    "perform analysis of csPot terms and return nicely formatted summary"
          
    ret = ""
    
    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'csPot')
    
    if not potList: return ret
    
    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()
    
    ret+= "%9s  %6s   %4s\n" % \
          (" " , "RMS", "Viols")
        
    rdcs = []
    for name in instanceNames:
        rdc = [x for x in potList if x.instanceName()==name][0]
        rdcs.append(rdc)

        #print cs.showViolations()
        
        #print cs.info()
        #tensor=rdc.tensor
    
        ret += "%-9s  %6.3f    %4d\n" % \
               (name , rdc.rms(), rdc.violations() )
        pass
    
    return ret

