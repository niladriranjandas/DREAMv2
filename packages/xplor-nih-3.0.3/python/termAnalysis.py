"""Tools to aid in reporting on and analyzing <m pot>.Pot terms
"""

#list of modules to be imported when generating complete documentation
# on analysis output
registeringModules="""anglePot ccrPotTools csaPotTools cstMagTools diffPotTools
                      distSymmTools
                      ensWeightsTools globDiffPotTools gyrPotTools h3JNCPot
                      hbPotTools hmxPotTools
                      jCoupPotTools
                      nbTargetPotTools
                      noePotTools orderPotTools
                      posRMSDPotTools
                      planeDistTools posDiffPotTools prePotTools
                      probDistPotTools
                      psolPotTools rdcCorrPotTools
                      rdcPotTools relaxRatioPotTools repelPotTools
                      residueAffPotTools  
                      sardcPotTools saTensorTools shapePotTools
                      solnScatPotTools
                      spartaPotTools torsionDBPotTools
                      varTensorTools xplorPotTools
                      ramaPushPot distSplinePotTools""".split()


def summarizePotentials(potList):
    "get rms, violations summaries from potential terms"
    terms = []

    for pot in potList:
        rms=-1
        viols=-1
        if pot.potName()=='PotList' and len(pot):
            lterms=summarizePotentials(pot)
            rms=0.
            viols=0
            count=0
            for (pName,name,lrms,lviols) in lterms:
                if lrms>=0:
                    count += 1
                    rms += lrms
                    viols += lviols
                    pass
                pass
            if count>0:
                rms /= count
                viols
                pass
            pass
        try:
            rms = pot.rms()
        except AttributeError:
            pass
        try:
            viols = pot.violations()
        except AttributeError:
            pass
        terms.append((pot.potName(),pot.instanceName(),rms,viols))
        pass
    return terms

def getPotTerms(potList,names):
    """return in a list all the potential terms in potList whose potType()
    matches names. names can be a string or a list of strings.
    """
    from simulationTools import potType, flattenPotList
    pl = flattenPotList(potList)
    if type(names)==type('string'): names = (names,)
    ret =[]
    for name in names:
        ret += [x for x in pl if potType(x)==name]
        pass
    return ret 
                

registeredTerms=[]
from collections import namedtuple
from functools import reduce
AnalyzeTermEntry=namedtuple('TermEntry',
                            'func title prefix moduleName doc'.split())

def registerTerm(analyzeFunc,termTitle,termPrefix,termDoc):
    """ register an analysis function to be called by the analyze() function
    [see below].

    One should specify an analysis function which takes a list of potential
    terms as an argument, and two strings to identify the potential term - long
    form, and short form, respectively.

    The termDoc argument documents the output of the returned string - in
    reStructuredText. This string is produced by the headerHelp helper
    program.

    analyzeFunc should return a string containing an analysis summary. A more
    detailed analysis can be printed to stdout.
    """
    moduleName=analyzeFunc.__module__
    if not moduleName in registeringModules:
        print("registerTerm: Warning: module is not cataloged:", moduleName)
    registeredTerms.append(
        AnalyzeTermEntry(analyzeFunc,termTitle,termPrefix,moduleName,termDoc) )
    return

def loadRegisteringModules():
    """Load all modules which call registerTerm- to populate the header help
    information.
    """
    for moduleName in registeringModules:
        exec("import "+moduleName)
        pass
    return

def canonicalTermName(name):
    """ Return the registered term name given a string. Currently, this simply
    makes the term name search case insensitive."""
    for term in registeredTerms:
        if term.prefix.lower() == name.lower():
            return term.prefix
        pass
    return None

def moduleName(termName):
    """ Return module name of associated term"""
    for term in registeredTerms:
        if term.prefix == termName:
            return term.moduleName
        pass
    return None

def getHeaderNames():
    return [termEntry.prefix for termEntry in registeredTerms]


def xplorDocMarkupToRST(text):
    """
    Convert Xplor-NIH-specific documentation markup to reStructuredText.

    Converted tags include
           <\m NAME> - creates a link for Xplor-NIH module named NAME
           <\s NAME> - creates a link for Python module named NAME
           <\l LINK NAME> - creates a link LINK with tag NAME
    """

    def escape(text): return text
    from xplorDoc import pyDocURL
    xplorNIHpyDocURL = "https://nmr.cit.nih.gov/xplor-nih/doc/current/python/ref"

    results = []
    here = 0
    import re
    pattern = re.compile(r'(<l\s+(\S+)\s+([^>]+)>|'
                         r'\<m\s+(\S+\w)\>|'
                         r'\<s\s+(\S+\w)\>)')
    while True:
        match = pattern.search(text, here)
        if not match: break
        start, end = match.span()
        results.append(escape(text[here:start]))

        (all, link, linkName, 
         modName, pyModName) = match.groups()
        if link:
            results.append('`%s <%s>`_' % (linkName, link))
        elif modName:
            results.append('`%s <%s/%s>`_' % (modName,
                                              xplorNIHpyDocURL,
                                              modName+'.html'))
        elif pyModName:
            results.append('`%s <%s/%s>`_' % (pyModName,
                                              pyDocURL,
                                              pyModName+'.html'))
        else:
            results.append(self.namelink(name, []))
            pass
        here = end
        pass
    results.append(escape(text[here:]))
    return ''.join(results)

def genHeaderHelp(term):
    """Return reStructured text string for the specified term.
    """
    
    for termEntry in registeredTerms:
        if termEntry.prefix==term:
            ret = "%s\n"%term 
            ret += "-"*len(term) + '\n'
            ret += xplorDocMarkupToRST(termEntry.doc)
            return ret
        pass
    return ""
        
    
from potList import PotList
def analyze(potList,extraTerms=PotList(),
            outFilename=0):
    """ pretty print appropriate terms from the given PotList and return
    them in a remarks string.
    The optional extraTerms is a PotList of terms which do not
    contribute to the total energy, but which should be analyzed all the same-
    use this for cross-validated potential terms. The potList and extraTerms
    arguments can have type of Pot, PotList, list or tuple.

    If outFilename is specified, the verbose violations info is written to
    the specified file.

    """

    try:
        len(potList)
    except:
        potList = [potList]
        pass

    if type(potList)==type([]) or type(potList)==type(tuple([])):
        terms=potList
        potList=PotList()
        for p in terms:
            potList.append(p)
            pass
        pass
    if type(extraTerms)==type([]) or type(extraTerms)==type(tuple([])):
        terms=extraTerms
        extraTerms=PotList()
        for p in terms:
            extraTerms.append(p)
            pass
        pass
    if extraTerms==None:
        extraTerms=PotList()
        pass

    if outFilename:
        outfile = open(outFilename,"w")
        print("\n  Violation Analysis \n", file=outfile) 
        outfile.close()
        pass

    ret=""
    totViols=0
    totEnergy = potList.calcEnergy()
    reports = potList.energyReports()

    terms = summarizePotentials(potList)

    for term in terms:
        (potType,name,rms,viols) = term
        energy = [x for x in reports if x[0]==name][0][1]
        rmsString=" "*8
        if rms>-1: rmsString = "%8.3f" % rms
        violString=' '*8
        if viols>-1:
            violString = "%8.1f" % viols
            totViols += viols
            pass
        ret += "summary %-10s %10.2f %s %s\n" % (term[1],
                                               energy,rmsString,violString)
        pass

    ret = "summary %-10s %10.2f %8s %8.1f\n" % ("total",
                                              totEnergy,'',totViols) + ret
    ret = "summary    Name       Energy      RMS     Violations\n" + ret
    ret = "-"*60 + '\n' + ret
    
    ret += "-"*60 + '\n'

    if len(extraTerms):
        ret += "\nCross-validated terms:\n"
        extraTerms.calcEnergy()
        reports = extraTerms.energyReports()

        terms = summarizePotentials(extraTerms)
    
    
        for term in terms:
            (potType,name,rms,viols) = term
            energy = [x for x in reports if x[0]==name][0][1]
            rmsString=" "*8
            if rms>-1: rmsString = "%8.3f" % rms
            violString=' '*8
            if viols>-1: violString = "%8.1f" % viols
            ret += "summary %-10s %10.2f %s %s\n" % (term[1],
                                                     energy,
                                                     rmsString,violString)
            pass

# why do I care about this?
#        terms = flattenPotList(extraTerms)
#
#        if len(terms) > len(extraTerms):
#            ret += "summary\n"
#        
#            lineBeg = "summary cross-validated terms: "
#            line=lineBeg
#            for term in terms:
#                if (len(line)+len(term.instanceName()))>71: 
#                    ret += line + '\n'
#                    line=lineBeg
#                    pass
#                line += "%s" % term.instanceName()
#                if term != terms[-1]: line += ", "
#                pass
#            if line != lineBeg: ret += line + '\n'
#            pass

        ret += "-"*60 + '\n'

    if outFilename:
        import sys
        old_stdout = sys.stdout
        sys.stdout = open(outFilename,"a")
        pass

    for analyzeTerm in registeredTerms:
        tmp = analyzeTerm.func(list(potList)+list(extraTerms))
        if tmp:
            ret += ' %s\n' % analyzeTerm.title
            ret += reduce(lambda x,y: x+'\n'+y,
                          [analyzeTerm.prefix+" "+x for x in tmp.splitlines()])
            ret += '\n' + "-"*60 + '\n'
            pass
        pass

    
    if outFilename:
        import sys
        sys.stdout.close()
        sys.stdout = old_stdout
        pass
    
    # first print overall energies, energies of each term

    #go through each potential type for analysis:
    # print out results if verbose flag is set.
    # bond
    # angle
    # impropers
    #

    # rms, # violations

    # NOE

    # rdc:
    #
    # deal with ensemble, non-ensemble cases
    return ret
    
#dictionary keyed by potName each member of which is a list of two-membered
# tuples: (name,function)
# where function is to be called on a Pot term, like this function(term) and
# return a floating value
extraStats={}
def registerExtraStats(potType,name,function,supportsList=False):
    """register extra terms, averages over selected structures will be
    reported by <m restraintStats>. 

    The four arguments are the potential type as given by the potName
    accessor, the name of the property to be averaged (a string), a
    function to be called on the pot term, like this: function(term) and
    which returns a floating value, and whether this function supports a
    list of terms. The function argument can also return a list of float
    values, but in this case registerExtraStats should be called only a
    single time for each potType.
    """
    if potType not in extraStats:
        extraStats[potType] = [(name,function,supportsList)]
    else:
        extraStats[potType].append( (name,function,supportsList) )
        pass
    return
