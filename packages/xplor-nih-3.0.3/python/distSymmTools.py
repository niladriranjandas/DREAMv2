
"""
Tools to aid in setup/analysis of the distance symmetry potential,
<m distSymmPot>.DistSymmPot.
"""

def create_DistSymmPot(name,
                       restraints=None,
                       sim=None,
                       sim2=None):
    """
    Create a DistSymmPot term with distances specified in a file or sequence.
    
    The form of the sequence should be:
    [[ (atom1, atom2), (atom3, atom4) , ...], ...]
    where atomN is a atom selection string which selects a single
    atom, and each pair of atoms specifies a distance. The sequence
    should contain at least two atom pairs.

    sim and sim2, if specified, specify the simulations of the first and
    second atoms, respectively.

    """
    from distSymmPot import DistSymmPot

    if not sim:
        from simulation import currentSimulation
        sim = currentSimulation()
        pass
    if not sim2:
        sim2 = sim
        pass

    pot = DistSymmPot(name,sim,sim2)

    import re
    if type(restraints)==type("string"):
        table=open(restraints).read()
        restraints=[]
        start=0
        match=re.search(r"assi(gn)? *[^(]",table[start:],re.IGNORECASE)
        while match:
            r=[]
            nextMatch=re.search(r"assi(gn)? *[^(]",table[start+4:],
                                re.IGNORECASE)
            
            subString=table[start:nextMatch.end() if nextMatch else -1 ]
            matches= re.findall(r"\([^)]*\)",subString)
            if not matches or len(matches)>4 or len(matches)%2==1:
                raise Exception("bad number of atom selections in: " +
                                subString)
            cnt=0
            while cnt<len(matches):
                (atom1,atom2) = (matches[cnt][1:-1],
                                 matches[cnt+1][1:-1])
                cnt += 2
                r.append( (atom1,atom2) )
                pass
            
            start+=len(subString)
            match=nextMatch
            
            restraints.append(r)
            pass
        pass

    if restraints!=None:
        for r in restraints:
            if len(r)<2 or len(r)%2==1:
                raise Exception("bad number of atom selections: " + str(r))
            pot.addRestraint(r)
            pass
        pass
    
    return pot

def genDimerRestraints(resids,
                       resids2=None,
                       segids=["A","B"],
                       extraSels=None,
                       atomName="CA"):
    """
    generate a set of restraints for DistSymmPot appropriate for a dimer
    with identical subunits in the specified segids, with the sequence
    of residues numbers in the resids argument corresponding obeying
    the symmetry relationship with the corresponding residue numbers
    in resids2. Restraints will be created, one for each specified residue
    using the specified atom name. 

    If resids2=None it is assumed that the two 
    subunits have the same residue numbering but different segids.

    If segids=None, no segid is used in the atom selection.

    The extraSels argument, if specified, should have two elements to be
    used in place of the segids argument in generating pairs of selection
    selection strings specifying unique atoms.

    If resids2 is specified len(resids2) should be equal to len(resids).
    In this case segids are omitted.

    """
    ret=[]

    if not resids2: 
        resids2=resids
    else:
        segids=[]

    if extraSels:
        extraSels = [ "("+ e + ") and " for e in extraSels ]
    else:
        if segids:
            extraSels = ( 'segid "%s" and ' % segids[0],
                          'segid "%s" and ' % segids[1] )
        else:
            extraSels = ('','')
            pass
        pass

    if len(resids) != len(resids2):
        raise Exception("len(resids) != len(resids2)")
    for i in range(len(resids)):
        t1=t2=t3=t4=None
        j = len(resids)-1-i
        t1 = ('%s resid %d and name %s'% (extraSels[0],resids[i], atomName) )
        t2 = ('%s resid %d and name %s'% (extraSels[1],resids2[j],atomName) )
        t3 = ('%s resid %d and name %s'% (extraSels[1],resids2[i],atomName) )
        t4 = ('%s resid %d and name %s'% (extraSels[0],resids[j], atomName) ) 
      
        ret.append( [(t1,t2),(t3,t4)] )
        pass
    return ret
        
    

def genPolyRestraints(resids,
                      resids2=None,
                      resids3=None,
                      segids=["A","B","C"],
                      atomName="CA"):
    """
    generate a set of restraints for DistSymmPot appropriate for a polymer
    with identical subunits in the specified segids, with the sequence
    of residues numbers in the resids argument corresponding obeying
    the symmetry relationship with the corresponding residue numbers
    in resids2. Restraints will be created, one for each specified residue
    using the specified atom name.

    These restraints will enforce the condition that the rotation+translation
    taking subunit 1 to subunit 2 is the same as that taking subunit 2 to
    subunit 3.

    If resids2=None it is assumed that the two 
    subunits have the same residue numbering but different segids.

    If segids=None, no segid is used in the atom selection.

    If resids2 and resids3 are specified len(resids2) and len(resids3) should
    be equal to len(resids). In this case segids are omitted.

    """
    ret=[]

    if not resids2: 
        resids2=resids
        resids3=resids
    else:
        segids=[]

    if (len(resids) != len(resids2)) or (len(resids2) != len(resids3)):
        raise Exception("len(resids) != len(resids2) or " +
                        "len(resids2) != len(resids3)")

    if len(segids)<2:
        return ret

    for i in range(len(resids)):
        t1=t2=t3=t4=None
        j = len(resids)-1-i
        if segids:
            t1 = ( 'segid "%s" and resid %d and name %s' % (segids[0],
                                                            resids[i],
                                                            atomName) )
            t2 = ( 'segid "%s" and resid %d and name %s' % (segids[1],
                                                            resids2[j],
                                                            atomName) )
            t3 = ( 'segid "%s" and resid %d and name %s' % (segids[1],
                                                            resids2[i],
                                                            atomName) )
            t4 = ( 'segid "%s" and resid %d and name %s' % (segids[2],
                                                            resids[j],
                                                            atomName) )
        else:
             t1 = ( 'resid %d and name %s' % ( resids[i], atomName) )
             t2 = ( 'resid %d and name %s' % ( resids2[j], atomName) )
             t3 = ( 'resid %d and name %s' % ( resids2[i], atomName) )
             t4 = ( 'resid %d and name %s' % ( resids3[j], atomName) ) 
      
        ret.append( [(t1,t2),(t3,t4)] )
        pass
    return ret
        
    

def analyze(potList):
    """perform analysis of DistSymmPot terms and return nicely formatted
    summary"""

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'DistSymmPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()


    ret += "%-9s  %6s  %6s\n" % \
           ( "", "RMS", "Viols")

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]


        print(term.showViolations())
        print(term.info())

        ret += "%-9s  %6.3f  %6d\n" % \
               (name , term.rms(), term.violations() )
        pass
    
    return ret


from simulationTools import registerTerm
registerTerm(analyze,"Distance Symmetry Restraints","DistSymm",
r"""
For each term the root mean square fit of calculated to experiment is printed,
along with the number of violations. 
""")
    
