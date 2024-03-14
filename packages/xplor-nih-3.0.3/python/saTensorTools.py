
def analyze(potList):
    """Perform analysis of SATensor terms and return nicely formatted summary.

    This prints information on an alignment tensor whose Da is the
    largest of any SARDCPot terms.
    """

    from simulationTools import getPotTerms
    saTensors = getPotTerms(potList,'SATensor')

    sardc= getPotTerms(potList,'SARDCPot')

    dmax=1
    if len(sardc)>0:
        dmax=-1e30
        for p in sardc:
            #Note this definition of dmax
            # this could be done on a SATensor-to-SATensor basis if
            # registerExptToSATensor were implemented
            dmax=max([dmax] + [r.Dmax() * p.avectorScale() for
                               r in p.restraints()])
            pass
        pass
    
    saTensors += [x.tensor for x in sardc]

    saTensors.sort(key=lambda x: x.instanceName())
    rSaTensors = saTensors
    saTensors=[]
    for t in rSaTensors:
        if not t.instanceName() in [x.instanceName() for x in saTensors]:
            saTensors.append(t)
            pass
        pass

    ret = ""
    for saTensor in saTensors:
        ret += "\nSATensor: %s" % saTensor.instanceName()
        
        from varTensorTools import saupeToVarTensor, VarTensor_analyze
        from varTensorTools import deletePseudoAtoms

        varTens = saupeToVarTensor(saupeMatrix(saTensor),dmax)
        ret += VarTensor_analyze([varTens])
        deletePseudoAtoms(varTens)
        del varTens
        pass
    from simulationTools import gcRegisteredTopoTerms
    gcRegisteredTopoTerms()
    #these get deleted when atoms are added or removed.
    import protocol
    protocol.initDihedrals(reload=True)

    return ret

def saupeMatrix(tensor,
                eIndex=None):
    """given a SATensor, return the associated Saupe matrix. If eIndex
    is specified, return the Saupe matrix associated with the specified
    ensemble member. NOT YET TESTED!
    """
    if eIndex==None:
        eIndex=tensor.esim.member().memberIndex()
        pass

    iTensor = tensor.irredTensor(eIndex)
    from math import sqrt
    Szz = iTensor[0]
    Sxx = 0.5*(-iTensor[0] + sqrt(3) * iTensor[3])
    Syy = 0.5*(-iTensor[0] - sqrt(3) * iTensor[3])
    Sxy = 0.5 * sqrt(3) * iTensor[4]
    Sxz = -0.5 * sqrt(3) * iTensor[1]
    Syz = -0.5 * sqrt(3) * iTensor[2]

    from mat3 import SymMat3

    ret = SymMat3(Sxx,
                  Sxy,Syy,
                  Sxz,Syz,Szz)

    return ret

import simulationTools
simulationTools.registerTerm(analyze,
                             "Steric Alignment Tensor Analysis","SATensor",
r"""
Print <m varTensor>.VarTensor information for the alignment tensor
corresponding to this term's representation.
""")

