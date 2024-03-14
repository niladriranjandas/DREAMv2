def create_GyrPot(name,sel="not PSEUDO",
                  targetVol=-1,
                  minResScaleFactor=12.28,
                  maxResScaleFactor=18.66,
                  softScale=1e-4,
                  hardScale=1
                  ):
    """
    create a <m gyrPot>.GyrPot with the appropriate target volume for the
    given atom selection

    The relationship between number of residues and gyration ellipsoid volume
    has been found to lie in the range minVol..maxVol, where

            minVol = minResScaleFactor * numResidues
            maxVol = maxResScaleFactor * numResidues

    where numResidues is returned by selectTools.numResidues(sel).

    This function creates two potential terms, soft and hard (whose
    names have _s and _h suffixes, repsectively). The first has a target
    value in the middle of minVol and maxVol. The second term has zero
    energy in the whole range minVol to maxVol, but it has a much larger force
    constant. 
    
    Alternately, if a positive value of targetVol is specified, a single
    energy term is included.
    """

    from selectTools import convertToAtomSel
    sel = convertToAtomSel(sel)

    from potList import PotList
    pl = PotList(name)

    from gyrPot import GyrPot
    hpot = GyrPot(name+"_h",sel)
    hpot.setScale(hardScale)

    import selectTools
    numResidues=selectTools.numResidues( sel )
    if targetVol<0:
        minVol = minResScaleFactor * numResidues
        maxVol = maxResScaleFactor * numResidues
        target = 0.5 * (maxVol + minVol)
        range =  0.5 * (maxVol - minVol)
        hpot.setVolTarget( target )
        hpot.setVolRange( range )

        spot = GyrPot(name+"_s",sel)
        spot.setScale(softScale)
        spot.setVolTarget( target )
        pl.append(spot)
        
    else:
        hpot.setVolTarget( targetVol )
        
        
        
        pass


    pl.append(hpot)

    return pl

def create_RgyrPot(name,
                   sel="not pseudo",
                   targetRg=-1,
                   rangeRg=0,
                   gyrWeights=None
                   ):
    """
    create a <m gyrPot>.GyrPot configured as a radius of gyration restraint
    with the specified target value of Rg. If rangeRg is positive, it is taken
    as a +/- range about targetRg for which the energy is zero.

    If targetRg is not specified it takes the value:

        2.2 * Nres^{0.38} - 1

    where Nres are the number of residues specified in sel.
    
    """

    from selectTools import convertToAtomSel
    sel = convertToAtomSel(sel)

    from gyrPot import GyrPot
    pot =  GyrPot(name,sel)
    pot.setTargetType('radius')
    if rangeRg>0: pot.setRRange(rangeRg)

    import selectTools
    numResidues=selectTools.numResidues( sel )
    if targetRg<0: targetRg = 2.2 * numResidues**0.38 - 1

    pot.setRTarget( targetRg )
        
    if gyrWeights:
        gyr.setGyrWeights( gyrWeights )
        pass

    return pot

def analyze(potList):
    """perform analysis of GyrPot terms and return nicely formatted
    summary"""

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'GyrPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()


    ret += "%-9s  %10s     %10s   %10s   %6s\n" % \
           ( "", "volume", "target vol", "delta", "Rg")

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]
#
#        print term.showViolations()
        print(term.info())
        #
        ret += "%-9s  %10.3f     %10.3f   %10.3f   %7.3f\n" % \
               (name , term.volume(), term.volTarget(),
                term.volume()- term.volTarget(),
                term.Rg()  )
        pass
    
    return ret


from simulationTools import registerTerm
registerTerm(analyze,"Gyration Tensor Volume terms","GyrPot",
r"""
For each term, print the computed and observed gyration volume, and the
difference. Also print the computed radius of gyration. The gyration volume is
defined as

.. math::

  V_{gyr} =   4/3\pi \sqrt{ |G| }

with

.. math::

  G = \frac{1}{N} \sum_{i=1}^N \Delta q_i\otimes \Delta q_i,

where :math:`\Delta q_i\otimes \Delta q_i` represents an outer product of the
position corresponding to atom :math:`i`.

Please see
    C.D. Schwieters and G.M. Clore, "A pseudopotential for improving
    the packing of ellipsoidal protein structures determined by
    NMR," J. Phys. Chem. B 112, 6070-6073 (2008).
""")
    
def Rgyr(potTerm): return potTerm.Rg()
def Vgyr(potTerm): return potTerm.volume()

from simulationTools import registerExtraStats
registerExtraStats("GyrPot","Rgyr",Rgyr)
registerExtraStats("GyrPot","Vgyr",Vgyr)
