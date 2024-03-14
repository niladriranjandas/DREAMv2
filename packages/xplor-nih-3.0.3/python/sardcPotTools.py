
"""tools to aid in setup/analysis of dipolar coupling potential term sardcPot

this module provides functions to simplify the creation, manipulation and
analysis of <m sardcPot>.SARDCPot potential terms.
"""

def create_SARDCPot(name,
                    file=0,
                    useSign=True,
                    useDist=False,
                    restraints="",
                    domainSel="not PSEUDO",
                    tensor=None
                    ):
    """create an <m sardcPot>.SARDCPot with given name, the filename of an rdc
    assignment table and/or a string of assignments.

    The file argument can optionally be a sequence of filenames.

    domainSel is an <m atomSel>.AtomSel which specifies atoms in an aligning
    subunit.

    If tensor is specified, it should be a <saTensor>.SATensor object.

    One can specify the useSign argument to specify whether the RDC sign will be
    significant in evaluating the energy. The useDist argument specifies whether
    the 1/r^3 RDC dependence will be taken from the coordinates.
    
    """

    from simulation import currentSimulation
    from saTensor import SATensor
    from sardcPot import SARDCPot

    if not tensor:
        tensor = SATensor("saTensor",domainSel)
        pass

    if file:
        files=[]
        if type(file)==type("string"):
            files.append( file)
        else:
            files=file
            pass
        for file in files:
            restraints += open(file).read()
            pass
        pass


    rdc = SARDCPot(name,tensor,"")

    addRestraints(rdc,restraints,useSign,useDist)


    if rdc.numRestraints()==0 and (file or restraints):
        print("create_SARDCPot: Warning: no restraints read.")
        pass

    from atomSel import AtomSel

    return rdc

def addRestraints(rdc,
                  restraints,
                  useSign=True,
                  useDist=False):
    """
    Add XPLOR XDIPO or SARDC style restraints from a string.

    One can specify the useSign argument to specify whether the RDC sign will be
    significant in evaluating the energy. The useDist argument specifies whether
    the 1/r^3 RDC dependence will be taken from the coordinates.
    """

    import potUtils
    restraints = potUtils.stripBracketedComments( restraints )

    import re
    while restraints:
        l_useDist = useDist
        l_useSign = useSign
        restraints = restraints.lstrip()
        while re.match("!.*",restraints):
            restraints = re.sub(r'![^\n]*','',restraints,1)
            restraints = restraints.lstrip()
            pass
            
        m=re.match("assi[a-z]*[^(]*",restraints,re.IGNORECASE)
        if m:
            (restraint,
             l_useDist,
             l_useSign,
             restraints) = readRestraint(restraints[len(m.group()):],
                                         l_useSign,l_useDist)
            rdc.addRestraints('assign ' + restraint,l_useDist,l_useSign)
            
            pass
        else:
            restraints = re.sub(r'\S*','',restraints,1)
        pass
    return 

def readRestraint(string,
                  useSign,
                  useDist):
    from parseTools import readFloat, readInt, findNested
    
    isXplorRestraint=False
    string = string.strip()
    i = findNested('(',')',0,string,0)
    sel1,string = string[1:i], string[i+1:]
    string = string.strip()
    i = findNested('(',')',0,string,0)
    sel2,string = string[1:i], string[i+1:]
    string = string.strip()
    if string[0]=='(':
        isXplorRestraint=True
        # read dummy sels
        i = findNested('(',')',0,string,0)
        sel2,string = string[1:i], string[i+1:]
        string = string.strip()
        i = findNested('(',')',0,string,0)
        sel2,string = string[1:i], string[i+1:]
        string = string.strip()

        i = findNested('(',')',0,string,0)
        sel1,string = string[1:i], string[i+1:]
        string = string.strip()
        i = findNested('(',')',0,string,0)
        sel2,string = string[1:i], string[i+1:]
        string = string.strip()
    

        pass

    (obs,string) = readFloat(string)
    (err,string) = readFloat(string)

    try:
        (err2,string) = readFloat(string)
    except ValueError:
        err2 = err
        pass

    from math import pi as PI
    
    gammas = {'H' : 26.752e7,
              'C' : 6.7283e7,
              'N' : -2.712e7}
    mu_0     = 4*PI*1.0e-7;
    h_planck = 6.6256e-34;
    h_bar    = h_planck/(2*PI);
    
    NUM_FAC = -h_bar*mu_0/(4*PI*PI*1.0e-30);

    from atomSel import AtomSel
    try:
        atom1=AtomSel(sel1)[0]
    except IndexError:
        print("no atom found for atom1 selection: %s" % sel1)
        raise
    try:
        atom2=AtomSel(sel2)[0]
    except IndexError:
        print("no atom found for atom2 selection: %s" % sel2)
        raise

    gamma1 = gammas[ atom1.atomName()[0] ]
    gamma2 = gammas[ atom2.atomName()[0] ]

    scale = gamma1 * gamma2 * NUM_FAC

    if isXplorRestraint:
        #deduce bond length from sel1, sel2
        idName = atom1.atomName()+"_"+atom2.atomName()
        try:
            dist = distances[idName]
        except KeyError:
            idName = atom2.atomName()+"_"+atom1.atomName()
            try:
                dist = distances[idName]
            except KeyError:
                dist=0.
                useDist=True
                pass
            pass
        pass
    else:
        (dist,string) = readFloat(string)
        (useSign,string) = readInt(string)
        useSign = not useSign #the sardc usesign field is defined backwards
        useDist = True if dist==.0 else False
        pass
    
    if not useDist:
        scale /= dist**3
        pass


    
    return ("(%s) (%s) %f %f %f %f" %(sel1,sel2,obs,err,err2,scale),
            useDist,
            useSign,
            string)
    

distances={ 'N_HN'  : 1.023,    # values taken from Huang, Grzesiek test script
            'CA_C'  : 1.520,
            'CA_CB' : 1.520,
            'C_N'   : 1.330,
            }

def Rfactor(pot,
            selection='all',
            normalize=False):
    """R-factor (in percent) for an infinite number of
    randomly distributed vectors.
    Eq. 3 from Clore+Garrett, JACS 121, 9008 (1999).

      R = 1/sqrt(2) sqrt( \sum (D_i^calc - D_i^obs)^2 / \sum D_i^obs^2 )

    The selection argument can be used to choose a subset of
    restraints whose atoms lie in selection.

    If normalize is set to True, all values of D_i^calc and D_i^obs are
    divided by <m sardcPot>.SARDCPot_Restraint's Dmax in the above equation. 
     """
    
    from selectTools import convertToAtomSel
    selection=convertToAtomSel(selection,pot.simulation().subSim())
    import atomSel

    diff2=0
    obs2=0
    for r in pot.restraints():
        norm = r.Dmax() if normalize else 1.0
        if atomSel.intersection(r.aSelection(),selection):
            calcd = r.calcd()
            obs = r.obs()
            diff = calcd - obs if r.useSign() else abs(calcd) - abs(obs)
            diff2 += (diff / norm)**2
            obs2 +=  (obs/norm)**2
            pass
        pass
    
    from math import sqrt

    if obs2>0.:
        Rinf = 100 * sqrt(diff2) / sqrt(2*obs2)
    else:
        Rinf=-1
        pass
        
    return Rinf

def chi2(pot,
         selection='all'):
    """
    Compute the Chi^2 value for the specified restraints. 

    The optional selection argument can be used to choose a subset of
    restraints whose atoms lie in selection.
    """

    if selection=='all':
        return pot.chisq()

    from selectTools import convertToAtomSel
    selection=convertToAtomSel(selection,pot.simulation().subSim())
    import atomSel

    chi2=0
    num=0
    for r in pot.restraints():
        if atomSel.intersection(r.aSelection(),selection):
            num+=1
            chi2 += ( (r.calcd() - r.obs()) / r.err() )**2
            pass
        pass

    if num>0:
        chi2 /= num
        pass

    return chi2

def avecScale(pot):
    "return pot.rdc.avectorScale()"
    return pot.avectorScale()
    

def rmsd(pot,
         selection='all'):
    """
    Compute the RMSD value for the specified restraints. 

    The optional selection argument can be used to choose a subset of
    restraints whose atoms lie in selection.
    """

    if selection=='all':
        return pot.rms()

    from selectTools import convertToAtomSel
    selection=convertToAtomSel(selection,pot.simulation().subSim())
    import atomSel
    from math import sqrt

    sum=0
    num=0
    for r in pot.restraints():
        if atomSel.intersection(r.aSelection(),selection):
            num+=1
            sum +=  (r.calcd() - r.obs())**2
            pass
        pass

    if num>0:
        sum /= num
        pass

    return sqrt(sum)
    

def makeTable(rdc):
    """Return the assignment table (a string) corresponding to the 
    restraints associated with the specified <m sardcPot>.SARDCPot. 

    """
    ret=""
    for restraint in rdc.rawRestraints():
        ret += "assign "
        ret += "( resid %4d and name OO)\n" % 99999
        ret += "\t( resid %4d and name Z )\n"  % 99999
        ret += "\t( resid %4d and name X )\n"  % 99999
        ret += "\t( resid %4d and name Y )\n"  % 99999
        ret += "\t( %s )\n"  % restraint.aSelection().string()
        ret += "\t( %s )  "  % restraint.bSelection().string()
        ret += "%7.4f %7.4f\n\n" % (restraint.obs(),
                                    restraint.err(),)
        pass
    return ret

def saupeMatrix(pot,
                eIndex=None):
    """given a SARDCPot term, return the associated Saupe matrix. If eIndex
    is specified, return the Saupe matrix associated with the specified
    ensemble member. NOT YET TESTED!
    """
    if eIndex==None:
        eIndex=pot.simulation().member().memberIndex()
        pass

    iTensor = pot.irredTensor(eIndex)
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

def writeNEF(rdc,name):
    """
    Return a formatted NEF record for the restraints in the given SARDCPot
    object as a string with the specified saveframe name.
    """

    from cif import Cif, CifDatablock, CifCategory
    block = CifDatablock()

    cat = CifCategory()
    cat["sf_category"] = "nef_rdc_restraint_list"
    cat["sf_framecode"] = name
    cat["tensor_magnitude"] = '.'#str(11.0000)
    cat["tensor_rhombicity"] = '.'#str(0.0670)

    block["nef_rdc_restraint_list"] = cat

    cat = CifCategory()
    for key in (
      "restraint_id",
#      "restraint_combination_id",
      "chain_code_1",
      "sequence_code_1",
      "residue_name_1",
      "atom_name_1",
      "chain_code_2",
      "sequence_code_2",
      "residue_name_2",
      "atom_name_2",
      "weight",
      "target_value",
      "use_sign",
      "lower_limit",
      "upper_limit",
      "scale",
      ):
        cat.addKey(key)
        pass
    
    from nefTools import toNefAtomname
    id=0
    for r in rdc.restraints():
        id += 1
        obs=r.obs()
        lower = obs - r.minusErr()
        upper = obs + r.plusErr()

        for a in r.aSelection():
            for b in r.bSelection():
                cat.addValue("restraint_id",str(id))
                cat.addValue("chain_code_1",a.segmentName())
                cat.addValue("sequence_code_1",str(a.residueNum()))
                cat.addValue("residue_name_1",a.residueName())
                cat.addValue("atom_name_1",toNefAtomname(a.atomName()))
                cat.addValue("chain_code_2",b.segmentName())
                cat.addValue("sequence_code_2",str(b.residueNum()))
                cat.addValue("residue_name_2",b.residueName())
                cat.addValue("atom_name_2",toNefAtomname(b.atomName()))
                cat.addValue("target_value",str(obs))
                cat.addValue("lower_limit",str(lower))
                cat.addValue("upper_limit",str(upper))
                cat.addValue("use_sign","True" if r.useSign() else "False")
                cat.addValue("weight","1.0")
                cat.addValue("scale","1.0")
                pass
            pass
        pass
    block["nef_rdc_restraint"] = cat

    cif=Cif()
    block.setIsSaveframe(True)
    cif[name] = block
    cif.setUseTrailingStop(True)
    cif.setUseTrailingSave(True)
    return cif.asString()


def analyze(potList):
    "perform analysis of SARDCPot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'SARDCPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret+= "%9s  %6s  %6s  %8s %5s  %4s\n" % \
          (" " , "RMS", "Chi^2", "R-factor" , "avector scale", "Viols")

    rdcs = []
    for name in instanceNames:
        rdc = [x for x in potList if x.instanceName()==name][0]
        rdcs.append(rdc)

        print(rdc.showViolations())

        print(rdc.info())

        ret += "%-9s  %6.3f  %6.3f %8.2f      %5.3f     %4d\n" % \
               (name , rdc.rms(), rdc.chisq(), Rfactor(rdc),
                rdc.avectorScale(),
                rdc.violations() )
        pass
    
    return ret

import saTensorTools

from simulationTools import registerTerm, registerExtraStats
registerTerm(analyze,"Steric Alignment Dipolar Coupling Analysis","SARDC",
r"""
For each term, report 
  RMS           - root mean square deviation between calculated and target
                  RDC values.
  Chi^2         - :math:`\chi^2` value - note that RDC error is significant in
                  computing this.
  R-factor      - quality metric- Eq. 3 from Clore+Garrett, JACS 121, 9008
                  (1999).
  avector scale - the overall scale factor (determined from fit to experiment)
                  applied to the alignment tensor.

  Viols         - number of violated restraints.

For more information, see <m sardcPot>.
""")

registerExtraStats("SARDCPot","chi^2",chi2)
registerExtraStats("SARDCPot","R-factor",Rfactor)
registerExtraStats("SARDCPot","avec scale",avecScale)
