
"""tools to create and operate the potential term restraining the shape of
protein and bond orientations using NMR relaxation data.
"""

import protocol

import	time, math, random
from	cdsMatrix		import *
from 	cdsVector 		import CDSVector_double as vector
from 	cdsVector		import norm
from	math			import sqrt, pi, fabs
from	atomSel		import AtomSel, intersection
from	mat3			import	Mat3, transpose
from	vec3			import	cross, dot
from	simulation		import currentSimulation
from 	relaxRatioPot 	import RelaxRatioPot
from	diffPotTools		import make_ratio
from	surfTessellation 	import SurfTessellation
from	selectTools		import convertToAtomSel, makeAtomSelStringFromResidList, breakAtomSelString 

protocol.addPseudoResName("TRRP")

def create_RelaxRatioPot(name,
                         data_in=None,
                         restraints="",
                         freq=None,
                         sel = "known and (not PSEUDO)", 
                         inc_sel = "known",
                         temperature = 293,
                         optTemp=False,
                         addAtoms=False,
                         mass=1000,
                         bond_type='NH',
                         link_to=None,
                         sigmaFactor=10000000,
                         CSA=-160 ):

    """Function which creates a <m relaxRatioPot>.RelaxRatioPot object

       Arguments:
       name        - an instanceName of a RelaxRatioPot object.
       data_in     - is the set relaxation data. Use the functions in
                     <m diffPotTools> for reading and processing the data.
       freq        - spectrometer frequency in MHz.
       sel         - an <m atomSel>.AtomSel specifying atoms to define
                     the molecular surface. This argument is not used
                     the link_to argument is specified.
       inc_sel     - an <m atomSel>.AtomSel specifying the relaxation
                     data to use for energy calculations, by atom.
       temperature - the nominal experimental temperature in Kelvins.
       optTemp     - whether or not to perform temperature optimization. 
                     [default: False]
       addAtoms    - synonym for optTemp - deprecated.
       mass        - mass of the pseudo atoms used for temperature
                     optimization. 
       bond_type   - the nuclei involved in relaxation. Currently only
                     'NH' is supported. 
       link_to     - specify that the surface tessellation should be
                     shared from the specified RelaxRatioPot or a standalone
                     <m surfTessellation>SurfTessellation object.
       sigmaFactor - value for the sigmaFactor of the RelaxRatioPot object.
       CSA         - value of Chemical Shielding tensor Anisotropy, in ppm

      Returns  a <m relaxRatioPot>.RelaxRatioPot object

    """
    sel=convertToAtomSel(sel)

    radii=[]
    for atom in sel:
        firstChar = atom.atomName()[0]
        radii.append( radiusMap[firstChar] )
        pass

    import os
    if restraints and os.path.exists(restraints):
        restraints = open(restraints).read()
        pass

    # checking input and making ratio of relaxation rates if necessary #
    if data_in==None and restraints=="":
        raise Exception("create_RelaxRatioPot: No relaxation data given: exiting")
    if freq==None:
        raise Exception("create_RelaxRatioPot: Please specify spectrometer "
                        "frequency: exiting")
    from atomSel import AtomSel, intersection
    inc_sel=convertToAtomSel(inc_sel)


    k=0

    # if link_to is not given explicitly then create it
    # if is given then use it and ignore sel
    if link_to == None:
        from surfTessellation import SurfTessellation
        tess = SurfTessellation(sel,radii)
    else:        
        try:
            tess = link_to.Tessellation() # treat as potential
        except:
            tess=link_to # treat as an rc pointer to tessellation object
            pass
        pass

    # create a potential instance #
    rR_pot = RelaxRatioPot(name,tess)

    rR_pot.setSigmaFactor(sigmaFactor)

    readRestraints(rR_pot,restraints,CSA)
    if data_in:
        for item in data_in:
            from diffPotTools import make_ratio
            make_ratio(item)
            pass
        for item in data_in:
            if item.Ratio != None:
                try:
                    aSel = intersection(item.Sel_N, inc_sel)
                    bSel = intersection(item.Sel_HN, inc_sel)
                    if len(aSel) and len(bSel):
                        rR_pot.addRestraint(intersection(item.Sel_N, inc_sel),
                                            intersection(item.Sel_HN,inc_sel),
                                            item.Ratio,item.Ratio_err,CSA)
                        pass
                    pass
                except:
                    #atom not in inc_sel
                    pass
                pass
            pass
        pass      

    rR_pot.setScale(1.0)
    rR_pot.setScaleNHgrad(1.0)
    rR_pot.setScaleSAgrad(1.0)
    rR_pot.setDiffShell(2.8)
    rR_pot.setDiffRrmsd(0.5)				
    rR_pot.setDiffRstep(30)	
    rR_pot.setFreq(freq)
    if bond_type=='NH':
       rR_pot.setGyroR(9.8656)
       #rR_pot.setdipC2(1.2975)	# for r_NH=1.0202 [A]
       rR_pot.setDipC2(1.1429)	# for r_NH=1.0420 [A]

    # calculations of dipolar coupling term                                                      #
    # dipC2=(mu_0*g_n*g_h*h_bar/(8*pi*r_NH^3))^2                                                 #
    # where                                                                                      #
    # magnetic permittivity of vacuum:  mu_0  =  4*pi*1e-7       [H/m]                           #
    # hydrogen gyro magnetic ratio:     g_n   =  2*pi*42.576*1e6 [Hz/T]                          #
    # 15N gyro magnetic ratio:          g_h   = -2*pi*4.3156*1e6 [Hz/T]                          #
    # Planck constant over 2*pi         h_bar =  1.05457*1e-34   [J s]                           #
    # effective nuclei pair distance    r_NH  =  1.0202          [A]   from neutron diffraction  #
    #                                   r_NH  =  1.0420          [A]   accounting for librations #
    #                                                                                            #
    # r_NH  =  1.0420 is used as the dfault value for calulations                 	             #
    # dipolar coupling used in calculations is in inverse nanoseconds                            #

    rR_pot.setMedianTmp(temperature)

    registerTerm(rR_pot)

    rR_pot.setDiffTmpF(rR_pot.get_TmpF( temperature ))
    print('\n'," ", rR_pot.instanceName(), ": potential restraining",bond_type, "relaxation rates ratio")
    print("  Assuming tumbling in water at", temperature, "K", "with hydration layer",rR_pot.diffShell(),"[A]")
    print("  Accepted",rR_pot.numRestraints(), end=' ')
    print("data points recorded at", rR_pot.freq(),"MHz")
    print("  Initial Energy=%5.2f"%rR_pot.calcEnergy(),"for initial Scale=%0.2f"%rR_pot.scale(),'\n')

    if addAtoms: optTemp = addAtoms
    if optTemp==True: 
        addTmPAtoms(pot=rR_pot, mass=mass)
        protocol.updatePseudoAtoms() 
        pass

    return rR_pot

def readRestraints(term,
                   restraints,
                   csa       ):
    """
    Load Xplor-NIH format restraint table into a <m relaxRatioPot>.RelaxRatioPot
    term.
    """

    import potUtils
    restraints = potUtils.stripBracketedComments( restraints )

    import re

    ret=""

    from selectTools import toAtomSelString
    while restraints:
        restraints = restraints.lstrip()
        while re.match("!.*",restraints):
            restraints = re.sub(r'![^\n]*','',restraints,1)
            restraints = restraints.lstrip()
            pass
            
        m=re.match("assi[a-z]*[^(]*",restraints,re.IGNORECASE)
        if m:
            restraints = readOneRestraint(term,restraints[len(m.group()):],csa)
            pass
        else:
            restraints = re.sub(r'\S*','',restraints,1)
            pass
        pass
    return ret


def readOneRestraint(term,string,csa):
    """
    read one Xplor-NIH style restraint from the string argument.
    """
    from parseTools import readFloat, readInt, findNested
    
    isXplorRestraint=False
    string = string.strip()
    i = findNested('(',')',0,string,0)
    sel1,string = string[1:i], string[i+1:]
    string = string.strip()
    i = findNested('(',')',0,string,0)
    sel2,string = string[1:i], string[i+1:]
    string = string.strip()

    (obs,string) = readFloat(string)
    (err,string) = readFloat(string)

    try:
        (err2,string) = readFloat(string)
    except ValueError:
        err2 = err
        pass


    term.addRestraint(sel1,sel2,obs,err,csa)

    return string


# add pseudo atoms for temperature optimization
default_segid=""
default_resid=1700
def addTmPAtoms(tpOatom=None,tpXatom=None,tpYatom=None,
                resid=None,segid=None,pot=None,mass=1000):
    """
    create psf and initial coordinates for pseudo atoms used for optimizing
    diffusion tensor temperature
    """
    print("  CREATING Pseudo Atoms for temperature optimization" '\n')

    if pot==None:
        print("  NOT CREATED: Please provide potential name",'\n')
        return 


    global default_segid, default_resid
    if not segid:
        segid = default_segid
        
    from xplorSimulation import getXplorSimulation
    xSim = getXplorSimulation()
    import protocol
    protocol.initParams("axis")

    from atomSel import AtomSel
    if not resid:
        while len( AtomSel('segid "%s" and resid %d' %
                           (segid,default_resid)))>0:
            default_resid +=1 
            pass
        resid = default_resid
        pass

    if len( AtomSel('segid "%s" and resid %s' % (segid,resid))) > 0:
        raise Exception('segid "%s" and resid %d already exist.' %
                        (segid,resid))

    cmd = psfTmPTemplate.replace('_n__','%-4d'%resid)
    cmd = cmd.replace('SGMT','%-4s'%segid)
    xSim.fastCommand(cmd)
    xSim.syncFrom()


    from simulation import currentSimulation
    currentSimulation().sync()


    from vec3 import unitVec, Vec3

    from atomAction import centerOfMass
    pCM = centerOfMass("known and not resname ANI")
    pO = pCM + Vec3(20,20,20)
    pX = pO + Vec3(1,0,0)
    pY = pO + Vec3(0,1,0)

    resSel = 'segid "%s" and resid %d' % (segid,resid)

    for (oname,xname,yname) in (("TPO","TPX","TPY"),):
        oAtom = AtomSel(resSel + " and name "+oname )[0]
        xAtom = AtomSel(resSel + " and name "+xname )[0]
        yAtom = AtomSel(resSel + " and name "+yname )[0]
        oAtom.setPos( pO )
        xAtom.setPos( pX )
        yAtom.setPos( pY )

    oAtom.setMass(mass) #
    xAtom.setMass(mass) #
    yAtom.setMass(mass) #	

    currentSimulation().sync()

    pot.setTmpAtoms(oAtom,xAtom,yAtom)

    print('\n' "  CREATED for",pot.instanceName(), "with mass=",mass,'\n')

    return (oAtom,xAtom,yAtom)


# the string _n__ is replaced by the residue number
# the string SGMT is replaced by the segment name
psfTmPTemplate = """
structure
PSF

       1 !NTITLE
 REMARKS   prePotTools.py: auto-generated structure parameters

       3 !NATOM
       1 SGMT _n__ TRRP TPO  OOO    0.000000E+00   10.0000           0
       2 SGMT _n__ TRRP TPX  XXX    0.000000E+00   10.0000           0
       3 SGMT _n__ TRRP TPY  YYY    0.000000E+00   10.0000           0

       2 !NBOND: bonds
       1       2       1       3

       0 !NTHETA: angles
       

       0 !NPHI: dihedrals


       0 !NIMPHI: impropers


       0 !NDON: donors


       0 !NACC: acceptors


       0 !NNB

       0       0       0       0

       1       0 !NGRP
       0       0       0

end
"""   

def analyze(potList):
    """
    Perform analysis of RelaxRatioPot terms and return nicely formatted
    summary in the structure file.
    """

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'RelaxRatioPot')

    if not potList: return ret

    #set to get rid of duplicates. list for ordering.
    instanceNames = list(set([x.instanceName() for x in potList]))
    instanceNames.sort()

    # temperature values
    ret+= "%52s\n"% ("\n Temperature for diffusion tensor calculations, [K]")
    ret+= "%8s  %9s  %11s  \n" % \
          (" " , "Optimized","Nominal")

    for name in instanceNames:
        dft = [x for x in potList if x.instanceName()==name][0]

        #dft.forceTessellation()  #- this doesn't work. Why?
        #dft.calcEnergy()
        
        ret += "%-9s  %7.3f    %7.3f%3s%3.1f\n" % \
               (name , dft.curr_Tmp(), dft.medianTmp(),"+/-",
                dft.rangeTmpFit())
        pass

    # rms values
    ret+= "%45s\n" % ("                     Norm.   N     N    Sigma")
    ret+= "%-9s %6s %8s %6s %5s  %6s\n" % ("name","Ch^2","Ch^2","exc.",
                                           "total","CutOff")

    for name in instanceNames:

        dft = [x for x in potList if x.instanceName()==name][0]

        sar = dft.showAllRestraints()
        print()
    
        if sar:
            print('Violated', end=' ')
            pass
        print('RelaxRatioPot restraints in potential term', end=' ')
        print(dft.instanceName(), "( threshold=%f )\n" %dft.threshold())

        if sar:
            print("V", end=' ')
            pass
        print("id %20s   %20s" % ("   atom A  ", "    atom B  "), end=' ')
        print("  %7s %7s %7s %7s %7s %7s %7s %7s" % (
            "  obs","   calcd","   diff","   err","  csa","   R1","   R2",
            "   noe"))
        print("-"*101)

        for (i,r) in enumerate(dft.restraints()):
           if r.violated() or sar:
              if sar:
                 if r.violated():
                    print("*", end=' ')
                 else:
                    print(" ", end=' ')
                    pass
                 pass

              print("%4d (%20s) (%20s)" % (i,
                                           r.aSelection()[0].string(),
                                           r.bSelection()[0].string()), end=' ')
              print("%7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f %7.2f" % (
                 r.obs(), r.calcd(), r.diff(), r.err(), r.csa(),
                 r.r1(0), r.r2(0), r.noe(0)))
              for j in range(1,len(r.aSelection())):
                 if sar:
                    print(" ", end=' ')
                    print("%4s (%20s) (%20s)" % ("",
                                                 r.aSelection()[j].string(),
                                                 r.bSelection()[j].string()), end=' ')
                    print("%7s %7s %7s %7s %7s %7.2f %7.2f %7.2f" % (
                       "","","","","",
                       r.r1(j), r.r2(j), r.noe(j)))
                 pass
              pass
           pass
        
                    
                

        ex_res=[r.aSelection()[0].residueNum()\
                for r in dft.restraints() if r.excluded()]

        print()
        print("For sigma cutoff (%7.2f) the excluded residues are:" % \
              dft.sigmaFactor(), end=' ')
        print(breakAtomSelString(makeAtomSelStringFromResidList(ex_res),
                                 lineLen=70,
                                 spacer="                        ",
                                 initialShift=24))

        in_res=[r for r in dft.restraints() if not r.excluded()]

        chi2=0
        for r in in_res:
            chi2 += (r.obs() - r.calcd())**2 / r.err()**2;
            pass

        if dft.numRestraints()!=len(ex_res):
            chi2_ex = chi2 /(dft.numRestraints()-len(ex_res))
        else:
            chi2_ex=-1
            pass

        ret += "%-9s %7.2f %7.2f %5i %5i %6.2f\n" % \
               (name, chi2, chi2_ex,
                len(ex_res), dft.numRestraints(), dft.sigmaFactor())

        print()
        print(dft.info())

        from diffPotTools import diffValuesEuler, tensParams
        (a,b,c,alpha,beta,gamma) = diffValuesEuler(dft.Diff_Tensor())
        (tau,anis,rhomb) = tensParams(a,b,c)
        print("Euler angles: alpha: %8.3f" % alpha)
        print("              beta:  %8.3f" % beta)
        print("              gamma: %8.3f" % gamma)
        print()
        print("  Tau_c:       %8.3f" % tau)
        print("  anis :       %8.3f" % anis)
        print("  rhombicity : %8.3f" % rhomb)
               
        pass



    return ret


def meanTemperature(term):
    """ term can be a single potential term, or a list of terms
    """
    if term.potName()=='RelaxRatioPot':
        q=term.curr_Tmp()
    else:
        # assume we have a list of terms
        from simulationTools import flattenPotList
        prelist = flattenPotList(term)
        
        snu = 0.0; cnt=0
        for p in prelist:
            cnt+=1
            snu += p.curr_Tmp()
            pass

        q = snu / cnt
        pass
    return q

def meanExResN(term):
    """ term can be a single potential term, or a list of terms.

    Return the number of excluded residues.
    """
    if term.potName()=='RelaxRatioPot':

        ret = len([r.aSelection()[0].residueNum()\
                   for r in term.restraints() if r.excluded()])
    else:
        # assume we have a list of terms
        from simulationTools import flattenPotList
        prelist = flattenPotList(term)

        if not prelist:
            return 0;
        
        sum = 0.0
        for p in prelist:
            snu += meanExResN(p)
            pass

        ret = snu / len(prelist)
        pass
    return ret

def makeTable(term):
    """Return the assignment table (a string) corresponding to the 
    restraints associated with the specified <m relaxRatioPot>.RelaxRatioPot. 

    """
    ret=""
    for restraint in term.restraints():
        ret += "assign "
        ret += "\t( %s )\n"  % restraint.aSelection().string()
        ret += "\t( %s )  "  % restraint.bSelection().string()
        ret += "%7.4f %7.4f\n\n" % (restraint.obs(),
                                    restraint.err())
        pass
    return ret
    

def meanChi2ExN(term):
    """ term can be a single potential term, or a list of terms

    Return the normalized chi^2 value.
    """
    if term.potName()=='RelaxRatioPot':
        in_res=[r for r in term.restraints() if not r.excluded()]

        if not in_res:
            return 0.
        chi2=0
        for r in in_res:
            chi2 += (r.obs() - r.calcd())**2 / r.err()**2;
            pass
        ret=chi2/len(in_res)
    else:
        # assume we have a list of terms
        from simulationTools import flattenPotList
        prelist = flattenPotList(term)

        if not prelist:
            return 0.
        
        sum = 0.0
        for p in prelist:
            snu += meanChi2ExN(p)
            pass

        ret = snu / len(prelist)
        pass
    return ret


import simulationTools
simulationTools.registerTerm(analyze,
                             "Relaxation Rates Ratio Potential(s)",
                             "RelaxRatioPot",
r"""
For each term, report
  the optimized diffusion tensor temperature
  the nominal temperature
  normalized chi^2
  unnormalized chi^2
  number of excluded restraints
  total number of restraints
  the sigmaFactor as documented here <m relaxRatioPot>
""")

simulationTools.registerExtraStats("RelaxRatioPot",
                                   "Opt. Temp. [K]",meanTemperature, True)
simulationTools.registerExtraStats("RelaxRatioPot",
                                   "Num. exc. res.",meanExResN, True)
simulationTools.registerExtraStats("RelaxRatioPot",
                                   "Norm. Chi^2",meanChi2ExN, True)

def topologySetup(ivm,list=[]):
    """
    configure the given <m ivm>.IVM object's topology setup using the
    freedom string for each RelaxRatioPot in list.
    This function should be called prior to
    <m ivm>.IVM.autoTorsion() or <m protocol>.torsionTopology()

    """

    try:
        len(list)
    except TypeError:
        list = [list]
        pass

    from simulation import currentSimulation
    global registeredTerms
    if not list:
        try:
            list = registeredTerms[currentSimulation().lookupID()]
        except:
            pass

        pass
    from simulationWorld import world as simWorld
    
    for p in list:
        if p.potName()!='RelaxRatioPot':
            continue
        if simWorld().logLevel()!='none':
            print("\n  Setting topology for %s pseudo atoms\n" % \
                  p.instanceName())
            pass
        #ivm.fix(p.get_O_TmP_atom()) #
        ivm.fix(p.get_X_TmP_atom())
        ivm.group((p.get_O_TmP_atom(),p.get_Y_TmP_atom()))
        ivm.hinge("bend", p.get_O_TmP_atom(),
                  p.get_X_TmP_atom(),
                  p.get_Y_TmP_atom())
        pass
    return

def massSetup(list=[],axisMass=1000):
    """
    appropriately set masses for pseudo atoms used in optimization of
    the apparent diffusion temperature.

    if list is not specified, then pseudoatoms associated with
    all registered RelaxRatioPot objects are configured. 
    """
    try:
        len(list)
    except TypeError:
        list = [list]
        pass

    from simulation import currentSimulation
    global registeredTerms
    if not list:
        try:
            list = list(registeredTerms[currentSimulation().lookupID()].values())
        except:
            pass

    for t in list:
        for a in [t.oAtom,t.xAtom,t.yAtom]:
            a().setMass(axisMass)
            pass
        pass
    return


registeredTerms={}
def registerTerm(term):
    """
    add the given RelaxRatioPot object to a list associated with its Simulation.
    These objects will be automatically processed by topologySetup and
    massSetup.
    """
    global registeredTerms

    #print term.getAtomSel().simulation()
    sim = term.getAtomSel().simulation() 
    id = sim.lookupID()       
    if not id in registeredTerms: registeredTerms[id] = []

    registeredTerms[id].append(term)

    return


def reset_relaxRatioPot_temp(pot=None, temperature=None):

    if pot==None:
        print('\n',"  No potential given.",'\n',"  EXITING",'\n')
        return 

    if pot.rangeTmpFit()==0:
        print('\n',"  Given Potential,",pot.instanceName(),", has zero range for temperature optimization.",'\n',"  EXITING",'\n')
        return 

    if temperature==None:
        print('\n',"  No temperature value given.",'\n',"  EXITING",'\n')
        return 

    if abs(temperature - pot.medianTmp()) > pot.rangeTmpFit():
        print('\n',"  Given temperature,",temperature,", is outside of the range ",pot.getMedianTmp(),"+/-",pot.rangeTmpFit()," for temperature optimization.",'\n',"  EXITING",'\n')
        return 

    cos_a=(temperature - pot.medianTmp())/pot.rangeTmpFit()

    oAtom = pot.get_O_TmP_atom()
    xAtom = pot.get_X_TmP_atom()
    yAtom = pot.get_Y_TmP_atom()

    from	vec3		import	Vec3

    from math import sqrt

    d_y=Vec3(cos_a,sqrt(1-cos_a*cos_a),0)

    yAtom.setPos( oAtom.pos() + d_y )
    
    print('\n',"  For ",pot.instanceName()," temperature set to ",pot.curr_Tmp(),"K",'\n')

    return

# atom radii to use for diffusion calculations
# first character of an atom's name determines its radius
radiusMap={}
radiusMap['H']=1.00
radiusMap['N']=1.75   
radiusMap['C']=1.85  
radiusMap['O']=1.60
radiusMap['S']=2.00 


# the top routine which used to filter out "bad" relaxation data
def filterRelaxData(data_in, Fr,
                    domain_sel="all",
                    fraction=0.13,
                    outFilename=None,
                    ):

    """Function which identifies the largest outliers in the
       input relaxation ratio data for rigid body docking.

       The procedure consists of determining a diffusion tensor which
       best-fits the rho (R2/R1) relaxation data given the HN
       orientations of each domain to be docked. During each iteration
       the data point which has the largest deviation

           (rho^{calcd}-rho^{obs}) / (rho^{calcd}+rho^{obs})

       is discarded. The process of fitting and discarding is repeated
       until the specified fraction of the input data has been
       discarded.

       In this implementation, the diffusion tensor is represented by
       three eigenvalues and three Euler angles. The Euler angles
       represent the orientation of the diffusion tensor relative to
       the input NH bond vectors and are fit, one per domain. One set
       of eigenvalues is set for all domains, with each domain getting
       a separate eigenvalue scale factor, corresponding to allowing
       separate correlation times for each domain.

       Arguments:
         data_in - contains relaxation data, as read by
                   <m diffPotTools>.readInRelaxData and possibly merged by
                   <m diffPotTools>.mergeRelaxData.

         Fr         - a sequence of spectrometer frequencies (in MHz)
                      one for each rigid body domain to be docked.

         domain_sel - a sequence of Xplor-NIH selection strings
                      corresponding to the Fr argument, used to
                      identify atoms in each rigid body domain to be
                      docked.
                      
         fraction   - the fraction of relaxation data to remove.
                      [Default: 13% - the "1.5 sigma rule"]

         outFilename- if specified, the R2/R1 ratios calculated with
                      the final tensor parameters are written out.


       Returns a text string of excluded residues formatted as an
       Xplor-NIH atom selection.

       Limitation: there cannot be any overlap in the residue numbers of
       the rigid domains.
       """ 

    # checking input and making ratio of relaxation rates if necessary #
    import types
    if isinstance(Fr,float) or isinstance(Fr,int):
        Fr=[Fr]
        pass

    if isinstance(domain_sel,bytes):
        domain_sel=[domain_sel]
        pass

    if len(Fr)!=len(domain_sel):
        raise Exception("Please specify as many spectrometer frequencies"
                        "as the number of domain selections.")

    from diffPotTools import make_ratio
    for item in data_in:
       make_ratio(item)
       pass    

    if len(domain_sel)<= 0:
       print('\n', "No domain selection : exiting", '\n')
       return
    else:
       print('\n', len(domain_sel),"Domain selections received")
       pass

    from atomSel import AtomSel, intersection

    # processing initial atom selections and matching them with
    #  relaxation data
    
    atoms_N=[]
    atoms_HN=[]
    in_relax=[]
    in_relax_err=[]
    sel_diff=[]
    #copies which contain all values
    in_relax0 = []
    in_relax_err0 = []
    atoms_N0 = []
    atoms_HN0 = []

    # spectral frequencies and, if necessary, initial guess vector 

    # dipolar coupling term in inverse nanoseconds		
    #d2=1.2975;	# for r_NH=1.0202 [A] 			 
    d2=1.1429		# for r_NH=1.0420 [A]			       
    g=9.8656;	# ratio of giromagnetic ratios for hydrogen and nitrogen  
    omegaLists=[]
    
    k=0
    for domain in domain_sel:
        t_sel=AtomSel(domain)
        t_atoms_N=[]
        t_atoms_HN=[]
        t_in_relax=[]
        t_in_relax_err=[]

        for datum in data_in:  
            tt_atoms_N=intersection(datum.Sel_N,t_sel) 
            tt_atoms_HN=intersection(datum.Sel_HN,t_sel)
            if ((datum.Ratio !=None) and
                (tt_atoms_N.size() > 0) and
                (tt_atoms_HN.size() > 0)):
                t_atoms_N+=tt_atoms_N
                t_atoms_HN+=tt_atoms_HN
                t_in_relax.append(datum.Ratio)
                t_in_relax_err.append(datum.Ratio_err)
                pass
            pass    
        
        if len(t_in_relax) > 0:
            atoms_N.append(t_atoms_N)
            atoms_HN.append(t_atoms_HN)
            in_relax.append(t_in_relax)
            in_relax_err.append(t_in_relax_err)
            atoms_N0.append(list(t_atoms_N))
            atoms_HN0.append(list(t_atoms_HN))
            in_relax0.append(list(t_in_relax))
            in_relax_err0.append(list(t_in_relax_err))
            sel_diff.append(t_sel)
            omegaLists.append( calcOmegaList(Fr[k], d2, g) )
            pass
        k+=1
        pass
    
    numDomains=len(in_relax)
    if numDomains == 0:
        print("All given selections are empty")
        print("or not matching relaxation data: exiting", '\n')
        return

    if numDomains==1:
        print("One non-empty selection:")
    else:
        print(numDomains,"non-empty selections:")
        pass

    all_data_len=0
    for k,itm in enumerate(in_relax):
        print("Domain #%d contains %d relaxation data points" % (k+1,len(itm)))
        all_data_len+=len(itm)
        pass


    # making initial guess tensors for selected domains
    print("Estimating initial guess tensors")
    valuesEulerList=[]
			
    from diffPotTools import create_DiffPot, diffValuesEuler, estimate_Dxyz
    for k,relaxData in enumerate(in_relax):
        dif_s = create_DiffPot('sim_DiffT',sel_diff[k])
        dT_e=dif_s.Diff_Tensor() 
        valuesEuler = diffValuesEuler(dT_e)
        Dxyz=estimate_Dxyz(in_relax[k], omegaLists[k][1])
        ing_dxyz=sum(valuesEuler[:3])/3 # ave of eigenvalues
        scale=Dxyz/ing_dxyz
        for i in range(3): valuesEuler[i] *= scale
        valuesEulerList.append(valuesEuler)
        pass
          
    #parameters to fit include one complete set of eigenvalues, one scale
    # factor for each domain above the first, and Euler angles for
    # each domain
    initialGuess=valuesEulerList[0][:3]

    trase_0=sum(initialGuess)

    for k in range(1,numDomains):
        trase_t=sum(valuesEulerList[k][:3])
        initialGuess.append(trase_t/trase_0)
        pass
    for (a,b,c,gamma,beta,alpha) in valuesEulerList:
        initialGuess += [alpha,beta,gamma]
        pass

    from math import ceil
    it_nmbr=int(ceil(fraction*all_data_len))

    from minimize import simplex
    omittedResids=[]
    for k in range(1,it_nmbr+1):
        print("\nIteration %i out of %i"%(k, it_nmbr))

        simplex_out=simplex(Target_filter_function,initialGuess,
                            args=[[in_relax,in_relax_err,atoms_N,atoms_HN,
                                   omegaLists,d2,g,numDomains]],
                            full_output=True,
                            xtol=1e-9,
                            ftol=1e-9,
                            maxIters=3000000,
                            maxfun=1000000) # simplex fitting routine #

        initialGuess=simplex_out[0]

        print("\t Chi^2 = %.2f" % (simplex_out[1]/(all_data_len-
                                                     len(omittedResids))))
        t_out_res_dat=find_max_dev_res(initialGuess,
                                       in_relax,
                                       in_relax_err,
                                       atoms_N,
                                       atoms_HN,
                                       omegaLists,d2,g,numDomains)

        # remove data point with largest deviation #

        domain_r=t_out_res_dat[1]
        pos_r=t_out_res_dat[2]

        in_relax[domain_r].pop(pos_r)
        in_relax_err[domain_r].pop(pos_r)
        atoms_N[domain_r].pop(pos_r)
        atoms_HN[domain_r].pop(pos_r)
           
        print("\nLargest deviation is for residue #%i"%t_out_res_dat[0], end=' ')
        print("in Domain #%i"%(domain_r+1), '\n')                 

        omittedResids.append(t_out_res_dat[0])

        pass
 
    if outFilename:
        ofile=open(outFilename,"w")
        ofile.write("# resid, R2/R1 calcd, R2/R1 obs, R2/R1 err, calcd-obs\n")

        out_sim=[]
        from diffPotTools import tensParams
        from relaxData import sim_relax_data
        for k in range(numDomains):
            # eigenvalue scale factor
            f_tor = 1 if k==0 else initialGuess[2+k]
            valuesEuler = [ f_tor*v for v in initialGuess[:3] ]
        
            index=2+numDomains+k*3
            valuesEuler += list(initialGuess[index:index+3])

            (alpha,beta,gamma) = valuesEuler[3:6]

            tau,anis,rhomb = tensParams(*valuesEuler[:3])

        
            # returns a list which contains the tuple
            # (resid, R1,R2,NOE, R2/R1)
            t_out_sim=sim_relax_data(atoms_N0[k],atoms_HN0[k],
                                     valuesEuler,
                                     omegaLists[k],
                                     d2,g)
            out_sim.append(t_out_sim)
            pass

        for m,t_out_sim in enumerate(out_sim):
            for k,(resid,R1,R2,noe,rho) in enumerate(t_out_sim):
                status = "excluded" if resid in omittedResids else ""
                ofile.write("%4d %5.2f %5.2f %5.2f %5.2f %s\n" %
                            (atoms_N0[m][k].residueNum(),
                             rho,in_relax0[m][k], in_relax_err0[m][k],
                             rho-in_relax0[m][k],status))
                pass
            pass
        pass

    # convert list of excluded residues into a text string selection in
    #Xplor-NIH format 
    s_o=makeAtomSelStringFromResidList(omittedResids) 

    return s_o

def Target_filter_function(in_param, ifp ):
    """target function for fitting in filtering routine
    """

    (in_relax,in_relax_err,atoms_N,atoms_HN,
     omegaLists,d2,g,numDomains) = ifp

    # simulation of relaxation data: sim_relax_data coded on C++ level
    # (see relaxData.cc)

    out_sim=[]
    
    from relaxData import sim_relax_data
    for k in range(numDomains):
        # eigenvalue scale factor
        f_tor = 1 if k==0 else in_param[2+k]
        valuesEuler = [ f_tor*v for v in in_param[0:3] ]
        
        index=2+numDomains+k*3
        valuesEuler += list(in_param[index:index+3])

        t_out_sim=sim_relax_data(atoms_N[k],atoms_HN[k],
                                 valuesEuler,omegaLists[k],d2,g)
        out_sim.append(t_out_sim)
        pass
     

    cost=0
    #chi^2 * N
    for m,t_out_sim in enumerate(out_sim):
        for k,item in enumerate(t_out_sim):
            cost += ((in_relax[m][k]-item[4])/in_relax_err[m][k])**2
            pass
        pass

#    print 'cost:', cost

    return cost


def find_max_dev_res(in_param,
                     in_relax,
                     in_relax_err,
                     atoms_N,
                     atoms_HN,
                     omegaLists,
                     d2,
                     g,
                     numDomains,
                     ):
    """function which finds a residue with maximum deviation.

    Returns a sequence of [resid, domain_id, data_index].

    """

    max_dev=0
    max_dev_res=0
    max_dom=-1;
    max_pos=-1;

    # simulation of relaxation data: sim_relax_data coded on C++ level
    # (see relaxData.cc)

    out_sim=[]
    
    from relaxData import sim_relax_data
    from diffPotTools import tensParams
    for k in range(numDomains):
        # eigenvalue scale factor
        f_tor = 1 if k==0 else in_param[2+k]
        valuesEuler = [ f_tor*v for v in in_param[:3] ]
        
        index=2+numDomains+k*3
        valuesEuler += list(in_param[index:index+3])

        
        (alpha,beta,gamma) = valuesEuler[3:6]

        tau,anis,rhomb = tensParams(*valuesEuler[:3])

        print('  domain: %d' % (k+1), end=' ')
        print('anisotropy: %.2f  rhombicity: %.2f  tau: %.2f' % (anis,
                                                                 rhomb,tau))
        print('    alpha: %.2f beta: %.2f gamma: %.2f' % (alpha,
                                                          beta,
                                                          gamma))
        
        # returns a list which contains the tuple
        # (resid, R1,R2,NOE, R2/R1)
        t_out_sim=sim_relax_data(atoms_N[k],atoms_HN[k],
                                 valuesEuler,
                                 omegaLists[k],
                                 d2,g)
        out_sim.append(t_out_sim)
        pass


    from math import fabs
    for m,t_out_sim in enumerate(out_sim):
        for k,(resid,R1,R2,noe,rho) in enumerate(t_out_sim):

           tt=fabs(2*(in_relax[m][k]-rho)/(in_relax[m][k]+rho))


           if tt > max_dev:
              max_dev=tt
              max_dev_res=atoms_N[m][k].residueNum()
              max_dom=m
              max_pos=k
              pass
           pass
        pass

    return [max_dev_res,max_dom,max_pos]

def calcOmegaList(Fr, d2, g): #
    """return frequencies at which the spectral density function must
    be evaluated at in computing R2/R1.

      [w_0, w_N, w_H-w_N, w_H, w_H+w_N]
    
     Fr - spectral frequency in MHz
     
     d2 - dipolar coupling term in inverse nanoseconds	(unused)
     g  - ratio of gyromagnetic ratios for hydrogen and nitrogen
     """

    from math import pi
    # omega values are in GHz #
    w_0=0			#omega(0)=w_0
    w_N=Fr*(1e-3)*2*pi/g;	#omega(1)=w_N
    w_H=Fr*(1e-3)*2*pi;	#omega(3)=w_H            
    				#omega(2)=w_H-w_N  
    				#omega(4)=w_H+w_N

    omega=[w_0, w_N, w_H-w_N, w_H, w_H+w_N]

    return omega





