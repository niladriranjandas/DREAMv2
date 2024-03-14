"""tools to create a potential term restraining the shape of protein using
   components of rotation diffusion tensor. Together with facilities to 
   fit NMR relaxation data. 
"""

import protocol
protocol.addPseudoResName("TEMP")

from math import sqrt

def create_DiffPot(name,
                   sel = "known and (not PSEUDO)", 
                   temperature = 293, addAtoms=False, mass=1000,
                   link_to=None):

    """ create an instance of an <m diffPot>.DiffPot object
    
        Arguments:																	
        name        - is an instance name of a diffPot object.											
        sel         - is selection of atoms which will be used to create
                      protein surface tessellation by surfTessellation object	
                      This parameter is used only if the created relaxRatioPot
                      instance is not using a surfTessellation 
                      borrowed from another potential, i.e. link_to=None, in
                      which case sel argument is omitted.
        temperature - the nominal experimental temperature in Kelvins 									
        addAtoms    - flag wich specifies whether the pseudo atoms for
                      temperature optimization are created or not. Default
                      setting does not assume temperature optimization. The
                      state addAtoms=True creates pseudo atoms
        mass        - value for the mass of the pseudo atoms										
        link_to     - used to link the potential to surfTessellation created
                      by another potential or to a standing alone
                      surfTessellation object. The default setting that
                      external surfTessellation is not given, and, therefore,
                      creates a new surfTessellation using sel atom selection.

        returns a diffPot object												
    """

    from simulation import currentSimulation
    from diffPot import DiffPot

    # if link_to is not given explicitly then create it
    # if is given then use it and ignore sel
    if link_to == None:
        from selectTools import convertToAtomSel
        sel=convertToAtomSel(sel)
        radii=[]
        for atom in sel:
            firstChar = atom.atomName()[0]
            radii.append( radiusMap[firstChar] )
            pass

        N_atoms = sel.size()

        from surfTessellation import SurfTessellation
        tess = SurfTessellation(sel,radii)

    else:        
        try:
            tess = link_to.Tessellation() # treat as potential
        except:
            tess=link_to # treat as an rc pointer to tessellation object
            pass
        pass

    D_pot = DiffPot(name,tess)
    N_atoms = D_pot.getAtomSel().size()

    # linear scaling calibrated using protein G: equal 200 for a 855 residue protein
    D_pot.setScale_Lin(200*(N_atoms/855)) 

    D_pot.setDiffTmpF(D_pot.get_TmpF( temperature ))
    print("  Assuming tumbling in water at", temperature, "K" '\n')

    D_pot.setMedianTmp(temperature)

    registerTerm(D_pot)

    if addAtoms==True: 
        addTmPAtoms(pot=D_pot, mass=mass)
        protocol.updatePseudoAtoms()  
        pass
    protocol.updatePseudoAtoms()

    return D_pot

# add pseudo atoms for temperature optimization
default_segid=""
default_resid=2700
def addTmPAtoms(tpOatom=None,tpXatom=None,tpYatom=None,
                resid=None,segid=None,pot=None,mass=1000):
    """
    create psf and initial coordinates for pseudo atoms used for optimizing
    correlation times tau_c, tau_t and tau_i
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
    protocol.updatePseudoAtoms(xSim)


    from simulation import currentSimulation
    currentSimulation().sync()


    from vec3 import unitVec, dot, cross, Vec3, norm

    from atomAction import centerOfMass
    pCM = centerOfMass("known and not PSEUDO")
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

    oAtom.setMass(mass)
    xAtom.setMass(mass)
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
       1 SGMT _n__ TEMP TPO  OOO    0.000000E+00   10.0000           0
       2 SGMT _n__ TEMP TPX  XXX    0.000000E+00   10.0000           0
       3 SGMT _n__ TEMP TPY  YYY    0.000000E+00   10.0000           0

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
    "perform analysis of DiffPot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'DiffPot')

    if not potList: return ret

    #set to get rid of duplicates. list for ordering.
    instanceNames = list(set([x.instanceName() for x in potList]))
    instanceNames.sort()

    ret+= "%8s  %9s  %11s  \n" % \
          (" " , "Optimized","Nominal")

    for name in instanceNames:
        dft = [x for x in potList if x.instanceName()==name][0]

        ret += "%-9s  %7.3f    %7.3f%3s%3.1f\n" % \
               (name , dft.curr_Tmp(), dft.getMedianTmp(),"+/-",dft.rangeTmpFit())
        pass
    
    return ret


def meanTemperature(term):
    """ term can be a single potential term, or a list of terms
    """
    if term.potName()=='DiffPot':
        q=term.curr_Tmp()
    else:
        # assume we have a list of terms
        from simulationTools import flattenPotList
        prelist = flattenPotList(term)
        
        from math import sqrt

        snu = 0.0; cnt=0
        for p in prelist:
            cnt+=1
            snu += p.curr_Tmp()
            pass
        
        q = snu / cnt
        pass
    return q


import simulationTools
simulationTools.registerTerm(analyze,
                        "Temperature for diffusion tensor calculations, [K]",
                             "DiffPot",
r"""
For each diffusion tensor, print

The fit temperature, and the allowed temperature range. Note that this
temperature is a fit parameter. Please see

    Y. Ryabov, J.-Y. Suh, A. Grishaev, G.M. Clore and
    C.D. Schwieters, "Using the Experimentally Determined
    Components of the Overall Rotational Diffusion Tensor To
    Restrain Molecular Shape and Size in NMR Structure
    Determination of Globular Proteins and Protein-Protein
    Complexes," J. Am. Chem. Soc. 131, 9522-9531 (2009)
""")

simulationTools.registerExtraStats("DiffPot", "Optimized Temperature, [K]",
                                   meanTemperature, True)

def topologySetup(ivm,list=[]):
    """
    configure the given <m ivm>.IVM object's topology setup using the
    freedom string for each DiffPot in list.
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
        if p.potName()!='DiffPot':
            continue
        if simWorld().logLevel()!='none':
            print("\n  Setting topology for %s pseudo atoms\n" % \
                  p.instanceName())
            pass
        #ivm.fix(p.get_O_TmP_atom())
        ivm.fix(p.get_X_TmP_atom())
        ivm.group((p.get_O_TmP_atom(),p.get_Y_TmP_atom()))
        ivm.hinge("bend", p.get_O_TmP_atom(),
                  p.get_X_TmP_atom(),
                  p.get_Y_TmP_atom())
        pass
    return


registeredTerms={}
def registerTerm(term):
    """
    add the given DiffPot object to a list associated with its Simulation.
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


def reset_DiffPot_temp(pot=None, temperature=None):

    if pot==None:
        print('\n',"  No potential given.",'\n',"  EXITING",'\n')
        return 

    if pot.rangeTmpFit()==0:
        print('\n',"  Given Potential,",pot.instanceName(),", has zero range for temperature optimization.",'\n',"  EXITING",'\n')
        return 

    if temperature==None:
        print('\n',"  No temperature value given.",'\n',"  EXITING",'\n')
        return 

    if abs(temperature - pot.getMedianTmp()) > pot.rangeTmpFit():
        print('\n',"  Given temperature,",temperature,", is outside of the range ",pot.getMedianTmp(),"+/-",pot.rangeTmpFit()," for temperature optimization.",'\n',"  EXITING",'\n')
        return 


    cos_a=(temperature - pot.getMedianTmp())/pot.rangeTmpFit()

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
radiusMap['P']=2.05 # guess CDS 2016/07/06

########################################################################
# the following are the subrutines used for fitting NMR relaxation data
########################################################################

class realx_datum():
    def __init__(self): # initializing constructor
        self.resid=None
        self.Sel_N=None
        self.Sel_HN=None
        self.R1=None
        self.R1_err=None
        self.R2=None
        self.R2_err=None
        self.Ratio=None
        self.Ratio_err=None
        self.NOE=None
        self.NOE_err=None
        pass
    pass

# routine for reading relaxation data files
def readInRelaxData(file,
                    pattern=['resid','ratio','ratio_err'],
                    verbose=False): 
    """ # pattern specifies the way the data given in the file
        # possible keys for pattern are :
        # resid 	- residue number
        # ratio 	- ratio of relaxation rates
        # ratio_err	- errors in ration of relaxation rates
        # R1		- longitudinal relaxation rate
        # R1_err	- error in longitudinal relaxation rate
        # R2		- transverse relaxation rate
        # R2_err	- error in transverse relaxation rate 
        # skip       - skip entry
        # T1		- longitudinal relaxation time
        # T1_err	- error in longitudinal relaxation time
        # T2		- transverse relaxation time
        # T2_err	- error in transverse relaxation time
        # NOE           - heteronuclear NOE
        # NOE_err       - error in heteronuclear NOE

    """

    restraints=""
    relax_data=[]

    if verbose:
        print('reading table: %s...' %file, end=' ')
        pass

    import simulationWorld
    # reading the files line by line
    if type(file)==type("string"):
        files = [ file ]
    else:
        files=file
        pass

    for file in files: 
        restraints = open(file).readlines() # reading file and splitting 
                                            # it into list of lines
    
        num=0
        for line in restraints: # analyze the lines and fill-in relax_data structure
    
            if line.startswith('#'):
                continue
            
            datum_line=str2num(line)
            RD=[]
            if len(datum_line) >= len(pattern):         
                RD=[realx_datum()]	 # create new realx_datum strcuture
                fill_in_relax_datum(RD[0],datum_line,pattern)
                num+=1
            else:
                if verbose:  print("discarding line:", datum_line)
                pass
    
            relax_data +=RD
            pass
        if simulationWorld.world().logLevel() != 'none':
            print('read relaxation data on %d residues from %s' % (num,file))
        pass

    if verbose:
        print('  read %d datapoints' % len(relax_data))

    return relax_data



# routine like above which starts from reading restraints
# for use in test file
def readInRST(restraints=None): 
    from atomSel import AtomSel
    relax_data=[]

    for line in restraints: # analyze the lines and fill-in relax_data structure
       RD=[]         
       RD=[realx_datum()]	
       RD[0].resid=int(line[0])
       RD[0].Sel_N=AtomSel('(resid '+str(RD[0].resid)+' and name N)')
       RD[0].Sel_HN=AtomSel('(resid '+str(RD[0].resid)+' and name HN)')
       RD[0].Ratio=line[1]
       RD[0].Ratio_err=line[2]

       relax_data +=RD

    return relax_data


def str2num(line=None): # extract all the numbers form string
                        # and returns a list of extracted numbers
    return_dat=[]

    if line==None:
        print("str2num: No line given")
        return return_dat
        pass

    for item in line.split():
        LL=[]
        try:
            LL=[float(item)]
        except:
            #not a float
            pass 
        try:
            LL=[int(item)]
        except:
            #not an int
            pass
        return_dat +=LL 
        pass
    return return_dat


# service subroutine which distributes 
# the list of given data into given 
# rerealx_datum structure according given pattern 
def fill_in_relax_datum(datum,datum_line,pattern): 

    from atomSel import AtomSel
    if len(datum_line) < len(pattern): # checking if data match the pattern
       print("fill_in_relax_datum: Pattern line is longer than data provided")
       return
       pass

    position=0

    for item in pattern:
      if item=='skip':
         pass
      if item=='resid':
         datum.resid=int(datum_line[position])
         datum.Sel_N=AtomSel('(resid '+str(datum.resid)+' and name N)')
         datum.Sel_HN=AtomSel('(resid '+str(datum.resid)+' and name HN)')
         pass
      if item=='ratio':
         datum.Ratio=datum_line[position]
         pass
      if item=='ratio_err':
         datum.Ratio_err=datum_line[position]
         pass
      if item=='R1':
         datum.R1=datum_line[position]
         pass
      if item=='R1_err':
         datum.R1_err=datum_line[position]
         pass
      if item=='R2':
         datum.R2=datum_line[position]
         pass
      if item=='R2_err':
         datum.R2_err=datum_line[position]
         pass
      if item=='T1':
         datum.R1=1/(datum_line[position])
         pass
      if item=='T1_err':
         datum.R1_err=(datum_line[position])*datum.R1*datum.R1
         pass
      if item=='T2':
         datum.R2=1/(datum_line[position])
         pass
      if item=='T2_err':
         datum.R2_err=(datum_line[position])*datum.R2*datum.R2
         pass
      if item=='NOE':
         datum.NOE=datum_line[position]
         pass
      if item=='NOE_err':
         datum.NOE_err=datum_line[position]
         pass


      position+=1
      pass

    return

def mergeRelaxData(data_merged):
    """ # subroutine with analyzes give list of realx_datum structures,	
        # sorts the list according to the residue number,			
        # finds items with the same residue numbers,				
        # and merges them in a single relaxa_datum.				
    """
 
    data_merged.sort(key=lambda x: x.resid)

    psn=1;
    resid_c=data_merged[0].resid

    while psn < len(data_merged):
       if data_merged[psn].resid==resid_c:
         R1_m=merge2values(data_merged[psn-1].R1,
         data_merged[psn-1].R1_err,data_merged[psn].R1,data_merged[psn].R1_err)
         R2_m=merge2values(data_merged[psn-1].R2,
         data_merged[psn-1].R2_err,data_merged[psn].R2,data_merged[psn].R2_err)
         Ratio_m=merge2values(data_merged[psn-1].Ratio,
         data_merged[psn-1].Ratio_err,data_merged[psn].Ratio,data_merged[psn].Ratio_err)

         data_merged[psn-1].R1=R1_m.v; data_merged[psn-1].R1_err=R1_m.e;
         data_merged[psn-1].R2=R2_m.v; data_merged[psn-1].R2_err=R2_m.e;
         data_merged[psn-1].Ratio=Ratio_m.v; data_merged[psn-1].Ratio_err=Ratio_m.e;
         
         data_merged=remove_position_from_list(data_merged,psn)
       else:
         resid_c=data_merged[psn].resid
         psn=psn+1 
         pass
       pass

    return data_merged

def merge2values(value1, err1, value2, err2):

    retV=v_e()

    if (value1 != None) & (value2 != None):
      retV.v=(value1+value2)/2
      pass 
    if (value1 != None) & (value2 == None):
      retV.v=value1
      pass 
    if (value1 == None) & (value2 != None):
      retV.v=value2
      pass 

    if (err1 != None) & (err2 != None):
      retV.e=0.5*sqrt(err1*err1+err2*err2)
      pass 
    if (err1 != None) & (err2 == None):
      retV.e=err1
      pass 
    if (err1 == None) & (err2 != None):
      retV.e=err2
      pass 

    return retV

class v_e():
    def __init__(self):
        self.v=None
        self.e=None
        pass
    pass

def show_relax_datum(datum):
    """
    debugging subroutine printing out the content of realx_datum structure
    """
    print("resid:", datum.resid)
    print("Sel_N: ",datum.Sel_N[0].residueName(),datum.Sel_N[0].residueNum(), datum.Sel_N[0].atomName())
    print("Sel_HN:",datum.Sel_HN[0].residueName(),datum.Sel_HN[0].residueNum(), datum.Sel_HN[0].atomName())
    print("ratio:    " , datum.Ratio)
    print("ratio_err:" , datum.Ratio_err)
    print("R1:       " , datum.R1)
    print("R1_err:   " , datum.R1_err)
    print("R2:       " , datum.R2)
    print("R2_err:   " , datum.R2_err,'\n')
    return  

def make_ratio(datum):
    """
    For given realx_datum structure constructs a ratio of relaxation
    rates R2/R1.
    """

    if datum.Ratio != None:
       #print "make_ratio: Relaxation rates ratio is already present for residue #", datum.resid, ": skipping"
       return
       pass

    if (datum.R1 !=None) & (datum.R2 !=None):
       datum.Ratio=datum.R2/datum.R1
       datum.Ratio_err=sqrt((datum.R2_err/datum.R1)**2+(datum.R1_err*datum.R2/(datum.R1**2))**2)
    else:
       print("make_ratio: Either R1 or R2 datum is absent for residue #", datum.resid, ". Not able to calculate the ratio : skipping")
       pass

    return 

# the top routine which fits relaxation data
def fitRelaxData(data_in=None, Fr=None, InG=None, mode='full', inc_sel='known', Err_calc=1):
    """Function which fits NMR NH relaxation data using Nelder-Mead downhill
        (amoeba) simplex method				      
        
       Arguments:
       data_in  - list of realx_datum structures which is generated by
                  readInRelaxData					       
       Fr       - is the spectrometer frequency in MHz													
       mode     - specifies the model for relaxation data fitting											
                    the available options for mode are:												
                    full  - is fully anisotropic tensor												
                    axial - is axially symmetric tensor 												
                    iso   - is spherically symmetric tensor												
       InG      - optional list of initial guess parameters for fitting
                  routine which should be compatiable with mode settings:											
           for full mode  InG=[Dx, Dy, Dz, alpha, betha, gamma]										
           for axial mode InG=[Dxy, Dz, betha, alpha]											
           for iso mode   InG-[Dxyz]														
           Dx, Dy, Dz are in inverse nanoseconds, alpha, betha, and gamma are
             in degrees					     
           if InG is not given then it is estimated for the loaded protein
           coordinates  using facilities in <m diffPot>
														
       inc_sel  - specifies the residues which will which will be included
                  into the fit when relaxation data are available in data_in.
                  Could be either a text string or an atom selection object.		
       Err_calc - number of bootstrap cycles for error calculations 									
           default settings Err_calc=1 which means that only the original data
           is used: no synthetic data sets prepared and no confidence
           intervals calculated. Bootstrap error estimations is lengthy. For
           robust estimation of confidence intervals the number of bootstrap
           cycles, Err_calc, should be equal to the number of R2/R1 relaxation
           data points in the data set used for fitting, len(in_relax).
           However to get an estimate of the order of magnitude for confidence
           intervals one may use about 20% of this number.   

       Returns:																	
           list_results[0], a python structure containing results of fitting
           for the original data set of relaxation data	
           for details see declaration of FitResults() class below.									
    """

    from math import sqrt
    import cdsMatrix
    # checking input and making ratio of relaxation rates if necessary
    if data_in==None:
       print('\n', "fitRelaxData: No relaxation data given: exiting", '\n')
       return
       pass
    if Fr==None:
       print('\n', "fitRelaxData: Please specify spectrometer frequency: exiting", '\n')
       return
       pass
    if (mode != 'full') & (mode != 'axial') & (mode != 'iso'):
       print('\n', "Invalid mode settings.",'\n',"Valid modes settings are: 'full', 'axial', or 'iso'",'\n', "exiting", '\n')
       return
       pass
    for item in data_in:
       make_ratio(item)
       pass    

    # preparing lists of atoms and relaxation rates ratios which will be passed to the relaxation data fitting subroutines
    atoms_N=[]
    atoms_HN=[]
    in_relax=[]
    in_relax_err=[]

    from atomSel import AtomSel, intersection

    try:
      inc_sel=AtomSel(inc_sel) #in this case inc_sel should be a text string
    except:
      #in this case inc_sel should be a selection object
      pass

    k=0
    for item in data_in:
       if item.Ratio != None: 
          atoms_N+=intersection(item.Sel_N,inc_sel)
          atoms_HN+=intersection(item.Sel_HN,inc_sel)
          try:
             tt=intersection(item.Sel_HN,inc_sel)
             tt[0].residueName()
             in_relax.append(item.Ratio)
             in_relax_err.append(item.Ratio_err)
          except:
             #atom not in inc_sel
             pass
          pass      
       pass  

    
    # preparing input for fitting routine: coupling, ratio of gyromagnetic ratios, spectral frequencies and, if necessary, initial guess vector

    # dipolar coupling term in inverse nanoseconds					#
    #d2=1.2975;	# for r_NH=1.0202 [A] 						# 
    d2=1.1429		# for r_NH=1.0420 [A]						#
    g=9.8656;		# ratio of giromagnetic ratios for hydrogen and nitrogen	#

    from math import pi
    # omega values are in GHz
    w_0=0			#omega(0)=w_0
    w_N=Fr*(1e-3)*2*pi/g;	#omega(1)=w_N
    w_H=Fr*(1e-3)*2*pi;	#omega(3)=w_H            
    				#omega(2)=w_H-w_N  
    				#omega(4)=w_H+w_N

    omega=[w_0, w_N, w_H-w_N, w_H, w_H+w_N]

    dT_e=None

    if InG==None:                      # estimation of initial guess parameters if they are not given
       print('\n', "   Initial guess parameters are not given.",'\n',"  Estimating using diffPot facilities.")
       dif_s = create_DiffPot('sim_DiffT')
       dT_e=dif_s.Diff_Tensor()       
       InG=diffValuesEuler(dT_e)
       Dxyz=estimate_Dxyz(in_relax, omega[1])
       print("  Scaling initial guess to Tau=%0.1f [ns] estimated from mean R2/R1 ratio" % (0.5/(3*Dxyz)))
       ing_dxyz=(InG[0]+InG[1]+InG[2])/3
       scale=Dxyz/ing_dxyz
       InG[0]=InG[0]*scale; InG[1]=InG[1]*scale; InG[2]=InG[2]*scale;
       pass

    # change InG to accommodate different fitting modes
    
    if mode=='axial':
       Dxy=(InG[0]+InG[1])/2
       Dz=InG[2]
       al_a=InG[3]
       be_a=InG[4]
       InG=[]; InG.append(Dxy); InG.append(Dz); InG.append(al_a); InG.append(be_a);
       pass
    if mode=='iso':
       Dxyz=(InG[0]+InG[1]+InG[2])/3
       InG=[]; InG.append(Dxyz); 
       pass

         # Fitting data
    k=1;
    import time
    t1=time.time()
    print("Start Fitting")

    in_relaxInG=in_relax; in_relax_errInG=in_relax_err;
    atoms_NInG=atoms_N; atoms_HNInG=atoms_HN;

    in_data_len=len(in_relaxInG)

    # fraction of the original data to be replaced with duplicates for
    # bootstrap error estimations
    fraction_remove=in_data_len*0.37 

    InGInG=InG;

    list_of_Fit_res_list=[]

    while k <=Err_calc: # loop of bootstrap error estimations

         if Err_calc > 1:
            print("Error estimation iteration: %i" % k)
            pass

         from minimize import simplex
         simplex_out=simplex(Target_relax_function,InGInG,
                             [[in_relaxInG,in_relax_errInG,atoms_NInG,atoms_HNInG,omega,d2,g,mode]],
                             full_output=True,xtol=1e-9,ftol=1e-9,maxIters=3000000,maxfun=1000000) # simplex fitting routine


         if Err_calc > 1:
            # generation of a synthetic data set for next iteration

            in_relaxInG=in_relax; in_relax_errInG=in_relax_err; atoms_NInG=atoms_N; atoms_HNInG=atoms_HN;
            in_relaxSrC=in_relax; in_relax_errSrC=in_relax_err; atoms_NSrC=atoms_N; atoms_HNSrC=atoms_HN;
            count = 1;
            add_relax=list(); add_relax_err=list(); add_atoms_N=list(); add_atoms_HN=list();

            import math, random
            while count <= fraction_remove:

              p2r=int(math.floor(len(in_relaxInG)*random.random())) # position to remove
              p4r=int(math.floor(len(in_relaxSrC)*random.random())) # position for replacment

              in_relaxInG=remove_position_from_list(in_relaxInG,p2r)
              in_relax_errInG=remove_position_from_list(in_relax_errInG,p2r)
              atoms_NInG=remove_position_from_list(atoms_NInG,p2r)
              atoms_HNInG=remove_position_from_list(atoms_HNInG,p2r)

              
              if count==1:
                add_relax=[in_relaxSrC[p4r]]
                add_relax_err=[in_relax_errSrC[p4r]]
                add_atoms_N=[atoms_NSrC[p4r]]
                add_atoms_HN=[atoms_HNSrC[p4r]]
              else:
                add_relax=add_relax+[in_relaxSrC[p4r]]
                add_relax_err=add_relax_err+[in_relax_errSrC[p4r]]
                add_atoms_N=add_atoms_N+[atoms_NSrC[p4r]]
                add_atoms_HN=add_atoms_HN+[atoms_HNSrC[p4r]]
                pass

              in_relaxSrC=remove_position_from_list(in_relaxSrC,p4r)
              in_relax_errSrC=remove_position_from_list(in_relax_errSrC,p4r)
              atoms_NSrC=remove_position_from_list(atoms_NSrC,p4r)
              atoms_HNSrC=remove_position_from_list(atoms_HNSrC,p4r)


              count=count+1
              pass

            in_relaxInG=in_relaxInG+add_relax
            in_relax_errInG=in_relax_errInG+add_relax_err
            atoms_NInG=atoms_NInG+add_atoms_N
            atoms_HNInG=atoms_HNInG+add_atoms_HN

            pass

         if k==1:
            InGInG=simplex_out[0]   # set up initial guess for further iterations
            pass

         # change out of simplex to accommodate different fitting modes
         if mode=='full':
            Fit_res_list=simplex_out[0]
            pass
         if mode=='axial':
            Fit_res_list=[]
            Fit_res_list.append(simplex_out[0][0])
            Fit_res_list.append(simplex_out[0][0])
            Fit_res_list.append(simplex_out[0][1])
            Fit_res_list.append(simplex_out[0][2])
            Fit_res_list.append(simplex_out[0][3])
            Fit_res_list.append(0)
            pass
         if mode=='iso':
            Fit_res_list=[]
            Fit_res_list.append(simplex_out[0][0])
            Fit_res_list.append(simplex_out[0][0])
            Fit_res_list.append(simplex_out[0][0])
            Fit_res_list.append(0)
            Fit_res_list.append(0)
            Fit_res_list.append(0)
            pass

         list_of_Fit_res_list.append([Fit_res_list,simplex_out[1]])

         k=k+1
         pass

    t2=time.time() 
    print("Fitting time:%7.3f" % (t2-t1))

    list_results=[]
  
# calculate parameters of diffusion tensor for all fitted data sets
    for item in list_of_Fit_res_list: 
         results = FitResults()
         fill_in_results(results, item, mode, len(in_relaxInG), dT_e)
         results.rms=rms_relax_function(results.Fit_vector,[in_relax,in_relax_err,atoms_N,atoms_HN,omega,d2,g,mode])/sqrt(results.Number_of_data_points)
         list_results=list_results+[results]
         pass

    for item in list_results:          
        list_results[0].Tau_mean=list_results[0].Tau_mean+item.Tau
        list_results[0].Anis_mean=list_results[0].Anis_mean+item.Anis
        list_results[0].Rhomb_mean=list_results[0].Rhomb_mean+item.Rhomb

        k2=0
        while k2 < 6:
           list_results[0].Fit_vector_mean[k2]=list_results[0].Fit_vector_mean[k2]+item.Fit_vector[k2]
           k2=k2+1
           pass
        
        pass

    nrmk=len(list_results)

    if nrmk > 1:

        list_results[0].Error_Show_Flag=1
    
        list_results[0].Tau_mean=list_results[0].Tau_mean/nrmk; 
        list_results[0].Anis_mean=list_results[0].Anis_mean/nrmk; 
        list_results[0].Rhomb_mean=list_results[0].Rhomb_mean/nrmk;

        k2=0
        while k2 < 6:
           list_results[0].Fit_vector_mean[k2]=list_results[0].Fit_vector_mean[k2]/nrmk
           k2=k2+1
           pass


        for item in list_results:
            list_results[0].Tau_err=(list_results[0].Tau_mean-item.Tau)*(list_results[0].Tau_mean-item.Tau)+list_results[0].Tau_err
            list_results[0].Anis_err=(list_results[0].Anis_mean-item.Anis)*(list_results[0].Anis_mean-item.Anis)+list_results[0].Anis_err
            list_results[0].Rhomb_err=(list_results[0].Rhomb_mean-item.Rhomb)*(list_results[0].Rhomb_mean-item.Rhomb)+list_results[0].Rhomb_err

            k2=0
            while k2 < 6:
                list_results[0].Fit_vector_err[k2]=list_results[0].Fit_vector_err[k2]+(list_results[0].Fit_vector_mean[k2]-
                                                   item.Fit_vector[k2])*(list_results[0].Fit_vector_mean[k2]-item.Fit_vector[k2])
                k2=k2+1
                pass
            pass

        k2=0
        from math import sqrt
        while k2 < 6:
           list_results[0].Fit_vector_err[k2]=sqrt(
               list_results[0].Fit_vector_err[k2]/(nrmk-1))
           k2=k2+1
           pass


        list_results[0].Tau_err=sqrt(list_results[0].Tau_err/(nrmk-1)); 
        list_results[0].Anis_err=sqrt(list_results[0].Anis_err/(nrmk-1)); 
        list_results[0].Rhomb_err=sqrt(list_results[0].Rhomb_err/(nrmk-1));

        pass

    return list_results[0]

def Target_relax_function(in_param, ifp ):
    """target function for fitting
    """

    in_relax=ifp[0]
    in_relax_err=ifp[1]

    if ifp[7]=='full':
       param_to_sim=[in_param[0],in_param[1],in_param[2],in_param[5],in_param[4],in_param[3]] #gamma, betha, alpha
       pass
    if ifp[7]=='axial':
       param_to_sim=[in_param[0],in_param[0],in_param[1],0,in_param[3],in_param[2]] #betha, alpha
       pass
    if ifp[7]=='iso':
       param_to_sim=[in_param[0],in_param[0],in_param[0],0,0,0]
       pass

    # simulation of relaxation data: sim_relax_data coded on C++ level
    # (see diffPot.cc)
    from relaxData import sim_relax_data
    out_sim=sim_relax_data(ifp[2],ifp[3],param_to_sim,ifp[4],ifp[5],ifp[6])     

    Target=0

    k=0
    for item in out_sim:
        tt=in_relax[k]-item[4]
        Target=Target+tt*tt/(in_relax_err[k]*in_relax_err[k])
        k=k+1
        pass

    return Target

def rms_relax_function(in_param, ifp ): # target function for fitting

    from math import sqrt
    in_relax=ifp[0]
    in_relax_err=ifp[1]

    param_to_sim=[in_param[0],in_param[1],in_param[2],in_param[5],in_param[4],in_param[3]] #gamma, betha, alpha

    # simulation of relaxation data: sim_relax_data coded on C++ level
    # (see diffPot.cc)
    from relaxData import sim_relax_data
    out_sim=sim_relax_data(ifp[2],ifp[3],param_to_sim,ifp[4],ifp[5],ifp[6])     

    Target=0


    k=0
    for item in out_sim:
       tt=in_relax[k]-item[4]
       Target=Target+tt*tt
       k=k+1
       pass 

    return sqrt(Target)


class FitResults():
    def __init__(self): # initializing constructor
        self.Error_Show_Flag=0
        self.InG_Tensor=None
        self.Model=None
        self.FitTensor=None
        self.Number_of_data_points=None
        self.DoF=None
        self.Chi_2=None
        self.Chi_2_Dof=None
        self.rms=None
        self.Fit_vector=None
        self.Tau=0
        self.Anis=0
        self.Rhomb=0
        self.Fit_vector_err=[0,0,0,0,0,0]
        self.Tau_err=0
        self.Anis_err=0
        self.Rhomb_err=0
        self.Fit_vector_mean=[0,0,0,0,0,0]
        self.Tau_mean=0
        self.Anis_mean=0
        self.Rhomb_mean=0
        pass
    pass

def diffValuesEuler(dT_e):
    """ given a diffusion tensor (in 100*[inverse ns]) return the
    eigenvalues and Euler angles corresponding to the rotation matrix
    formed from the eigenvectors 
    in a list:

     [ a,b,c, gamma, beta, alpha ]

     The eigenvalues are scaled by a factor of 0.1 such that they have
    units of (10 * [inverse ns]).
    """
    
    from cdsMatrix import SymMatrix_double, eigen
    m = SymMatrix_double([dT_e[0,0],
                          dT_e[0,1],dT_e[1,1],
                          dT_e[0,2],dT_e[1,2],dT_e[2,2]])
    e=list(eigen(m)) 
    e.sort(key=lambda p: p.value())

    ret = [0.01*pair.value() for pair in e]

    (x,y,z) = [pair.vector() for pair in e]
    
    from vec3 import cross
    z=cross(x,y) # enforce right-handed triple of x,y, and z

    from mat3 import Mat3
    rot_Mat=Mat3(x[0],y[0],z[0],
                 x[1],y[1],z[1],
                 x[2],y[2],z[2]) 
   
    from relaxData import Euler_from_rotmat
    (alpha,beta,gamma)=Euler_from_rotmat(rot_Mat)

    ret += [gamma,beta,alpha]

    return ret

def estimate_Dxyz(in_data,
                  omega_1):
    """
    Estimates eigenvalue of diffusion tensor for given set of initial data
     using assumption of isotropic tumbling model
    """

    aveRatio=sum(in_data) / len(in_data)

    from math import sqrt
    Dxyz=omega_1/(6*sqrt((3/2)*aveRatio))
 
    return Dxyz

def showFitResults(results):
    print('\n',"Diffusion tensor model:  ",results.Model)
    print("Number of data points:   ",results.Number_of_data_points)
    print("Chi^2 value:              %0.2f"%results.Chi_2)
    print("Chi^2/degrees of freedom: %0.2f"%results.Chi_2_Dof)
    print("rms value:                %0.2f"%results.rms)
    if results.Error_Show_Flag==0:
         print('\n',"Correlation time:         %0.2f [ns]"%results.Tau)
    else:
         print('\n',"Correlation time:         %0.2f +/-"%results.Tau,"%0.2f [ns]"%results.Tau_err)
         pass
    if results.Model != 'iso':
      if results.Anis >=1:
         print(" Prolate Tensor")
      else:
         print(" Oblate Tensor")
         pass
      if results.Error_Show_Flag==0:
         print("Anisotropy:               %0.2f"%results.Anis)
         print("Rhombicity:               %0.2f"%results.Rhomb)
      else:
         print("Anisotropy:               %0.2f +/-"%results.Anis,"%0.2f"%results.Anis_err)
         print("Rhombicity:               %0.2f +/-"%results.Rhomb,"%0.2f"%results.Rhomb_err)
         pass
      pass

    print('\n',"Diffusion Tensor:",'\n',results.FitTensor)
    print('\n',"Eigen Values, in [10^-7 s^-1]:")

    if results.Error_Show_Flag==0:
      print("    Dx:                   %0.2f"%(100*results.Fit_vector[0]))
      print("    Dy:                   %0.2f"%(100*results.Fit_vector[1]))
      print("    Dz:                   %0.2f"%(100*results.Fit_vector[2]))

    else:
      print("    Dx:                   %0.2f"%(100*results.Fit_vector[0]),"+/- %0.2f"%(100*results.Fit_vector_err[0]))
      print("    Dy:                   %0.2f"%(100*results.Fit_vector[1]),"+/- %0.2f"%(100*results.Fit_vector_err[1]))
      print("    Dz:                   %0.2f"%(100*results.Fit_vector[2]),"+/- %0.2f"%(100*results.Fit_vector_err[2]))

    if results.Model != 'iso':
       print('\n',"Euler angles, in degrees:")
       if results.Error_Show_Flag==0:
         print("    alpha              %7.2f"%(results.Fit_vector[3]))
         print("    betha              %7.2f"%(results.Fit_vector[4]))
         if results.Model != 'axial':
            print("    gamma              %7.2f"%(results.Fit_vector[5]))
       else:
         print("    alpha              %7.2f"%(results.Fit_vector[3]),"+/- %5.2f"%(results.Fit_vector_err[3]))
         print("    betha              %7.2f"%(results.Fit_vector[4]),"+/- %5.2f"%(results.Fit_vector_err[4]))
         if results.Model != 'axial':
            print("    gamma              %7.2f"%(results.Fit_vector[5]),"+/- %5.2f"%(results.Fit_vector_err[5]))
         pass 
       pass     
    print('')
    return

def remove_position_from_list(l_st,p_sn):
    ret_lst=[]
    ret_lst=l_st[:p_sn]+l_st[(p_sn+1):len(l_st)]
    return ret_lst

def fill_in_results(results, in_Fit_res_list, mode, dlina, dT_e):

    Fit_res_list=in_Fit_res_list[0]
    chi_2=in_Fit_res_list[1]

    # making diffusion tensor from fitting results and fill in results
    # structure
    from mat3 import Mat3
    ei_g_m=Mat3(100*Fit_res_list[0],0,0,
                0,100*Fit_res_list[1],0,
                0,0,100*Fit_res_list[2])

    from relaxData import rotmat_from_Euler
    #                     gamma,            beta,           alpha 
    r_m=rotmat_from_Euler(Fit_res_list[5],Fit_res_list[4],Fit_res_list[3])

    r_M=Mat3(r_m[0,0],r_m[0,1],r_m[0,2],
             r_m[1,0],r_m[1,1],r_m[1,2],
             r_m[2,0],r_m[2,1],r_m[2,2])
    
    r_m_t=Mat3(r_m[0,0],r_m[1,0],r_m[2,0],
               r_m[0,1],r_m[1,1],r_m[2,1],
               r_m[0,2],r_m[1,2],r_m[2,2])

    Fitted_diffT=r_M*ei_g_m*r_m_t

    results.InG_Tensor=dT_e
    results.Model=mode
    results.FitTensor=Fitted_diffT
    results.Number_of_data_points=dlina
    results.Chi_2=chi_2
    results.Fit_vector=Fit_res_list
    if mode=='full':
       N_param=6
       pass
    if mode=='axial':
       N_param=4
       pass
    if mode=='iso':
       N_param=1
       pass
    results.DoF=results.Number_of_data_points-N_param
    results.Chi_2_Dof=results.Chi_2/results.DoF

    f_eig=sorted(Fit_res_list[0:3])

    if results.Model != 'iso':
        results.Tau, results.Anis, results.Rhomb = tensParams(*f_eig)
    else:
        results.Anis=2*f_eig[2]/(f_eig[1] + f_eig[0])  # iso
        pass
    pass

    return

def tensParams(Dx,Dy,Dz):
    """
    Given eigenvalues of a diffusion tensor, return a tuple
    containing (tau, anisotropy , rhombicity), where tau is the
    correlation time in ns. 
    """
    tau = 0.5 / (Dx+Dy+Dz)

    if (Dz - Dy) >= (Dy - Dx):
        anis=2*Dz/(Dy + Dx)  # prolate
        rhomb=1.5*(Dy - Dx)/(Dz - 0.5*(Dy + Dx))
    else:
        anis=2*Dx/(Dz + Dy)  # oblate
        rhomb=1.5*(Dy - Dz)/(Dx - 0.5*(Dz + Dy))
        pass

    return (tau,anis,rhomb)



def calc_relax_data(Fr=None, InG=None, res_sel=None):
    """ # Subrutine which simulates NMR relaxation data						
        # INPUT:											
        # Fr is spectrometer frequency in GHz							
        # InG is a vector of diffusion tensor parameters which will be used			
        #     for relaxation rates calculation							
        #     InG=[Dx, Dy, Dz, alpha, betha, gamma]						
        #     for axial mode InG=[Dxy, Dz, betha, alpha]					
        #     for iso mode   InG=[Dxyz]								
        #     Dx, Dy, Dz are in inverse nanoseconds, alpha, betha, and gamma are in degrees
        #     InG is in the same format as results.Fit_vector					
        # res_sel is selection of residues in a format 'resid N:M' where integers N and M	
        #     specify the range of residues used to simulate relaxation data			
        # RETURN:											
        #     simulated_data list of realx_datum() structures containing R1, R2, and R2/R1	
        #     for details see declaration of FitResults() class below.			
    """

    from atomSel import AtomSel, intersection
    if InG==None:
       print('\n', "calc_relax_data: No parameters of diffusion tensor given: exiting", '\n')
       return
       pass
    if Fr==None:
       print('\n', "calc_relax_data: No spectrometer frequency given: exiting", '\n')
       return
       pass
    if res_sel==None:
       print('\n', "calc_relax_data: No residue selection given: exiting", '\n')
       return
       pass



    # preparing input for data simulation routine					

    # dipolar coupling term in inverse nanoseconds					#
    #d2=1.2975;	# for r_NH=1.0202 [A] 						# 
    d2=1.1429		# for r_NH=1.0420 [A]						#
    g=9.8656;		# ratio of giromagnetic ratios for hydrogen and nitrogen	#

    from math import pi
    # omega values are in GHz
    w_0=0			#omega(0)=w_0
    w_N=Fr*(1e-3)*2*pi/g;	#omega(1)=w_N
    w_H=Fr*(1e-3)*2*pi;	#omega(3)=w_H            
    				#omega(2)=w_H-w_N  
    				#omega(4)=w_H+w_N

    omega=[w_0, w_N, w_H-w_N, w_H, w_H+w_N]

    param_to_sim=[InG[0],InG[1],InG[2],InG[5],InG[4],InG[3]] 	# gamma, betha, alpha

    # preparing lists of atoms which will be passed to relaxation data
    # simulation routine 

    atoms_N=[]
    atoms_HN=[]
 


    try:
       sel=intersection(AtomSel(res_sel), AtomSel('name CA'))
    except:
       sel=intersection(res_sel, AtomSel('name CA'))
       pass


    for item in sel:
   
       try:
          tt=AtomSel('(resid '+str(item.residueNum())+' and name N)')
          tt[0].residueName()
          tt=AtomSel('(resid '+str(item.residueNum())+' and name HN)')
          tt[0].residueName()
          atoms_N.append(AtomSel('(resid '+str(item.residueNum())+' and name N)'))
          atoms_HN.append(AtomSel('(resid '+str(item.residueNum())+' and name HN)'))
       except:
          # atom N and/or atom H not present
          pass
              
       pass  

    simulated_data=[]

    from relaxData import sim_relax_data
    out_sim=sim_relax_data(atoms_N,atoms_HN,param_to_sim,omega,d2,g) 

    for item in out_sim:
       sd=realx_datum()
       sd.resid=int(item[0])

       sd.Sel_N=AtomSel('(resid '+str(sd.resid)+' and name N)')
       sd.Sel_HN=AtomSel('(resid '+str(sd.resid)+' and name HN)')

       sd.Ratio=item[4]
       sd.R1=item[1]
       sd.R2=item[2]

       simulated_data.append(sd)

       pass

    return simulated_data

def calcR2R1ratio(Freq=None,
                  Ing=None,
                  Rsel='resid 1:1000'):
    """  Subrutine which simulates NMR relaxation data						
         INPUT:											
         Freq is spectrometer frequency in GHz							
         Ing is a vector of diffusion tensor parameters which will be used			
             for relaxation rates calculation							
             Ing=[Dx, Dy, Dz, alpha, betha, gamma]						
             for axial mode Ing=[Dxy, Dz, betha, alpha]					
             for iso mode   Ing=[Dxyz]								
             Dx, Dy, Dz are in inverse nanoseconds, alpha, betha, and
             gamma are in degrees
             Ing is in the same format as results.Fit_vector					
         res_sel is selection of residues in a format 'resid N:M' where
             integers N and M specify the range of residues used to
             simulate relaxation data			
         RETURN:											
             out_dat list in format [resNum, R2/R1]			
    """
    simulated_data=calc_relax_data(Fr=Freq, InG=Ing, res_sel=Rsel)

    out_dat=[]

    for item in simulated_data:
        out_dat.append([item.resid, item.R2/item.R1])           

    return out_dat





