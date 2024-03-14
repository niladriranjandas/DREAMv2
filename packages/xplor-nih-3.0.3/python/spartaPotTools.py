"""Tools to aid in setup/analysis of the Sparta potential energy term.

This module provides functions to simplify the creation, manipulation and
analysis of <m spartaPot>.SpartaPot potential terms.
"""

def create_SpartaPot(name,filename=None,
                     restraints="",
                     selection="all",
                     format="plain",
                     saveSet=None,
                     spartaObj=None,
                     annSet="20110324",
                     defaultError=1,
                     useTableErrors=False,
                     ambiguousGlyHA=False,
                     segid=None,
                     verbose=0):
    """

    Create one or more instance of a <m spartaPot>.SpartaPot given a
    restraint table filename or string. The format of this string can
    be plain (internal Xplor-NIH format), 'TALOS', 'PIPP', or 'NMRSTAR'.
    Optionally, a <m sparta>.Sparta object can be supplied. The saveSet
    argument can be specified to choose a specific set of chemical shifts
    from a NMRSTAR table.    

    If useTableErrors is True, errors for restraints are taken from a table
    from Sparta (per residue, per atom type data) averaged over
    psi/phi surfaces. Otherwise, errors not specified in a restraint table are
    set to defaultError

    If ambiguousGlyHA is True, HA1 and HA2 stereo specific chemical
    shifts will be made ambiguous.

    """
    from simulationWorld import SimulationWorld_world
    if verbose==0 and SimulationWorld_world().logLevel()!='none':
        verbose="notSparta"
        pass

    if not spartaObj:
        from spartaTools import create_Sparta
        spartaObj=create_Sparta(selection=selection,annSet=annSet,
                                verbose=verbose if verbose!="notSparta" else 0)
        pass

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
    
    from spartaPot import SpartaPot
    if type(restraints)==type('string'):
        p=SpartaPot(name,spartaObj)
        p.addRestraints(restraints)
        
    else:        
        from potList import PotList
        pl = PotList(name)
        
        atoms = list(restraints.keys())
        atoms.sort()
        for atom in atoms:
            table = restraints[atom]
            
            p = SpartaPot(name+"-%s"%atom,spartaObj)
            p.setDefaultError(defaultError)
            p.addRestraints(table)
            if verbose:
                print("  atom: %4s  num restraints: %d" % (atom,
                                                           p.numRestraints()))
                pass
            if annSet=="20110324":
                p.setScale( forceConstants[atom] )
                pass
            p.setShowAllRestraints(True)

            if useTableErrors:
                from selectTools import threeToOne
                for r in p.restraints():
                    
                    residueName = threeToOne(r.atom.residueName())
                    try:
                        r.setErr( errTable[residueName][atom] )
                    except KeyError:
                        print("Key Error: residueName: ", residueName, end=' ')
                        print("atom: ", atom)
                        raise
                    pass
                pass
            pl.append(p)

            pass
        p = pl
        pass

    return p

forceConstants={}

forceConstants["CA"] = 0.950**(-2) # for the 20110324 annSet
forceConstants["CB"] = 1.140**(-2)
forceConstants["C"]  = 1.082**(-2)
forceConstants["N"]  = 2.529**(-2)
forceConstants["HN"] = 0.510**(-2)
forceConstants["HA"] = 0.256**(-2)

# per-residue type errors averaged over all phi/psi values -
# from Yang Shen

errTable={ #proline N data is bogus
    'A' : {'CA' : 0.7272, 'CB' : 0.9862, 'C' : 1.1803, 'N' : 2.1977, "HA" : 0.2480, 'HN' : 0.4665 },
    'C' : {'CA' : 1.7251, 'CB' : 2.1403, 'C' : 1.3336, 'N' : 2.7183, "HA" : 0.3277, 'HN' : 0.5272 },
    'D' : {'CA' : 0.8881, 'CB' : 1.0826, 'C' : 1.1026, 'N' : 2.4507, "HA" : 0.2307, 'HN' : 0.4883 },
    'E' : {'CA' : 0.8603, 'CB' : 0.9597, 'C' : 1.0043, 'N' : 2.2795, "HA" : 0.2260, 'HN' : 0.4585 },
    'F' : {'CA' : 1.1600, 'CB' : 1.1719, 'C' : 1.1134, 'N' : 2.4639, "HA" : 0.3195, 'HN' : 0.5049 },
    'G' : {'CA' : 0.8390,                'C' : 1.1789, 'N' : 2.5161, "HA" : 0.3095, 'HN' : 0.5341 },
    'H' : {'CA' : 1.2755, 'CB' : 1.5310, 'C' : 1.1312, 'N' : 2.6716, "HA" : 0.3100, 'HN' : 0.5327 },
    'I' : {'CA' : 1.2634, 'CB' : 1.2393, 'C' : 1.1133, 'N' : 2.6296, "HA" : 0.2788, 'HN' : 0.4567 },
    'K' : {'CA' : 0.9389, 'CB' : 0.9840, 'C' : 1.0389, 'N' : 2.3347, "HA" : 0.2359, 'HN' : 0.4620 },
    'L' : {'CA' : 0.8278, 'CB' : 1.1123, 'C' : 1.0904, 'N' : 2.3406, "HA" : 0.2453, 'HN' : 0.4505 },
    'M' : {'CA' : 1.0578, 'CB' : 1.5094, 'C' : 1.1385, 'N' : 2.3856, "HA" : 0.2849, 'HN' : 0.4562 },
    'N' : {'CA' : 0.8898, 'CB' : 1.1543, 'C' : 1.1457, 'N' : 2.4978, "HA" : 0.2395, 'HN' : 0.4951 },
    'P' : {'CA' : 0.8233, 'CB' : 0.9489, 'C' : 1.2279, 'N' : 2.5000, "HA" : 0.2571,               },
    'Q' : {'CA' : 0.8814, 'CB' : 1.1230, 'C' : 0.9872, 'N' : 2.2881, "HA" : 0.2353, 'HN' : 0.4576 },
    'R' : {'CA' : 0.9908, 'CB' : 1.0391, 'C' : 1.0835, 'N' : 2.4214, "HA" : 0.2588, 'HN' : 0.4679 },
    'S' : {'CA' : 0.9949, 'CB' : 0.9625, 'C' : 1.1049, 'N' : 2.5803, "HA" : 0.2613, 'HN' : 0.5011 },
    'T' : {'CA' : 1.0680, 'CB' : 1.0481, 'C' : 1.0670, 'N' : 2.9168, "HA" : 0.2837, 'HN' : 0.4921 },
    'V' : {'CA' : 1.0039, 'CB' : 0.9294, 'C' : 1.0982, 'N' : 2.6730, "HA" : 0.2817, 'HN' : 0.4549 },
    'W' : {'CA' : 1.3104, 'CB' : 1.2560, 'C' : 1.2012, 'N' : 2.6095, "HA" : 0.3195, 'HN' : 0.5507 },
    'Y' : {'CA' : 1.2017, 'CB' : 1.1846, 'C' : 1.1235, 'N' : 2.4618, "HA" : 0.3149, 'HN' : 0.5035 }}




            
            

# Total Deuterium Isotope Shifts for 13Ca and 13Cb Nuclei
# Table 9.1 from Book "Protein NMR Spectroscopy", 2nd Edition, 2006
deuteriumIsotopeShifts={}
d=deuteriumIsotopeShifts
d['A'] = { 'CA' :  -0.680 , 'CB' :  -1.000}  # ALA
d['R'] = { 'CA' :  -0.690 , 'CB' :  -1.110}  # ARG
d['N'] = { 'CA' :  -0.550 , 'CB' :  -0.710}  # ASN
d['D'] = { 'CA' :  -0.550 , 'CB' :  -0.710}  # ASP
d['C'] = { 'CA' :  -0.550 , 'CB' :  -0.710}  # CYS
d['c'] = { 'CA' :  -0.550 , 'CB' :  -0.710}  # cys
d['Q'] = { 'CA' :  -0.690 , 'CB' :  -0.970}  # GLN
d['E'] = { 'CA' :  -0.690 , 'CB' :  -0.970}  # GLU
d['G'] = { 'CA' :  -0.780 , 'CB' :   0.000}  # GLY
d['H'] = { 'CA' :  -0.550 , 'CB' :  -0.710}  # HIS
d['#'] = { 'CA' :  -0.550 , 'CB' :  -0.710}  # HIH
d['h'] = { 'CA' :  -0.550 , 'CB' :  -0.710}  # his
d['I'] = { 'CA' :  -0.770 , 'CB' :  -1.280}  # ILE
d['L'] = { 'CA' :  -0.620 , 'CB' :  -1.260}  # LEU
d['K'] = { 'CA' :  -0.690 , 'CB' :  -1.110}  # LYS
d['M'] = { 'CA' :  -0.690 , 'CB' :  -0.970}  # MET
d['F'] = { 'CA' :  -0.550 , 'CB' :  -0.710}  # PHE
d['P'] = { 'CA' :  -0.690 , 'CB' :  -1.110}  # PRO
d['S'] = { 'CA' :  -0.550 , 'CB' :  -0.710}  # SER
d['T'] = { 'CA' :  -0.630 , 'CB' :  -0.810}  # THR
d['W'] = { 'CA' :  -0.550 , 'CB' :  -0.710}  # TRP
d['Y'] = { 'CA' :  -0.550 , 'CB' :  -0.710}  # TYR
d['V'] = { 'CA' :  -0.840 , 'CB' :  -1.200}  # VAL

del d

def applyIsotopeCorrection(p,correct="computed"):
    """
    Given a <m spartaPot>.SpartaPot object, update the CA and CB chemical
    shift tables with appropriate deuterium chemical shifts. If you have
    measured CA and CB chemical shifts with fully deuterated samples, you
    should call this function before comparing or refining against SPARTA
    predicted chemical shift values.

    The correct argument is used to specify whether the tabulated correction
    is added to the computed chemical shift value (default), or subtracted
    from the observed values.
    """
    import potList
    import sparta
    if correct=="computed":
        if type(p)==type(potList.PotList()):
            p=p[0]
            pass
        sparta=p.sparta()
        rcTab = sparta.RC_Tab
        for (resName,vals) in list(deuteriumIsotopeShifts.items()):
            for (atomName,val) in list(vals.items()):
                rcShift = float(sparta.getRC(resName,atomName))
                rcShift += val
                
                for i in range(1,rcTab.numEntries()+1):
    #                print rcTab.entry(i)
                    if rcTab.entry(i)["RESNAME"]==resName:
                        rcTab.setEntry(i,atomName, str(rcShift))
                        pass
                    pass
                pass
            pass
        pass
    else:
        threeOneMap = {}
        import selectTools
        for (one,three) in list(selectTools.residueMapProtein.items()):
            threeOneMap[three] = one
            pass
        restraints=[]
        for term in p:
            for r in term.rawRestraints():
                restraints.append(r)
                pass
            pass
        for r in restraints:
            resOne = threeOneMap[r.atom.residueName()]
            if r.atom.atomName()=="CA":
                print('correcting CA')
                r.setObs( r.obs() - deuteriumIsotopeShifts[resOne]["CA"] )
            elif r.atom.atomName()=="CB":
                r.setObs( r.obs() - deuteriumIsotopeShifts[resOne]["CB"] )
                pass
            pass
        pass
            
    return

def checkCSReference(p,
                     verbose=False):
    """
    Given a <m spartaPot>.SpartaPot object or a list thereof, check the
    reference for CA/CB, C and HA chemical shifts.

    If verbose is True, status messages will be printed for each test.

    A dictionary of computed referencing corrections whose keys are the
    atom names 'CA', 'CB', 'HA' and 'C'.
    """
    ret = { "CA":0, "CB":0, "HA":0, "C":0 }
    #min num shifts used for fitting and offset checking
    MIN_FIT_DATASIZE=10 
    #correction is needed if offset > OFFSET_THRESH*esitmated_error
    OFFSET_THRESH=5.

    clipper = SShiftClipper()
 
    print("\nChecking for Chemical Shift Outliers ...")

    shifts={}
    p.calcEnergy()
    if hasattr(p,"restraints"):
        restraints = p.restraints()
    else:
        restraints = reduce(lambda x,y: x+y,
                            [list(p.restraints()) for p in p])
        pass
    from sparta import GDB
    import os
    filename = os.path.join(os.environ['SPARTA'],"randcoil.tab")
    rcDB=GDB(filename)
    
    for r in restraints:
        atomName = r.atom.atomName()
        atomName=r.atom.atomName().strip("123")
        resname = r.atom.residueName()
        shift = r.obs()-float(rcDB.getEntry("RESCODE",resname,1)[atomName])
        resid = r.atom.residueNum()
        tmp = clipper.clip(r.atom,shift)
        if tmp!=shift:
            print("%s Secondary Shift: %8.3f    Limit: %8.3f" % \
                  (r.atom.string(), shift, tmp), end=' ')
            if abs(shift) > 2*abs(tmp) and  atomName!="HA":
                print(" !", end=' ')
            else:
                if abs(shift) > 3*abs(tmp) and atomName=="HA":
                    print(" !", end=' ')
                    pass
                pass
            print()
            pass

        key = r.atom.segmentName()+'!'+str(resid)
        if not key in list(shifts.keys()): shifts[key]={}
        shifts[key][atomName] = shift
        shifts[key]["resname"] = resname
        pass
        
 
    CAmCB_CABp=[]
    CAmCB_CABn=[]
    CAp=[]; CBp=[]
    CAn=[]; CBn=[]
    CAmCB_HAn=[] ; CAmCB_HAp=[]
    CAmCB_Cn=[] ; CAmCB_Cp=[]
    HAn=[] ; HAp=[] ; 
    Cn=[] ; Cp=[] ; 
    for entry in list(shifts.values()):
        if not "CA" in entry: continue
        if not "CB" in entry: continue

        resname=entry['resname']
        if resname=="PRO" or resname=="CYS":
            continue

        CA=entry["CA"]
        CB=entry["CB"]

        if CA-CB>=0:
            CAmCB_CABp.append(CA-CB);
            CAp.append(CA);	
            CBp.append(CB);
            if 'HA' in entry:
                HA=entry["HA"]
                CAmCB_HAp.append(CA-CB);
                HAp.append(HA);
                pass
            if 'C' in entry:
                C =entry["C"]
                CAmCB_Cp.append(CA-CB)
                Cp.append(C)
                pass
            pass
        else:
            CAmCB_CABn.append(CA-CB);
            CAn.append(CA);
            CBn.append(CB);	
            if 'HA' in entry:
                HA=entry["HA"]
                CAmCB_HAn.append(CA-CB)
                HAn.append(HA)
                pass
            if 'C' in entry:
                C =entry["C"]
                CAmCB_Cn.append(CA-CB)
                Cn.append(C)
                pass
            pass
        pass

    print("\nChecking Chemical Shift Referencing ...")
 
    #float slope, intercept, rms;
    #float intercept_p, rms_p, intercept_n, rms_n;
 
    # //CA/CB
    (slope, intercept_n, rms_n) = linearFit(CAmCB_CABn, CAn)
    (slope, intercept_p, rms_p) = linearFit(CAmCB_CABp, CAp)
  
    #weighted intercept
    intercept_n -=0.094;
    intercept_p-=0.459;
    intercept=0; rms=0.0;
    from math import sqrt
  
    if len(CAp)>=MIN_FIT_DATASIZE and len(CAn)>=MIN_FIT_DATASIZE:
        # weighted average when both groups have sufficient data
        intercept = (len(CAn)*intercept_n +
                     len(CAp)*intercept_p)/(len(CAn)+len(CAp))
        rms = sqrt( ( rms_n**2*len(CAn)**2 +
                      rms_p**2*len(CAp)**2 )/
                    (len(CAn)+len(CAp))**2 )
    elif len(CAp)<MIN_FIT_DATASIZE and len(CAn)>=MIN_FIT_DATASIZE:
        # when one group is too small
        intercept = intercept_n; rms = rms_n
    elif len(CAn)<MIN_FIT_DATASIZE and len(CAp)>=MIN_FIT_DATASIZE:
        # when one group is too small
        intercept = intercept_p; rms = rms_p
    else:
        if len(CAn)+len(CAp) > 0:
            raise Exception("Not enough CA shifts to check " +
                            "CA/CB chemical shift referencing offset")
        pass
    

    if abs(intercept)>OFFSET_THRESH*rms and rms>0:
        if len(CAn)+len(CAp)>10:
            # do correction when only number of shifts > 10
            print(" Estimated Referencing Offset for CA/CB: ", end=' ')
            print(" %6.3f +/- %6.3fppm Size: %3d" % ( intercept,
                                                        rms,
                                                        len(CAn)+len(CAp) ))
            ret["CA"] = intercept;
            ret["CB"] = intercept;
            pass
        pass
    elif verbose:
        print("  CA and CB referencing check out ok:", end=' ')
        print("Offset:  %6.3f +/- %6.3fppm Size: %3d" % ( intercept,
                                                            rms,
                                                            len(CAn)+len(CAp) ))
        pass
    
    # HA
    if len(HAp)>0:
        (slope, intercept_p, rms_p) = linearFit(CAmCB_HAp, HAp)
        intercept_p+=0.008;
        pass
    if len(HAn)>0:
        (slope, intercept_n, rms_n) = linearFit(CAmCB_HAn, HAn)
        intercept_n-=0.089;
        pass
    intercept=0; rms=0.0;

    if len(HAp)>=MIN_FIT_DATASIZE and len(HAn)>=MIN_FIT_DATASIZE:
        # weighted average when both groups have sufficient data	
        intercept = ( len(HAn)*intercept_n+
                      len(HAp)*intercept_p)/(len(HAn)+len(HAp));
        rms = sqrt(( rms_n**2*len(HAn)**2 + rms_p**2*len(HAp)**2 )/
                     (len(HAn)+len(HAp))**2 )
    elif len(HAp)<MIN_FIT_DATASIZE and len(HAn)>=MIN_FIT_DATASIZE:
        # when one group is too small
        intercept = intercept_n; rms = rms_n;
    elif len(HAn)<MIN_FIT_DATASIZE and len(HAp)>=MIN_FIT_DATASIZE:
        # when one group is too small
        intercept = intercept_p; rms = rms_p;
    elif (len(HAn)+len(HAp))>0:
        print("Not enough HA shifts to check", end=' ')
        print("HA chemical shift referencing offset.")

    if abs(intercept)>OFFSET_THRESH*rms and rms>0:
        if len(HAn)+len(HAp)>10:
            # do correction when only number of shifts > 10
            print(" Estimated Referencing Offset for HA: ", end=' ')
            print(" %6.3f +/- %6.3fppm Size: %3d" % ( intercept,
                                                        rms,
                                                        len(HAn)+len(HAp) ))
            ret["HA"] = intercept;
            pass
        pass
    elif verbose:
        print("  HA referencing checks out ok:", end=' ')
        print("      Offset:  %6.3f +/- %6.3fppm Size: %3d" % ( intercept,
                                                            rms,
                                                            len(HAn)+len(HAp) ))
        
    
 
    (slope, intercept_p, rms_p) = linearFit(CAmCB_Cp, Cp)
    (slope, intercept_n, rms_n) = linearFit(CAmCB_Cn, Cn)
    intercept_n+=0.272;
    intercept_p-=0.226;
    intercept=0; rms=0.0;

    if len(Cp)>=MIN_FIT_DATASIZE and len(Cn)>=MIN_FIT_DATASIZE:
        # weighted average when both groups have sufficient data	
        intercept = (len(Cn)*intercept_n+
                     len(Cp)*intercept_p)/(len(Cn)+len(Cp))
        rms = sqrt(( rms_n**2*len(Cn)**2 + rms_p**2*len(Cp)**2 )/
                   (len(Cn)+len(Cp))**2)
    elif len(Cp)<MIN_FIT_DATASIZE and len(Cn)>=MIN_FIT_DATASIZE:
        # when one group is too small
        intercept = intercept_n; rms = rms_n;
    elif len(Cn)<MIN_FIT_DATASIZE and len(Cp)>=MIN_FIT_DATASIZE:
        # when one group is too small
        intercept = intercept_p; rms = rms_p;
    elif len(Cn)+len(Cp)>0:
        raise exception("Not enough C' shifts to check " +
                        "C' chemical shift referencing offset")

    if abs(intercept)>OFFSET_THRESH*rms and rms>0: 
        if len(Cn)+len(Cp)>10:
            # do correction when only number of shifts > 10
            print(" Estimated Referencing Offset for C': ", end=' ')
            print(" %6.3f +/- %6.3fppm Size: %3d" % ( intercept,
                                                        rms,
                                                        len(Cn)+len(Cp) ))
            ret["C"] = intercept;
            pass
        pass
    elif verbose:
        print("  C' referencing checks out ok:", end=' ')
        print("      Offset:  %6.3f +/- %6.3fppm Size: %3d" % ( intercept,
                                                          rms,
                                                          len(Cn)+len(Cp) ))
    return ret

def linearFit(x,y):
    """
    Fit to a linear equation given input data in vectors x and y.
    Return tuple of slope, intercept, rms.
    """
    nx = len(x)
    ny = len(y)
    
    if nx != ny: raise Exception("input array are of unequal lengths")
    if nx == 0:  raise Exception("input arrays are of zero length")

    from cdsVector import CDSVector_double as vector
    x = vector(x)
    y = vector(y)
    xy = vector(nx)
    xx = vector(ny)
    
    for i in range(nx):    #exclude index 0?
        xy[i] = x[i]*y[i]
        xx[i] = x[i]**2
        pass
    
    from cdsVector import sum
    slope = (sum(xy) - sum(y)*sum(x)/nx) / (sum(xx) - sum(x)*sum(x)/nx)
    intercept = (sum(y) - slope*sum(x))/nx
    
    d  = vector(nx)
    d2  = vector(nx)
    sum_d2 = 0.
    for i in range(nx): # exclude index 0?
        fx = intercept+slope*x[i];
        d[i] = y[i] - fx 
        d2[i] = (y[i] - fx)*(y[i] - fx) 
        sum_d2 += (y[i] - fx)*(y[i] - fx)
        pass

    from math import sqrt
    stdDev = sqrt( sum(d**2) / (nx-1) )
    rms = stdDev/sqrt( nx );
    
    return (slope, intercept, rms)


class SShiftClipper:
    def __init__(s,verbose=False):
        from sparta import GDB
        import os
        filename = os.path.join(os.environ['SPARTA'],"sslimit.tab")
        s.clipTab = GDB(filename)
        s.verbose=verbose
        return
    def clip(s,
             atom,
             shift):
        """
        Check if shift of specified <m atom>.Atom is an outlier.
        Return a clipped chemical shift value
        """
        entry = s.clipTab.getEntry("RESCODE",atom.residueName(),1);
        atomName=atom.atomName().strip("123")
        
        minVal = float( entry[atomName+"_MIN"] )
        maxVal = float( entry[atomName+"_MAX"] )

        if shift<minVal and shift<999:
            if s.verbose:
                print(" Clipped %s Secondary Shift: %8.3f" % (atom.string(),
                                                              shift), end=' ')
                print(" Limit: %8.3f" % minVal)
                pass
            return minVal;
        if shift>maxVal and shift<999:
            if s.verbose:
                print(" Clipped %s Secondary Shift: %8.3f" % (atom.string(),
                                                              shift), end=' ')
                print(" Limit: %8.3f" % maxVal)
                pass
            return maxVal
        return shift;
    pass


def analyze(potList):
    """
    Perform analysis of <m spartaPot>.SpartaPot  terms and return nicely
    formatted summary.
    """

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'SpartaPot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    ret += "%-15s  %6s  %5s  %5s\n" % ("", "RMS(ppm)", "Number",
                                         "Num Viols")

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]

        print(term.showViolations())

        print(term.info())

        ret += "%-15s  %6.3f  %5d  %5d\n" % \
               (name , term.rms(), term.numRestraints(), term.violations() )

        pass
    
    return ret
#
# print print analysis
import simulationTools
from functools import reduce
simulationTools.registerTerm(analyze,
                             "SpartaPot Analysis","SpartaPot",
r"""
For each term report::

  RMS           - root mean square deviation between calculated and target
                  RDC values.
  Number        - total number of chemical shift restraints
  Num Viols     - number of violated chemical shift restraints
""")
        
