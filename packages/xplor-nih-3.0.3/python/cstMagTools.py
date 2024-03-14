
"""
Tools to aid in setup/analysis of the chemical shielding tensor potential
term, <m cstMagPot>.CSTMagPot

Reference:
  B.J. Wylie, C.D. Schwieters, E. Oldfield and C.M. Rienstra,
  ``Protein Structure Refinement Using C-13 alpha Chemical
  Shift Tensors,'' J. Am. Chem. Soc. 131, 985-992 (2009).
"""

def create_CSTMagPot(name,
                     surfName,
                     file=None,
                     restraints="",
                     segid=None):
    """
    Create a CSTMagPot term with magnitudes appropriate for the specified
    surfaces referenced by surfName. Restraints can be specified by
    filename or string using the file or restraints arguments, respectively.

    create_CSTMagPot supports a simplified assignment table in which
    name, assign pairs can be replaced with a single statement ca_assign of
    the form:

    ca_assign resid s11 s22 s33
    
    Dihedral angles appropriate for the C-alpha atom phi and psi angles are
    selected by this statement. The segid id argument can be used to specify
    the segment name for the ca_assign statements.
    """

    p1Values, p2Values, surf1, surf2, surf3 = readSurface2D(surfName)

    from cstMagPot import CSTMagPot
    pot = CSTMagPot(name,p1Values,p2Values,surf1,surf2,surf3)

    if file:
        restraints += open(file).read()
        pass

    # implement the ca_assign statement
    from atomSel import AtomSel
    nstr=''
    for line in restraints.split('\n'):
        line = line.strip()
        if line.lower().startswith('ca_assign'):
            (resid,s11,s22,s33) = line.split()[1:5]
            resid=int(resid)
            resname = AtomSel("resid %d" % resid)[0].residueName()
            name = "%s%d" %(resname,resid)
            selCm  = "resid %d and name C"  % (resid-1)
            selN   = "resid %d and name N"  % resid
            selCA  = "resid %d and name CA" % resid
            selC   = "resid %d and name C" % resid
            selNp  = "resid %d and name N"  % (resid+1)
            if segid!=None:
                selCm  += ' and segid "%s"' % segid
                selN   += ' and segid "%s"' % segid
                selCA  += ' and segid "%s"' % segid
                selC   += ' and segid "%s"' % segid
                selNp  += ' and segid "%s"' % segid
                pass
                
            nstr += "\nname %s\n" % name
            nstr += "assign (%s) (%s)\n" %(selCm,selN)
            nstr += "       (%s) (%s)\n" %(selCA,selC)
            nstr += "       (%s) (%s)\n" %(selN,selCA)
            nstr += "       (%s) (%s)  " %(selC,selNp)
            nstr += "%s %s %s" % (s11,s22,s33)
        else:
            nstr += line + '\n'
            pass
        pass

    print(nstr)
    pot.addRestraints( nstr )

    return pot

import os
from os import environ as env

def readSurface2D(surfID):
    """read three 2D surfaces from a format of
    phi1 phi2 s11 s22 s33
    all other nonempty lines should start with a #
    """
    
    datafile = env['CHEMSHIFTTENS'] + "/" + surfID + ".txt"

    if not os.access(datafile,os.F_OK):
        raise Exception("could not access " + datafile +
                        ". Is the surfID specified correctly?")

    file=open(datafile,errors='replace')
    ddata={}
    for line in file.readlines():
        line=line.strip()
        if not line or line.startswith("#"):
            continue
        (p1,p2,s1,s2,s3) = [float(x) for x in line.split()[:5]]
        if not p1 in list(ddata.keys()):
            ddata[p1]={}
        ddata[p1][p2] = (s1,s2,s3)
        pass

    p1Values=list(ddata.keys())
    p1Values.sort(key=lambda x: float(x))

    p2Values=list(ddata[p1Values[0]].keys())
    p2Values.sort(key=lambda x: float(x))

    import cdsMatrix
    s11Mat = cdsMatrix.CDSMatrix_double(len(p1Values),len(p2Values))
    s22Mat = cdsMatrix.CDSMatrix_double(len(p1Values),len(p2Values))
    s33Mat = cdsMatrix.CDSMatrix_double(len(p1Values),len(p2Values))
    ok=True
    from sys import stderr
    for i in range(len(p1Values)):
        p1 = p1Values[i]
        for j in range(len(p2Values)):
            p2 = p2Values[j]
            try:
                s11Mat[i,j] = ddata[p1][p2][0]
                s22Mat[i,j] = ddata[p1][p2][1]
                s33Mat[i,j] = ddata[p1][p2][2]
            except KeyError:
                print(datafile + ": missing data for %s %s" % (p1,p2), file=stderr)
                ok=False
            pass
        pass
    if not ok:
        raise Exception("please fill in missing datapoints")

    p1Values = [float(x) for x in p1Values]
    p2Values = [float(x) for x in p2Values]

    p1Values.append( p1Values[0]+360 )
    p2Values.append( p2Values[0]+360 )

    return p1Values, p2Values, s11Mat, s22Mat, s33Mat,

def calcSlopeIntercept(terms,verbose=False):
    """
    given a sequence of CSTMagPot terms, determine the slope and magnitude
    to minimize (sigma_calcd - sigma_obsd)^2 and set cstOffset and cstScale
    for each of the terms.
    """

    try:
        len(terms)
    except:
        terms = [terms]
        pass
    

    # foreach term, collect values of raw surface (could first set
    # offset=0, scale=1
    # and collect observed values for each of s11, s22, s33
    #
    # do svd on matrix of M = ( 1, sig11_calcd(phi1,psi1) )
    #                         ( 1, sig22_calcd(phi1,psi1) )
    #                         ( 1, sig33_calcd(phi1,psi1) )
    #                         ( 1, sig11_calcd(phi2,psi2) )
    #                         ( 1, sig22_calcd(phi2,psi2) )
    # ( offset, slope) =  v * diag * UT * b
    #
    # where M = u * diag * vT
    # and b is the column vector: ( sig11_obs_1 )
    #                           : ( sig22_obs_1 )
    #                           : ( sig33_obs_1 )
    #                           : ( sig11_obs_2 )
    #                           : ( sig22_obs_2 )

    from simulationTools import getPotTerms
    terms = getPotTerms(terms,'CSTMagPot')

    neq = 0
    for term in terms:
        term.setCSTOffset(0)
        term.setCSTScale(1)
        neq += 3*term.numRestraints()
        pass

    from cdsMatrix import RMat, svd, transpose
    from cdsVector import CDSVector_double
    m = RMat( neq , 2)
    b = CDSVector_double( neq )

    row=0
    for term in terms:
        for r in term.restraints():
            m[row+0,0] = 1
            m[row+0,1] = r.calcd1
            m[row+1,0] = 1
            m[row+1,1] = r.calcd2
            m[row+2,0] = 1
            m[row+2,1] = r.calcd3
            b[row+0] = r.obs1
            b[row+1] = r.obs2
            b[row+2] = r.obs3
            row += 3
            pass
        pass

    svdResults = svd(m,'S')

    v = transpose( svdResults.vT )
    uT = transpose( svdResults.u )
    diag = RMat(2,2,0)
    for i in range(2):
        diag[i,i] = 1.0/svdResults.sigma[i]
        pass

    ( offset, slope) = v * diag * uT * b

    for term in terms:
        term.setCSTOffset(offset)
        term.setCSTScale(slope)
        pass

    return
    

def analyze(potList):
    """perform analysis of CSTMagPot terms and return nicely formatted
    summary"""

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'CSTMagPot')

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
registerTerm(analyze,"Chemical Shielding Tensor Magnitudes","CSTMag",
r"""
For each term the root mean square fit of calculated to experiment is printed,
along with the number of violations. Please see

    B.J. Wylie, C.D. Schwieters, E. Oldfield and C.M. Rienstra,
    "Protein Structure Refinement Using C-13 alpha Chemical
    Shift Tensors," J. Am. Chem. Soc. 131, 985-992 (2009).
""")
    
