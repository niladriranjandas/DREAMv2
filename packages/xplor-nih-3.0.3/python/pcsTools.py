################################################################################
# Methods for dealing with PCS data
#
# These methods come from PyParaTools: http://comp-bio.anu.edu.au/mscook/PPT/
# 
# This code was written by Mitchell Stanton-Cook while a PhD student
# in Dr Thomas Huber Lab. No warranty, support or guarantee of correctness 
# is provided. Bugs are likely.
#
# If you find this code useful please acknowledge me.
# 
# DO NOT REMOVE THIS MESSAGE
#
# Mitchell Stanton-Cook (m.stantoncook@gmail.com)
################################################################################


from math     import cos
from math     import acos
from math     import sin
from math     import asin
from math     import radians
from math     import degrees
from math     import pi, sqrt
from mat3 import Mat3

from sys      import exit

from cminpack import leastsq

def getTagResidSels(pcs):
    """
    Generate a list of <m atomSel>.AtomSels for each unique tag residue. This
    can be used to conveniently configure selection pairs for non-bonded
    interactions.

    Here it is assumed that the paramagnetic center(s) are specified by the
    ``a'' member(s) of each restraint's selPairs member.

    This function returns a named tuple with sels and union members, with the
    second argument the union of the list.
    """
    from atomSel import AtomSel
    tagSel=AtomSel("")
    from atomSel import AtomSel, union
    for r in pcs.restraints():
        # for now, assume the tag is specified by sel1
        for pair in r.selPairs():
            tagSel = union(tagSel, pair.a )
        pass
    
    tagResidueSel=AtomSel("")
    for tagAtom in tagSel:
        tagResidueSel = union(tagResidueSel,
                       AtomSel('segid "%s" and resid %d' % (tagAtom.segmentName(),
                                                            tagAtom.residueNum())))
        pass
    from selectTools import getSegsResidues
    segDict = getSegsResidues(tagResidueSel)
    sels=[]
    for segid in list(segDict.keys()):
        sels += [ AtomSel('segid "%s" and resid %d' % (segid, tuple[0])) for
                  tuple in segDict[segid] ]
        pass

    from collections import namedtuple
    Ret = namedtuple('getTagResidSels',['sels','union'])
    return Ret(sels,tagResidueSel)




def ZYZRot(A, B, G):
    """Returns a ZYZ rotation matrix when given 3 Euler angles (in degrees).

    """
    ca = cos(radians(A))
    cb = cos(radians(B))
    cg = cos(radians(G))
    sa = sin(radians(A))
    sb = sin(radians(B))
    sg = sin(radians(G))
    rot = Mat3( -sg*sa + cb*ca*cg , sg*ca + cb*sa*cg , -cg*sb,
                -cg*sa - cb*ca*sg , cg*ca - cb*sa*sg ,  sg*sb,
                 sb*ca            , sb * sa          ,  cb)
    return rot

def RotX90():
    """Returns the rotation matrix for a 90 degree rotation about X.

    """
    rot = Mat3(1, 0, 0,
               0, 0, 1,
               0,-1, 0)
    return rot

def RotY90():
    """Returns the rotation matrix for a 90 degree rotation about Y.

    """
    rot = Mat3(0, 0,-1,
               0, 1, 0,
               1, 0, 0)
    return rot

def RotZ90():
    """Returns the rotation matrix for a 90 degree rotation about Z.

    """
    rot = Mat3( 0, 1, 0,
               -1, 0, 0,
                0, 0, 1)
    return rot

def correctAngle(cosv, sinv):
    """Return an angle in correct quadrant given a cosine and sine of the angle.

    """
    if (cosv <= pi/2.0):
        if (sinv < 0.0):
            sinv = sinv + 2*pi
            return sinv
        else:
            return sinv
    else:
        if(sinv > 0.0):
            return cosv
        else:
            return -1*(cosv) +2*pi

def EulerFromZYZRot(rotMat):
    """Return the 3 Euler angles (A, B, G) in ZYZ from a given rotation matrix.

    """
    try:
        b_c = acos(rotMat[2,2])
        a_c = acos(rotMat[2,0]/sin(b_c))
        g_c = acos(-1*rotMat[0,2]/sin(b_c))
        a_s = asin(rotMat[2,1]/sin(b_c))
        g_s = asin(rotMat[1,2]/sin(b_c))
        aE = correctAngle(a_c, a_s)
        bE = b_c
        gE = correctAngle(g_c, g_s)
        return aE, bE, gE
    except ValueError:
        print(80*'-')
        print("ERROR: One of A,B or G is equal to 0, 90, 180, 270 or 360 deg")
        print("This is a known issue and will be fixed in future versions")
        exit(0)
        print(80*'-')

def ToUTR(Ax, Rh, a,b,g):
    """Reconfigure a X-tensor into Unique Tensor Rep

    """
    toDegrees = True

    Dx = -(Ax)/3.0 + (Rh)/2.0
    Dy = -(Ax)/3.0 - (Rh)/2.0
    Dz = 2.0/3.0*((Ax))
    aDx, aDy, aDz = abs(Dx), abs(Dy), abs(Dz)

    # Determine the UTR case
    if (aDz >= aDy) and (aDy >= aDx):
        toDegrees = False
    if (aDz >= aDx)and (aDx >= aDy):
        g = g + 90.0
        Dy, Dx = Dx, Dy
    if (aDy >= aDz) and (aDz >= aDx):
        Dy, Dz = Dz, Dy
        a,b,g = EulerFromZYZRot((RotX90() * ZYZRot(a,b,g)))
    if (aDy >= aDx) and (aDx >= aDz):
        g = g + 90.0
        Dy, Dx = Dx, Dy
        Dz, Dx = Dx, Dz
        a,b,g = EulerFromZYZRot((RotY90() * ZYZRot(a,b,g)))
    if(aDx >= aDz) and (aDz >= aDy):
        g = g + 90.0
        Dy, Dx = Dx, Dy
        Dy, Dz = Dz, Dy
        a,b,g = EulerFromZYZRot((RotX90() * ZYZRot(a,b,g)))
    if(aDx >= aDy) and (aDy >= aDz):
        Dz, Dx = Dx, Dz
        a,b,g = EulerFromZYZRot((RotY90()* ZYZRot(a,b,g)))

    if toDegrees != False:
        a,b,g = degrees(a), degrees(b), degrees(g)

    #Axial and Rhombic are now in UTR
    Ax = Dz - (Dx + Dy)/2.0
    Rh = Dx - Dy

    # Make Euler angles in 0-360 after manipulation.
    a = a%360
    b = b%360
    g = g%360

    # Do manipulations such that A,B,G in 0-180
    if a >= 0.0 and a < 180.0:
        if b >= 0.0  and  b < 180.0:
            if g >= 0.0  and  g < 180.0:
                pass
            else:
                g = g + 180.0
        else:
            if g >= 0.0 and g < 180.0:
                b =  b + 180.0
                g = -g +180
            else:
                b = b + 180.0
                g = -g
    else:
        if b >= 0 and  b < 180.0:
            if g >=0  and  g < 180.0:
                a =  a + 180.0
                b = -b + 180.0
                g = -g + 180.0
            else:
                a =  a + 180.0
                b = -b + 180.0
                g = -g
        else:
            if g >= 0 and  g < 180.0:
                a =  a + 180.0
                b = -b
                g =  g
            else:
                a = a + 180.0
                b = -b
                g = g + 180.0

    # Important. Fix to 0-360 to get in UTR (really 0-180).
    a = a%360
    b = b%360
    g = g%360

    return Ax,Rh, a,b,g

def buildAxes(tp):

    nAx, nRh, na, nb, ng = ToUTR(tp[3], tp[4], tp[5], tp[6], tp[7])
    rot = ZYZRot(na, nb, ng)
    return (tp[0], tp[1], tp[2],
            rot[0,0], rot[0,1], rot[0,2],
            rot[1,0], rot[1,1], rot[1,2],
            rot[2,0], rot[2,1], rot[2,2])

from vec3 import Vec3, norm
from cdsVector import sqrt

def PCS1M1S(p0, meas, x,y,z):
    """Optimize for the X-tensor given a single model


    """
    xm, ym, zm, ax, rh, a, b, g = p0
    rot  = ZYZRot(a, b, g)
    X   = x - xm
    Y   = y - ym
    Z   = z - zm
    x_t = rot[0,0]*X + rot[0,1]*Y + rot[0,2]*Z
    y_t = rot[1,0]*X + rot[1,1]*Y + rot[1,2]*Z
    z_t = rot[2,0]*X + rot[2,1]*Y + rot[2,2]*Z
    r2 = (x_t*x_t)+(y_t*y_t)+(z_t*z_t)

    r5 = (r2*r2) * sqrt(r2)
    tmp = 1.0/r5
    err = meas - (tmp*(ax * (3.0*z_t*z_t -r2) + rh*1.5*(x_t*x_t - y_t*y_t)))
    return err

def XTensorFit(init, meas, x,y,z,
               tol=1e-8):
    """ Fit the X-tensor by non-linear leastsquares optimization

    init = mx, my, mz, Ax, Rh, A, B, G
    """
    ret = leastsq(PCS1M1S,init,funcArgs=[meas, x, y, z], tol=tol)
    suc = (ret['info'])
    xm,ym,zm,ax,rh,a,b,g = list(ret['x'])
    a,b,g = a%360, b%360, g%360
    return xm,ym,zm, ax, rh, a, b, g, suc


def calcXTensor(vTensor,
                expts=0,
                coords=None,
                weights=None,
                maxDisplacement=15,
                tol=1e-8,
                ):
    """
    Calculate and return the deltachi-tensor using non-linear least-squares
    optimization procedure.

    As side effects, the <m varTensor>.VarTensor and paramagnetic center are
    set to the calculated values.

    This routine requires that the first atom selection for each
    restraint (excluding dummy pseudo atom restraints) be the
    paramagnetic center.

    This is a modified version of varTensorTools.calcTensor.

    This function returns a list of :

       [x1, y1, z1, nAx, nRh, na, nb, ng]

    where x1,y1,z1 are the position of the paramagnetic center, nAx and nRh
    are the axial and rhombic components of the tensor and na, nb, and ng give
    Euler angles (in degrees).

    author: Mitchell Jon Stanton-Cook (m.stantoncook@gmail.com)

    """
    from vec3        import Vec3
    from math        import sqrt
    from random      import randrange

    # To store measured pcs and coordinates
    meas, cx, cy, cz = [], [], [], []

    if not expts: expts = vTensor.expts()

    rdcs = [x for x in expts if x.potName()=='RDCPot1' and x.scale()>1e-15]

    # first call to restraints() in multi-threaded region
    list(map(lambda pot: pot.rawRestraints(), rdcs))

    if not rdcs: return

    ens = rdcs[0].simulation()

    # Might want to check this...                                       #MJSC
    Da = ens.sharedObj(0.)
    R = ens.sharedObj(0.)
    xPos= ens.sharedObj((0,0,0))
    yPos= ens.sharedObj((0,0,0))
    zPos= ens.sharedObj((0,0,0))
    savedX = ens.sharedObj(1)
    metalAtom=None
   
    cDa, cRh = vTensor.Da(), vTensor.Da()*vTensor.Rh()
    from varTensorTools import eulerAngles
    eulerInDeg = eulerAngles(vTensor)[0]
 
    if ens.singleThread():

        coordsSets=[]
        if coords:
            if not weights: weights = [1./len(coords)]*len(coords)
            for i in range(len(coords)):
                coordsSets.append( (coords[i],weights[i]) )
                pass
            pass
        else:
            for cnt in range(ens.size()):
                mem=ens.members(cnt)
                coordsSets.append( (mem.atomPosArr(), mem.weight()) )
                pass
            pass

        nRestraints=0
        for rdc in rdcs:
            for restraint in rdc.rawRestraints():
                #Build the measured for PCSFitting                      #MJSC
                meas.append(restraint.obs())
                if len(restraint.selPairs()) > 1:
                    raise Exception(">1 selPair not supported for PCS.\n" +
                                    "restraint: " +
                                    restraint.name())

                for selPair in restraint.selPairs():
                    if len(selPair.a) > 1:
                        raise Exception(">1 atom not supported for first " +
                                        "selection for PCS.\n" +
                                        "restraint: " +
                                        restraint.name())
                    if metalAtom == None:
                        metalAtom=selPair.a[0]
                    elif metalAtom != selPair.a[0]:
                        raise Exception("different metal center specified.\n" +
                                        "restraint: " + restraint.name() +
                                        "\nprevious: " + metalAtom.string())
                               
                    for bAtom in selPair.b: nRestraints += 1
                pass
            pass

        sim = ens.subSim()
        saveCoords = sim.atomPosArr()


        for rdc in rdcs:
            for restraint in rdc.rawRestraints():
                selPair = restraint.selPairs()[0]
                aAtom = selPair.a[0]
                for bAtom in selPair.b:
                    # Store the spin coordinates                #MJSC
                    for (coords, weight) in coordsSets:
                        bPos = coords[bAtom.index()]
                        cx.append(bPos[0])
                        cy.append(bPos[1])
                        cz.append(bPos[2])
                        pass
                    pass
                pass
            pass

        # Convert extracted data to an array                            #MJSC
        from cdsVector import CDSVector_double
        meas, cx, cy, cz  = (CDSVector_double(meas), CDSVector_double(cx),
                             CDSVector_double(cy), CDSVector_double(cz))
 
        sim.setAtomPosArr( saveCoords )

        if len(meas)<8:
            raise Exception("insufficient restraints found: %d" % len(b))

        # Extract the current saved deltachi-tensor paramaters         #MJSC
        cpCen = metalAtom.pos()
        cxPos = vTensor.xAtom().pos()
        cyPos = vTensor.yAtom().pos()
        czPos = vTensor.zAtom().pos()
        vals = Mat3( cxPos[0]-cpCen[0], cxPos[1]-cpCen[1], cxPos[2]-cpCen[2],
                     cyPos[0]-cpCen[0], cyPos[1]-cpCen[1], cyPos[2]-cpCen[2],
                     czPos[0]-cpCen[0], czPos[1]-cpCen[1], czPos[2]-cpCen[2])
#	angs = EulerFromZYZRot(vals)
        current = [cpCen[0], cpCen[1], cpCen[2],                         
                   cDa, cRh,
                   eulerInDeg[0], eulerInDeg[1], eulerInDeg[2]]
#                   degrees(angs[0]), degrees(angs[1]), degrees(angs[2])]

        print("Current tensor paramaters are: ", current)                
        init   = list(current)
        tensor = XTensorFit(init, meas, cx, cy, cz,tol=tol)

        # Check that we optimize a good solution                        #MJSC
        suc   = tensor[8]
        prange = 10
        count = 0
        maxit = 10
        while suc > 4:
            print("Opt failure. Trying with different metal coords")
            tensor = XTensorFit(init, meas,                   \
                                 cx+randrange(-prange, prange), \
                                 cy+randrange(-prange, prange), \
                                 cz+randrange(-prange, prange))
            suc   = tensor[8]
            count = count+1
            if count > maxit:
                print("Performed max re-optimizations.")
                print("Setting current tensor to previous best tensor")
                tensor = current
                break
        
        # Handle cases where we get large metal position displacement   #MJSC
        count    = 0
        dispx    = (cpCen[0] - tensor[0])**2
        dispy    = (cpCen[1] - tensor[1])**2
        dispz    = (cpCen[2] - tensor[2])**2
        cur_disp = sqrt(dispx + dispy + dispz)
        suc      =  99
        while cur_disp > maxDisplacement:
            print("> maxDisplacement. Trying with different metal coords")
            count = count +1
            if count > maxit:
                print("Performed max re-optimizations.")
                print("Setting current tensor to previous best tensor")
                tensor = current
                break
            while suc > 4:
                print("Opt failure. Trying with different metal coords")
                tensor = XTensorFit(init, meas, \
                                     cx+randrange(-prange, prange), \
                                     cy+randrange(-prange, prange), \
                                     cz+randrange(-prange, prange))
                suc      = tensor[8]
                count    = count +1
                dispx    = (cpCen[0] - tensor[0])**2
                dispy    = (cpCen[1] - tensor[1])**2
                dispz    = (cpCen[2] - tensor[2])**2
                cur_disp = sqrt(dispx + dispy + dispz)
                if count > maxit:
                    break

        # Convert to solution to UTR                                    #MJSC
        nAx, nRh, na, nb, ng = ToUTR(tensor[3], tensor[4],                 \
                                     tensor[5], tensor[6], tensor[7])
        savedX.set([tensor[0], tensor[1], tensor[2], nAx, nRh, na, nb, ng])
        print("Optimized tensor paramaters: ",                             \
                   [tensor[0], tensor[1], tensor[2], nAx, nRh, na, nb, ng])
        
        # Generate principle axes - positions of OO, OO2, X, Y, Z       #MJSC
        
        x1, y1, z1, x2, y2, z2, x3, y3, z3, x4, y4, z4 = buildAxes(tensor)
        noopos = [x1, y1, z1]
        nxpos  = [x2, y2, z2]
        nypos  = [x3, y3, z3]
        nzpos  = [x4, y4, z4]

        Da.set( nAx )
        R.set( nRh/nAx )
        xPos.set( nxpos )
        yPos.set( nypos )
        zPos.set( nzpos )

        pass

    # Not sure about multi-thread???                                    #MJSC
    ens.multiThread()

    if vTensor.oAtom().isValid():
        oPos = vTensor.oAtom().pos()
    else:
        oPos = Vec3(0,0,0)
        vTensor.oAtom().setPos( oPos )
        pass

    metalAtom.setPos(Vec3(x1,y1,z1))

    xPos = Vec3(xPos.get())
    yPos = Vec3(yPos.get())
    zPos = Vec3(zPos.get())

    vTensor.o2Atom().setPos( oPos )
    vTensor.xAtom().setPos( xPos + oPos)
    vTensor.yAtom().setPos( yPos + oPos )
    vTensor.zAtom().setPos( zPos + oPos )

    vTensor.setDa( Da.get() )
    vTensor.setRh( R.get() )

    #update pot values
    for rdc in rdcs: rdc.calcEnergy()

    return [x1, y1, z1, nAx, nRh, na, nb, ng]

def spaceSeparatedToRestraint(inString,
                              paraSel,
                              atomName="HN",
                              residCol=0,
                              pcsCol=1,
                              errCol=None,
                              defaultErr=0.1,
                              segid=None,
                              ):
    """
    Convert string restraint table (inString) consisting of columns
    for resid and pcs values to an Xplor-NIH readable restraint
    table. Column numbers start with 0.  paraSel is an atom selection
    specifying the paramagnetic center.
    """

    lines = [line for line in inString.split('\n') if not line.startswith('#')]

    ret=""
    for line in lines:
        cols = line.split()
        if not cols:
            continue
        try:
            resid = int(cols[residCol])
            pcs = float(cols[pcsCol])
            err = cols[errCol] if errCol!=None else defaultErr
        except IndexError:
            print("Warning: could not read line: %s" % line)
            continue
            
        atomSel = "name %s and resid %d" % (atomName,resid)
        if segid != None:
            atomSel += ' and segid "%s"' % segid
            pass
        ret += "assign () () () () (%s) (%s) %f %f ! %s\n" % (paraSel,atomSel,
                                                              pcs,err,line)
        pass
    return ret
