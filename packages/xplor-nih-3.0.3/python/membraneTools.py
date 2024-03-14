#============================================================================
# Ez potential for quick embedding of transmembrane domains.
# Reference: 
#   Senes, A.; Chadi, D. C.; Law, P. B.; Walters, R. F. S.; Nanda, V.; DeGrado, W. F. J. Mol. Biol. 2007, 366, 436-448.


from math import exp,pow,fabs
from pyPot import PyPot
from vec3 import Vec3
from atomSel import AtomSel

#Ez parameters for different type of lipids.
#Data from:http://www.brocku.ca/researchers/peter_rand/lipid/default.html

lipidType={'eggPC':(37.0,24.9,25.9), \
           'DLPC':(31.6,27.4,18.9), \
           'DMPC-27':(35.7,26.5,23.1), \
           'DPPC-25':(47.1,16.7,30.4), \
           'DPPC-50':(35.9,31.1,24.1), \
           'DOPC':(35.9,28.1,25.4), \
           'DSPC':(47.7,19.6,31.9), \
           'SOPC':(40.6,24.0,28.7), \
           '16,22PC':(38.3,25.2,27.4), \
           'eggPE':(33.8,19.1,25.0), \
           'eggPE-trans':(37.4,14.6,27.7), \
           'MMPE':(40.8,21.0,29.6), \
           'DMPE':(40.4,22.7,28.8), \
           'POPE':(41.8,11.4,28.7), \
           'DOPS':(39.6,13.9,26.2), \
           'DGDG':(38.8,14.4,22.5), \
           'DOPE|DOPC':(37.7,19.4,28.6,3.0), \
           'POPE|SOPC':(41.7,14.7,30.3,9.0), \
           'POPE|SOPC':(41.5,19.7,30.2,2.0), \
           'DGDG|SOPC':(38.9,18.3,25.1,0.8), \
           'DGDG|POPE':(38.9,15.1,25.0,1.0), \
           'eggPC|chol':(42.0,23.5,33.0,1.0), \
           'DPPC|cholLow':(43.1,22.9,34.2,1.0), \
           'DPPC|cholHigh':(50.8,29.2,34.0,8.0), \
           'eggPC|DAG':(36.6,26.4,25.9,0.125)
        }


class EzPot(PyPot):
    def __init__( self, name, sel="" ):
        PyPot.__init__( self, name)
        self.scale = 10.0
        self.wallforce = 20.0
        self.wall=40.0
        self.T=30.0
        self.E0 = 1.
        self.Zmid = 0.
        self.type = 1
        self.XYcenter=0
        CACB="((name CB and not resName Gly) or (name CA and resName Gly))"
        if not sel:
            self.sel=AtomSel(CACB)
        else:
            self.sel=AtomSel(sel+" and " + CACB)
        self.rr=10
        self.Table={'ALA':(1,-0.29,10.22,4.67),\
                   'ASP':(1,1.19,14.25,8.98),\
                   'GLU':(1,1.30,14.66,4.16),\
                   'PHE':(1,-0.80,19.67,7.12),\
                   'GLY':(1,-0.01,13.86,6.00),\
                   'HIS':(1,0.75,12.26,2.77),\
                   'ILE':(1,-0.56,14.34,10.69),\
                   'LYS':(1,1.66,11.11,2.09),\
                   'LEU':(1,-0.64,17.34,8.61),\
                   'MET':(1,-0.28,18.04,7.13),\
                   'ASN':(1,0.89,12.78,6.28),\
                   'PRO':(1,0.83,18.09,3.53),\
                   'GLN':(1,1.21,10.46,2.59),\
                   'ARG':(1,1.55,9.34,4.68),\
                   'SER':(1,0.10,13.86,6.00),\
                   'THR':(1,0.01,13.86,6.00),\
                   'VAL':(1,-0.47,11.35,4.97),\
                   'TRP':(2,-0.85,11.65,7.20),\
                   'TYR':(2,-0.42,13.04,6.20),\
                   'CYS':(1,0.0,0.01,0.0)
               }

 
    def setScale(self,scale):
        self.scale=scale
    def setWall(self,wall):
        self.wall=wall
    def setThickness(self,T):
        self.T=T
    def setLipid(self,Type):
        par=lipidType[Type]
        self.setWall(par[0]/2+par[1])
        self.setThickness(par[2])
    def setXYCenter(self,s):
        if s:self.XYcenter=s
        
        
    def setPara ( self, type, E0, Zmid, n ) :
        self.type=type
        self.E0=E0
        self.Zmid=Zmid
        self.n=n

# Calculate energy of specific residue
    def resEzEnergy(self, z, resitype) :
        self.setPara(self.Table[resitype][0],self.Table[resitype][1],self.Table[resitype][2],self.Table[resitype][3])
        if fabs(z)<=self.wall:
            zz=fabs(z)*15.0/(self.T/2)
            if self.type==1:
                temp=1+pow(zz/self.Zmid,self.n)
                energy=self.E0/temp
            else:
                temp=-pow(zz-self.Zmid,2)/(2*self.n*self.n)
                energy=self.E0*exp(temp)
        else:
            energy=self.wallforce*fabs(z)*fabs(z)
        return (self.scale*energy) 
    
    def calcEnergy (self) :
        totE = 0.
#        print self.T,self.wall
        for i in self.sel:
            E=self.resEzEnergy(i.pos()[2],i.residueName())
            totE=totE+E
        return totE
      
    def calcEnergyAndDerivList(self,derivs):
        energy=self.calcEnergy()
        for i in self.sel:
            temp=self.resEzDerive(i.pos()[2],i.residueName())
            if self.XYcenter:
                center=self.calcXYCenter()
                derivs[i]=Vec3(center[0],center[1],temp[2])
            else:
                derivs[i]=Vec3(temp)
        return energy
  
#calculate gradient for single residue
    def resEzDerive(self,z,resitype):
        dB=[0.0,0.0,0.0]      
        self.setPara(self.Table[resitype][0],self.Table[resitype][1],self.Table[resitype][2],self.Table[resitype][3])
        if fabs(z)<=self.wall:
            zz=fabs(z)*15.0/(self.T/2)
            sign=fabs(z)/z
            if self.type==1:
                temp=1+pow(zz/self.Zmid,self.n)
                d=-sign*self.E0*self.n*pow(zz,self.n-1)/(pow(self.Zmid,self.n)*pow(temp,2))
            else:
                temp=-pow(zz-self.Zmid,2)/(2*self.n*self.n)
                energy=self.E0*exp(temp)
                d=-sign*energy*(zz-self.Zmid)/pow(self.n,2)
        else:
            d=2.0*self.wallforce*z
        dB[2]=self.scale*round(d,5)
        return dB
    def calcXYCenter(self):
        XX,YY=0.0,0.0
        l=0
        for i in self.sel:
            XX=XX+i.pos()[0]
            YY=YY+i.pos()[1]
            l=l+1
        return [round(XX/l,5),round(YY/l,5)]   

if __name__ == '__main__' :
    mp = EzPot ( "ezPot" )
    print("Ez potential Loaded.")
#========================================================================


def calcTiltAngle(zMin=0,zMax=0,res=""):
    """
    """
    from math import floor, ceil,pi,asin
    from vec3 import cross
    from numpy import linalg, array, mean,std
   
    pos=[]
    zpos=[]
    pos_center=[0,0,0]
    tilt=[]
    zMin=float(zMin)
    zMax=float(zMax)
    bbSel="(name CA or name C or name N or name O or name HN)"
    #Generate atom list based on user's choice. 
    if zMin != zMax:
        sel=AtomSel(res+" AND "+bbSel) if len(res)>0 else AtomSel(bbSel)
        for i in sel:
            if i.pos()[2]>=zMin and i.pos()[2]<=zMax:
                pos.append(i.pos())
                zpos.append(i.pos()[2])
        zMin=max(min(zpos),zMin)
        zMax=min(max(zpos),zMax)
    elif len(res)>0:
        sel=AtomSel(res+" AND "+bbSel)
        for i in sel:
            pos.append(i.pos())
            zpos.append(i.pos()[2])
        zMin=min(zpos)
        zMax=max(zpos)
    else: 
        print("Atom selection failed for tilt angle calculation.")
        return

    #Find geometric center in each bin and fit to 3D line and calculate tilt angles. 
    bin_size=4
    bin_num=int(ceil((zMax-zMin)/bin_size))
    pos_sum=[]
    bin_cnt=[]
    for i in range(bin_num):
        pos_sum.append(Vec3(0.,0.,0.))
        bin_cnt.append(0)
    for i in pos:
        bin_=int(floor((i[2]-zMin)/bin_size))
        pos_sum[bin_]+=i
        bin_cnt[bin_]+=1

    pos_ave=[]
    for i in range(bin_num):
        pos_ave.append(pos_sum[i]/bin_cnt[i])
    if len(pos_ave)>=4:
        for np in range(4,len(pos_ave)+1):
            selectedPoints=pos_ave[:np]
            pos_center=array(selectedPoints).mean(axis=0)
            uu,dd,vv=linalg.svd(selectedPoints-pos_center)
            tilt.append(asin(linalg.norm(cross(Vec3(vv[0]),Vec3([0,0,1]))))*180/pi)
    else:
        print("Selceted helix is too short, recommended >8 residues .")
    
    return round(mean(tilt),2),round(std(tilt),2)

def calcTiltAngleNH(res=""):
    from vec3 import norm,cross
    from numpy import linalg
    from math import asin,pi
    sel=AtomSel(res+" and (name N or name HN)")
    NHvecs,temp=[],[]
    for i in sel:
        if temp:
            if temp[0]==i.residueNum():
                N=temp[2] if temp[1]=="N" else i.pos()
                H=i.pos() if temp[1]=="N" else temp[2]
                NHvecs.append((N-H)/norm(N-H))
                temp=[]
            else: temp=[i.residueNum(),i.atomName(),i.pos()]
        else:
            temp=[i.residueNum(),i.atomName(),i.pos()]
    sumvec=Vec3(0,0,0)
    for v in NHvecs:
        sumvec=sumvec+v
    return round(asin(linalg.norm(cross(sumvec/len(NHvecs),Vec3(0,0,1))))*180/pi,2)


def fitTiltAngleIH(pdb="",alphaSel=""):
    import protocol
    from xplorSimulation import XplorSimulation
    from selectTools import minResid, maxResid
    from atomSelAction import Fit


    bbSel = "(name C or name N or name CA)"
    numRes = 8 # number of residues to use in determination of axis
    xsim = XplorSimulation(clone=False)
    protocol.loadPDB(string=idealHelixModel,simulation=xsim,deleteUnknownAtoms=True)

    bbAlpha = AtomSel("(%s) and (%s)" % (bbSel,alphaSel))
    if not bbAlpha:
        protocol.loadPDB(pdb,deleteUnknownAtoms=True)
        bbAlpha = AtomSel("(%s) and (%s)" % (bbSel,alphaSel))

    startResid = minResid(bbAlpha)
    alphaStretch = AtomSel("resid %d:%d and (%s)" % (startResid,
                                                     startResid+numRes-1,
                                                     bbSel                ))
    fitter = Fit(alphaStretch,
                 AtomSel("resid 1:%d and (%s)" % (numRes, bbSel), xsim ))
    
    axis0 = fitter.rotation() * Vec3(0,0,1)

    from math import acos, pi
    from vec3 import dot
    return round(acos( dot(axis0,[0,0,1]) ) * 180 / pi,2)



def setCenter(selection="not PSEUDO",
              zshift=0):
    """
    Translate the coordinates of all atoms such that the centroid of the
    specified <m atomSel>.AtomSel is set to the origin. The zshift argument
    specifies an offset for the z component of the centroid.
    """
    from atomSelAction import Translate
    sel=AtomSel(selection)
    pos=Vec3(0,0,0)
    print(len(sel))
    for i in sel:
        pos+=i.pos()
    pos=pos/len(sel)
    pos[2]=pos[2]-zshift
    AtomSel("All").apply(Translate(-pos))

def setCenterXY(selection="not PSEUDO"):
    """
    Translate the coordinates of all atoms such that the x- and y- components
    of the centroid of the specified <m atomSel>.AtomSel are set to the origin.
    The z-coordinates are not translated.
    """
    from atomSelAction import Translate
    sel=AtomSel(selection)
    pos=Vec3(0,0,0)
    for i in sel:
        pos+=i.pos()
    pos=pos/len(sel)
    AtomSel("All").apply(Translate(-Vec3(pos[0],pos[1],0)))

def flipAll(Axis="X",theta=180):
    """Needs to be documented and moved somewhere more appropriate.
    """
    from atomSelAction import Rotate
    from mat3 import Mat3, transpose
    from math import sin,cos,pi
    ct=cos(theta*pi/180)
    st=sin(theta*pi/180)
    if Axis=="X":
        AtomSel("all").apply( Rotate(transpose(Mat3(1,0,0,
                                                    0,ct,-st,
                                                    0,st,ct ))))
    if Axis=="Y":
        AtomSel("all").apply( Rotate(transpose(Mat3(ct,0,st,
                                                    0,1,0,
                                                    -st,0,ct ))))
    if Axis=="Z":
        AtomSel("all").apply( Rotate(transpose(Mat3(ct,-st,0,
                                                    st,ct,0,
                                                    0,0,1))) )

#==========coordinate of an ideal helix for tilt anlge fiting=============
idealHelixModel=r'''
ATOM      1  N   ALA     1      -0.802  -1.224 -15.529  1.00  0.00      ALPH
ATOM      2  CA  ALA     1      -0.113  -2.156 -14.653  1.00  0.00      ALPH
ATOM      3  HA  ALA     1      -0.869  -2.766 -14.158  1.00  0.00      ALPH
ATOM      4  CB  ALA     1       0.789  -3.068 -15.486  1.00  0.00      ALPH
ATOM      5  HB1 ALA     1       0.927  -4.016 -14.965  1.00  0.00      ALPH
ATOM      6  HB2 ALA     1       0.326  -3.250 -16.456  1.00  0.00      ALPH
ATOM      7  HB3 ALA     1       1.758  -2.588 -15.630  1.00  0.00      ALPH
ATOM      8  C   ALA     1       0.667  -1.374 -13.593  1.00  0.00      ALPH
ATOM      9  O   ALA     1       0.601  -1.691 -12.407  1.00  0.00      ALPH
ATOM     10  N   ALA     2       1.387  -0.364 -14.061  1.00  0.00      ALPH
ATOM     11  HN  ALA     2       1.434  -0.112 -15.028  1.00  0.00      ALPH
ATOM     12  CA  ALA     2       2.178   0.466 -13.169  1.00  0.00      ALPH
ATOM     13  HA  ALA     2       2.908  -0.176 -12.677  1.00  0.00      ALPH
ATOM     14  CB  ALA     2       2.924   1.525 -13.986  1.00  0.00      ALPH
ATOM     15  HB1 ALA     2       2.221   2.044 -14.636  1.00  0.00      ALPH
ATOM     16  HB2 ALA     2       3.392   2.241 -13.310  1.00  0.00      ALPH
ATOM     17  HB3 ALA     2       3.691   1.042 -14.591  1.00  0.00      ALPH
ATOM     18  C   ALA     2       1.264   1.086 -12.110  1.00  0.00      ALPH
ATOM     19  O   ALA     2       1.581   1.062 -10.922  1.00  0.00      ALPH
ATOM     20  N   ALA     3       0.148   1.624 -12.580  1.00  0.00      ALPH
ATOM     21  HN  ALA     3      -0.102   1.638 -13.547  1.00  0.00      ALPH
ATOM     22  CA  ALA     3      -0.814   2.247 -11.686  1.00  0.00      ALPH
ATOM     23  HA  ALA     3      -0.313   3.072 -11.181  1.00  0.00      ALPH
ATOM     24  CB  ALA     3      -1.981   2.808 -12.503  1.00  0.00      ALPH
ATOM     25  HB1 ALA     3      -2.132   3.856 -12.250  1.00  0.00      ALPH
ATOM     26  HB2 ALA     3      -1.754   2.720 -13.566  1.00  0.00      ALPH
ATOM     27  HB3 ALA     3      -2.886   2.244 -12.278  1.00  0.00      ALPH
ATOM     28  C   ALA     3      -1.271   1.227 -10.642  1.00  0.00      ALPH
ATOM     29  O   ALA     3      -1.312   1.529  -9.451  1.00  0.00      ALPH
ATOM     30  N   ALA     4      -1.603   0.040 -11.127  1.00  0.00      ALPH
ATOM     31  HN  ALA     4      -1.567  -0.198 -12.098  1.00  0.00      ALPH
ATOM     32  CA  ALA     4      -2.055  -1.028 -10.251  1.00  0.00      ALPH
ATOM     33  HA  ALA     4      -2.958  -0.682  -9.748  1.00  0.00      ALPH
ATOM     34  CB  ALA     4      -2.397  -2.263 -11.085  1.00  0.00      ALPH
ATOM     35  HB1 ALA     4      -2.919  -2.990 -10.462  1.00  0.00      ALPH
ATOM     36  HB2 ALA     4      -3.037  -1.974 -11.919  1.00  0.00      ALPH
ATOM     37  HB3 ALA     4      -1.478  -2.708 -11.469  1.00  0.00      ALPH
ATOM     38  C   ALA     4      -0.978  -1.311  -9.202  1.00  0.00      ALPH
ATOM     39  O   ALA     4      -1.276  -1.418  -8.014  1.00  0.00      ALPH
ATOM     40  N   ALA     5       0.252  -1.425  -9.680  1.00  0.00      ALPH
ATOM     41  HN  ALA     5       0.487  -1.337 -10.649  1.00  0.00      ALPH
ATOM     42  CA  ALA     5       1.376  -1.694  -8.800  1.00  0.00      ALPH
ATOM     43  HA  ALA     5       1.191  -2.651  -8.309  1.00  0.00      ALPH
ATOM     44  CB  ALA     5       2.658  -1.805  -9.627  1.00  0.00      ALPH
ATOM     45  HB1 ALA     5       3.128  -2.772  -9.441  1.00  0.00      ALPH
ATOM     46  HB2 ALA     5       2.416  -1.717 -10.685  1.00  0.00      ALPH
ATOM     47  HB3 ALA     5       3.344  -1.007  -9.342  1.00  0.00      ALPH
ATOM     48  C   ALA     5       1.460  -0.596  -7.738  1.00  0.00      ALPH
ATOM     49  O   ALA     5       1.610  -0.886  -6.551  1.00  0.00      ALPH
ATOM     50  N   ALA     6       1.361   0.641  -8.201  1.00  0.00      ALPH
ATOM     51  HN  ALA     6       1.239   0.867  -9.168  1.00  0.00      ALPH
ATOM     52  CA  ALA     6       1.424   1.783  -7.306  1.00  0.00      ALPH
ATOM     53  HA  ALA     6       2.393   1.762  -6.809  1.00  0.00      ALPH
ATOM     54  CB  ALA     6       1.314   3.075  -8.118  1.00  0.00      ALPH
ATOM     55  HB1 ALA     6       0.274   3.395  -8.155  1.00  0.00      ALPH
ATOM     56  HB2 ALA     6       1.917   3.852  -7.648  1.00  0.00      ALPH
ATOM     57  HB3 ALA     6       1.676   2.899  -9.132  1.00  0.00      ALPH
ATOM     58  C   ALA     6       0.320   1.661  -6.253  1.00  0.00      ALPH
ATOM     59  O   ALA     6       0.571   1.845  -5.062  1.00  0.00      ALPH
ATOM     60  N   ALA     7      -0.877   1.353  -6.730  1.00  0.00      ALPH
ATOM     61  HN  ALA     7      -1.072   1.205  -7.699  1.00  0.00      ALPH
ATOM     62  CA  ALA     7      -2.018   1.205  -5.843  1.00  0.00      ALPH
ATOM     63  HA  ALA     7      -2.172   2.158  -5.335  1.00  0.00      ALPH
ATOM     64  CB  ALA     7      -3.267   0.882  -6.668  1.00  0.00      ALPH
ATOM     65  HB1 ALA     7      -3.229   1.426  -7.612  1.00  0.00      ALPH
ATOM     66  HB2 ALA     7      -3.302  -0.189  -6.868  1.00  0.00      ALPH
ATOM     67  HB3 ALA     7      -4.156   1.178  -6.112  1.00  0.00      ALPH
ATOM     68  C   ALA     7      -1.714   0.127  -4.801  1.00  0.00      ALPH
ATOM     69  O   ALA     7      -1.946   0.328  -3.609  1.00  0.00      ALPH
ATOM     70  N   ALA     8      -1.198  -0.991  -5.286  1.00  0.00      ALPH
ATOM     71  HN  ALA     8      -1.012  -1.146  -6.257  1.00  0.00      ALPH
ATOM     72  CA  ALA     8      -0.859  -2.101  -4.412  1.00  0.00      ALPH
ATOM     73  HA  ALA     8      -1.773  -2.423  -3.914  1.00  0.00      ALPH
ATOM     74  CB  ALA     8      -0.317  -3.263  -5.246  1.00  0.00      ALPH
ATOM     75  HB1 ALA     8      -1.150  -3.837  -5.652  1.00  0.00      ALPH
ATOM     76  HB2 ALA     8       0.288  -2.872  -6.064  1.00  0.00      ALPH
ATOM     77  HB3 ALA     8       0.296  -3.907  -4.618  1.00  0.00      ALPH
ATOM     78  C   ALA     8       0.142  -1.626  -3.357  1.00  0.00      ALPH
ATOM     79  O   ALA     8      -0.023  -1.903  -2.169  1.00  0.00      ALPH
ATOM     80  N   ALA     9       1.157  -0.916  -3.827  1.00  0.00      ALPH
ATOM     81  HN  ALA     9       1.284  -0.694  -4.793  1.00  0.00      ALPH
ATOM     82  CA  ALA     9       2.184  -0.399  -2.938  1.00  0.00      ALPH
ATOM     83  HA  ALA     9       2.657  -1.250  -2.448  1.00  0.00      ALPH
ATOM     84  CB  ALA     9       3.239   0.347  -3.757  1.00  0.00      ALPH
ATOM     85  HB1 ALA     9       2.923   1.381  -3.902  1.00  0.00      ALPH
ATOM     86  HB2 ALA     9       4.191   0.330  -3.226  1.00  0.00      ALPH
ATOM     87  HB3 ALA     9       3.355  -0.136  -4.728  1.00  0.00      ALPH
ATOM     88  C   ALA     9       1.533   0.490  -1.877  1.00  0.00      ALPH
ATOM     89  O   ALA     9       1.828   0.363  -0.690  1.00  0.00      ALPH
ATOM     90  N   ALA    10       0.661   1.373  -2.343  1.00  0.00      ALPH
ATOM     91  HN  ALA    10       0.427   1.469  -3.311  1.00  0.00      ALPH
ATOM     92  CA  ALA    10      -0.033   2.283  -1.449  1.00  0.00      ALPH
ATOM     93  HA  ALA    10       0.718   2.891  -0.946  1.00  0.00      ALPH
ATOM     94  CB  ALA    10      -0.947   3.201  -2.263  1.00  0.00      ALPH
ATOM     95  HB1 ALA    10      -1.960   2.798  -2.263  1.00  0.00      ALPH
ATOM     96  HB2 ALA    10      -0.952   4.196  -1.818  1.00  0.00      ALPH
ATOM     97  HB3 ALA    10      -0.580   3.263  -3.288  1.00  0.00      ALPH
ATOM     98  C   ALA    10      -0.803   1.474  -0.402  1.00  0.00      ALPH
ATOM     99  O   ALA    10      -0.737   1.773   0.789  1.00  0.00      ALPH
ATOM    100  N   ALA    11      -1.515   0.467  -0.885  1.00  0.00      ALPH
ATOM    101  HN  ALA    11      -1.564   0.232  -1.855  1.00  0.00      ALPH
ATOM    102  CA  ALA    11      -2.297  -0.387  -0.006  1.00  0.00      ALPH
ATOM    103  HA  ALA    11      -3.031   0.241   0.500  1.00  0.00      ALPH
ATOM    104  CB  ALA    11      -3.036  -1.436  -0.838  1.00  0.00      ALPH
ATOM    105  HB1 ALA    11      -4.110  -1.254  -0.783  1.00  0.00      ALPH
ATOM    106  HB2 ALA    11      -2.710  -1.373  -1.876  1.00  0.00      ALPH
ATOM    107  HB3 ALA    11      -2.816  -2.430  -0.448  1.00  0.00      ALPH
ATOM    108  C   ALA    11      -1.374  -1.016   1.040  1.00  0.00      ALPH
ATOM    109  O   ALA    11      -1.687  -1.016   2.229  1.00  0.00      ALPH
ATOM    110  N   ALA    12      -0.254  -1.536   0.558  1.00  0.00      ALPH
ATOM    111  HN  ALA    12      -0.007  -1.530  -0.411  1.00  0.00      ALPH
ATOM    112  CA  ALA    12       0.716  -2.167   1.436  1.00  0.00      ALPH
ATOM    113  HA  ALA    12       0.223  -3.004   1.929  1.00  0.00      ALPH
ATOM    114  CB  ALA    12       1.885  -2.701   0.606  1.00  0.00      ALPH
ATOM    115  HB1 ALA    12       2.492  -1.866   0.254  1.00  0.00      ALPH
ATOM    116  HB2 ALA    12       2.496  -3.361   1.221  1.00  0.00      ALPH
ATOM    117  HB3 ALA    12       1.499  -3.255  -0.250  1.00  0.00      ALPH
ATOM    118  C   ALA    12       1.167  -1.160   2.497  1.00  0.00      ALPH
ATOM    119  O   ALA    12       1.215  -1.482   3.683  1.00  0.00      ALPH
ATOM    120  N   ALA    13       1.487   0.038   2.032  1.00  0.00      ALPH
ATOM    121  HN  ALA    13       1.446   0.293   1.065  1.00  0.00      ALPH
ATOM    122  CA  ALA    13       1.932   1.095   2.926  1.00  0.00      ALPH
ATOM    123  HA  ALA    13       2.840   0.749   3.420  1.00  0.00      ALPH
ATOM    124  CB  ALA    13       2.260   2.348   2.112  1.00  0.00      ALPH
ATOM    125  HB1 ALA    13       1.562   3.144   2.373  1.00  0.00      ALPH
ATOM    126  HB2 ALA    13       3.278   2.669   2.335  1.00  0.00      ALPH
ATOM    127  HB3 ALA    13       2.175   2.123   1.049  1.00  0.00      ALPH
ATOM    128  C   ALA    13       0.855   1.349   3.982  1.00  0.00      ALPH
ATOM    129  O   ALA    13       1.157   1.438   5.172  1.00  0.00      ALPH
ATOM    130  N   ALA    14      -0.377   1.461   3.510  1.00  0.00      ALPH
ATOM    131  HN  ALA    14      -0.613   1.388   2.542  1.00  0.00      ALPH
ATOM    132  CA  ALA    14      -1.500   1.705   4.399  1.00  0.00      ALPH
ATOM    133  HA  ALA    14      -1.322   2.652   4.906  1.00  0.00      ALPH
ATOM    134  CB  ALA    14      -2.787   1.818   3.578  1.00  0.00      ALPH
ATOM    135  HB1 ALA    14      -2.745   1.124   2.739  1.00  0.00      ALPH
ATOM    136  HB2 ALA    14      -3.642   1.577   4.208  1.00  0.00      ALPH
ATOM    137  HB3 ALA    14      -2.888   2.837   3.202  1.00  0.00      ALPH
ATOM    138  C   ALA    14      -1.571   0.587   5.442  1.00  0.00      ALPH
ATOM    139  O   ALA    14      -1.720   0.854   6.634  1.00  0.00      ALPH
ATOM    140  N   ALA    15      -1.462  -0.641   4.956  1.00  0.00      ALPH
ATOM    141  HN  ALA    15      -1.341  -0.849   3.985  1.00  0.00      ALPH
ATOM    142  CA  ALA    15      -1.512  -1.799   5.831  1.00  0.00      ALPH
ATOM    143  HA  ALA    15      -2.480  -1.796   6.331  1.00  0.00      ALPH
ATOM    144  CB  ALA    15      -1.394  -3.075   4.995  1.00  0.00      ALPH
ATOM    145  HB1 ALA    15      -1.593  -3.942   5.625  1.00  0.00      ALPH
ATOM    146  HB2 ALA    15      -2.117  -3.045   4.181  1.00  0.00      ALPH
ATOM    147  HB3 ALA    15      -0.386  -3.149   4.585  1.00  0.00      ALPH
ATOM    148  C   ALA    15      -0.407  -1.686   6.883  1.00  0.00      ALPH
ATOM    149  O   ALA    15      -0.652  -1.893   8.070  1.00  0.00      ALPH
ATOM    150  N   ALA    16       0.786  -1.360   6.408  1.00  0.00      ALPH
ATOM    151  HN  ALA    16       0.977  -1.193   5.441  1.00  0.00      ALPH
ATOM    152  CA  ALA    16       1.930  -1.216   7.294  1.00  0.00      ALPH
ATOM    153  HA  ALA    16       2.091  -2.176   7.784  1.00  0.00      ALPH
ATOM    154  CB  ALA    16       3.172  -0.867   6.471  1.00  0.00      ALPH
ATOM    155  HB1 ALA    16       2.953  -0.016   5.825  1.00  0.00      ALPH
ATOM    156  HB2 ALA    16       3.992  -0.613   7.141  1.00  0.00      ALPH
ATOM    157  HB3 ALA    16       3.455  -1.724   5.860  1.00  0.00      ALPH
ATOM    158  C   ALA    16       1.618  -0.160   8.356  1.00  0.00      ALPH
ATOM    159  O   ALA    16       1.856  -0.380   9.543  1.00  0.00      ALPH
ATOM    160  N   ALA    17       1.091   0.963   7.891  1.00  0.00      ALPH
ATOM    161  HN  ALA    17       0.901   1.134   6.924  1.00  0.00      ALPH
ATOM    162  CA  ALA    17       0.745   2.053   8.787  1.00  0.00      ALPH
ATOM    163  HA  ALA    17       1.659   2.375   9.287  1.00  0.00      ALPH
ATOM    164  CB  ALA    17       0.190   3.225   7.974  1.00  0.00      ALPH
ATOM    165  HB1 ALA    17      -0.731   3.583   8.434  1.00  0.00      ALPH
ATOM    166  HB2 ALA    17       0.923   4.031   7.953  1.00  0.00      ALPH
ATOM    167  HB3 ALA    17      -0.017   2.895   6.956  1.00  0.00      ALPH
ATOM    168  C   ALA    17      -0.247   1.550   9.836  1.00  0.00      ALPH
ATOM    169  O   ALA    17      -0.082   1.808  11.028  1.00  0.00      ALPH
ATOM    170  N   ALA    18      -1.258   0.840   9.357  1.00  0.00      ALPH
ATOM    171  HN  ALA    18      -1.386   0.634   8.386  1.00  0.00      ALPH
ATOM    172  CA  ALA    18      -2.278   0.298  10.238  1.00  0.00      ALPH
ATOM    173  HA  ALA    18      -2.757   1.136  10.745  1.00  0.00      ALPH
ATOM    174  CB  ALA    18      -3.329  -0.442   9.410  1.00  0.00      ALPH
ATOM    175  HB1 ALA    18      -2.832  -1.063   8.664  1.00  0.00      ALPH
ATOM    176  HB2 ALA    18      -3.929  -1.074  10.065  1.00  0.00      ALPH
ATOM    177  HB3 ALA    18      -3.974   0.280   8.911  1.00  0.00      ALPH
ATOM    178  C   ALA    18      -1.617  -0.604  11.282  1.00  0.00      ALPH
ATOM    179  O   ALA    18      -1.907  -0.499  12.473  1.00  0.00      ALPH
ATOM    180  N   ALA    19      -0.737  -1.469  10.798  1.00  0.00      ALPH
ATOM    181  HN  ALA    19      -0.506  -1.547   9.829  1.00  0.00      ALPH
ATOM    182  CA  ALA    19      -0.032  -2.389  11.674  1.00  0.00      ALPH
ATOM    183  HA  ALA    19      -0.776  -3.013  12.169  1.00  0.00      ALPH
ATOM    184  CB  ALA    19       0.887  -3.285  10.841  1.00  0.00      ALPH
ATOM    185  HB1 ALA    19       0.293  -4.050  10.341  1.00  0.00      ALPH
ATOM    186  HB2 ALA    19       1.404  -2.681  10.095  1.00  0.00      ALPH
ATOM    187  HB3 ALA    19       1.619  -3.761  11.493  1.00  0.00      ALPH
ATOM    188  C   ALA    19       0.733  -1.592  12.733  1.00  0.00      ALPH
ATOM    189  O   ALA    19       0.674  -1.912  13.919  1.00  0.00      ALPH
ATOM    190  N   ALA    20       1.435  -0.571  12.265  1.00  0.00      ALPH
ATOM    191  HN  ALA    20       1.479  -0.317  11.298  1.00  0.00      ALPH
ATOM    192  CA  ALA    20       2.212   0.275  13.157  1.00  0.00      ALPH
ATOM    193  HA  ALA    20       2.953  -0.355  13.649  1.00  0.00      ALPH
ATOM    194  CB  ALA    20       2.938   1.346  12.342  1.00  0.00      ALPH
ATOM    195  HB1 ALA    20       2.747   2.326  12.777  1.00  0.00      ALPH
ATOM    196  HB2 ALA    20       4.010   1.146  12.352  1.00  0.00      ALPH
ATOM    197  HB3 ALA    20       2.576   1.329  11.314  1.00  0.00      ALPH
ATOM    198  C   ALA    20       1.286   0.876  14.216  1.00  0.00      ALPH
ATOM    199  O   ALA    20       1.603   0.859  15.405  1.00  0.00      ALPH
ATOM    200  N   ALA    21       0.161   1.395  13.747  1.00  0.00      ALPH
ATOM    201  HN  ALA    21      -0.089   1.405  12.779  1.00  0.00      ALPH
ATOM    202  CA  ALA    21      -0.813   2.001  14.639  1.00  0.00      ALPH
ATOM    203  HA  ALA    21      -0.326   2.835  15.145  1.00  0.00      ALPH
ATOM    204  CB  ALA    21      -1.988   2.540  13.822  1.00  0.00      ALPH
ATOM    205  HB1 ALA    21      -2.469   3.352  14.368  1.00  0.00      ALPH
ATOM    206  HB2 ALA    21      -1.625   2.911  12.864  1.00  0.00      ALPH
ATOM    207  HB3 ALA    21      -2.709   1.740  13.651  1.00  0.00      ALPH
ATOM    208  C   ALA    21      -1.252   0.973  15.684  1.00  0.00      ALPH
ATOM    209  OT1 ALA    21      -1.297   1.273  16.875  1.00  0.00      ALPH
END 
'''


  
