
"""
3x3 matrix class

In addition to the classes and functions described here, there is the
function

  rotVector(vec, angle) - return matrix corresponding to rotation about 
                          3-vector vec by angle (in radians)

"""

from cdsMatrix import *
class Mat3(CDSMatrix_double):
    def __init__(s,
                 m11=0,m12=0,m13=0,
                 m21=0,m22=0,m23=0,
                 m31=0,m32=0,m33=0):
        """
        Mat3 has a specialized constructor which allows all matrix elements
        to be enumerated as arguments. example:
        m=Mat3(1,2,3,
               4,5,6,
               7,8,9)
        """
        CDSMatrix_double.__init__(s,3,3)
        s.fromList([[m11,m12,m13],
                    [m21,m22,m23],
                    [m31,m32,m33]])
        return
#    def __str__(s):
#        return s.printAsString()
#    def __repr__(s):
#        return s.printAsString()
    def __mul__(s,x):
        import vec3
#        if type(x)==type(vec3.Vec3(0,0,0)):
#            import cdsVector
#            v=cdsVector.CDSVector_double()
#            x = v.fromList(list(x))
#            pass
        return CDSMatrix_double.__mul__(s,x)
    pass

def rotationAmount(R):
    """given rotation matrix R, determine and return the magnitude of the
    associated rotation, in radians.
    """
    from math import acos, pi
    cosTheta=min(1,max(-1,0.5*(trace(R)-1)))
    return acos(cosTheta)


        
class SymMat3(SymMatrix_double):
    def __init__(s,
                 a=0,b=0,c=0,d=0,e=0,f=0):
        """
        SymMat3 has a specialized constructor which allows lower triangular
        elements to be enumerated as arguments. example:
        m=SymMat3(m11,
                  m12,m22,
                  m13,m23,m33)
        """
           
        SymMatrix_double.__init__(s,[a,b,c,d,e,f])
        return
    pass
    
