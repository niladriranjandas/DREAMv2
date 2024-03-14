
"""
Compute h3JNC values.

References:

 M. Barfield, ``Structural dependencies of interresidue scalar
 coupling (h3)J(NC') and donor $^1$H chemical shifts in the hydrogen
 bonding regions of proteins,'' J. Am. Chem. Soc. 124, 4158-4168 (2002).
 
 H.J. Sass, F.F. Schmid and S. Grzesiek, ``Correlation of
 protein structure and dynamics to scalar couplings across hydrogen
 bonds,'' J. Am. Chem. Soc. 129, 5898-5903 (2007).
""" 

#Grzesiek values
#simple set
a=-357
r0=0
b=3.2
c=1
d=e=f=g=0
#a=-366
#b=3.2
#c=1
#d=0.47
#e=0.7
#f=0.11
#g=0

#Barfield values
a=1
b=3.2
c=-1.31
d=0.62
e=0.92
f=0.14
g=0
r0=1.760

class H3JNC:
    def __init__(s,h,o):
        from selectTools import convertToAtom

        s.h=convertToAtom(h)
        s.o=convertToAtom(o)
        if s.h.atomName().startswith('O'):
            s.h, s.o = (s.o, s.h)
            pass
        s.c=convertToAtom('name C and bondedto atom "%s" %d %s'%
                          (s.o.segmentName(),
                           s.o.residueNum(), s.o.atomName()))
        s.n=convertToAtom('name N and bondedto atom "%s" %d %s'%
                          (s.c.segmentName(),
                           s.c.residueNum(), s.c.atomName()))
#        for a in (s.h,s.o,s.c,s.n):
#            print a.string()
        return

                           
    def calc(s):
        from vec3 import norm
        from bondAngle import BondAngle
        from dihedral import Dihedral
        from math import cos, exp, sin
        
        s.rHO = norm(s.h.pos() - s.o.pos())
        s.angle=BondAngle(s.h,s.o,s.c)
        theta= s.angle.value()
        s.cosTheta=cos(theta)
        s.sinTheta=sin(theta)
        s.cosTheta2=s.cosTheta**2
        s.sinTheta2=1-s.cosTheta2
        s.dihed=Dihedral(s.h,s.o,s.c,s.n)
        rho=s.dihed.value()
        s.cosRho=cos(rho)
        s.sinRho=sin(rho)
        s.cosRho2=s.cosRho**2

        s.drHO = s.rHO - r0
        ret = a * exp(-b * s.drHO) * (c*s.cosTheta2 +
                                      (d*s.cosRho2+
                                       e*s.cosRho+f)*s.sinTheta2) + g
        return ret
    def derivs(s):
        """
        Return the gradients of the h3JNC value with respect to
        h, o, c and n atoms
        """
        from math import cos, exp, sin
        dJdTheta = a * \
                   exp(-b * s.drHO) *(-2*c*s.cosTheta*s.sinTheta +
                                    2*(d*s.cosRho2+e*s.cosRho+f)*
                                      s.sinTheta*s.cosTheta)
        dJdRho = a * exp(-b * s.drHO) * (-2*d*s.cosRho*s.sinRho-
                                         e*s.sinRho)*s.sinTheta2
        dJdrHO = -a * b* exp(-b * s.drHO) * (c*s.cosTheta2 +
                                           (d*s.cosRho2+
                                            e*s.cosRho+f)*s.sinTheta2)
        from vec3 import unitVec
        dh = dJdrHO * unitVec(s.h.pos()-s.o.pos())
        do = -dh

        dTheta_dh, dTheta_do, dTheta_dc = s.angle.derivs()
        dh += dJdTheta * dTheta_dh
        do += dJdTheta * dTheta_do
        dc = dJdTheta * dTheta_dc

        dPhi_dh, dPhi_do, dPhi_dc, dPhi_dn = s.dihed.derivs()
        dh += dJdRho * dPhi_dh
        do += dJdRho * dPhi_do
        dc += dJdRho * dPhi_dc
        dn = dJdRho * dPhi_dn
        
        
        return dh,do,dc,dn
