"""Routines to calculate ring pucker and distortion.

It currently supports only 5-membered rings.

"""

#
# 5-membered ring info taken from Marzec and Day,
#  J. Biomol. Struct Dyn 10, 1091 (1993).
#

from vec3 import Vec3, dot, cross, unitVec
from math import sin, cos, sqrt, pi, atan2

import atomSel

#def pucker(atoms):
#    """calculate ring pucker coordinates. Returns amplitude, phase as
#    a tuple"""
#    if len(atoms)!=5:
#        raise "pucker: ring size %d not yet supported." % len(atoms)
#
#    pos = map(lambda a:a.pos(),atoms)
#    rc = reduce(lambda x,y:x+y,pos) / len(pos)
#    r = map(lambda q:q-rc,pos)
#
#    A = reduce(lambda y,i:y+sin(2*i*pi/5)*r[i],range(5),Vec3(0,0,0))
#    B = reduce(lambda y,i:y+cos(2*i*pi/5)*r[i],range(5),Vec3(0,0,0))
#
#    zHat = unitVec(cross(A,B))
#    z = map(lambda q:dot(q,zHat),r)
#
#    sinSum = reduce(lambda y,j:y+z[j]*sin(4*j*pi/5),range(5))
#    cosSum = reduce(lambda y,j:y+z[j]*cos(4*j*pi/5),range(5))
#    
#    P = pi/2 + atan2( -sinSum, cosSum )
#    if cos(P-pi/2)>0.9:
#        q = sqrt(2./5) /cos(P-pi/2) * cosSum
#    else:
#        q = -sqrt(2./5) /sin(P-pi/2) * sinSum
#        pass
#
#
#    return (q,P)

from dihedral import Dihedral
from functools import reduce

def pucker(angles):
    """Calculate ring pucker amplitude and phase.

    The amplitude and phase are returned in a (amplitude, phase)
    tuple, and are as defined by:

    C. Altona and M. Sundaralingam, JACS 94, 8205 (1972).  

    The angles argument is a sequence with the endocyclic torsion angles (in
    radians).  In the case of the sugar ring of nucleotides and nucleosides,
    angles should be:

    [nu2, nu3, nu4, nu0, nu1]

    where nu2 is the torsion around the C2'-C3' bond, nu3, around
    C3'-C4', nu4, around C4'-O4', etc.  Alternatively, angles can be a
    sequence with <m atom>.Atom instances representing the covalent
    sequence of ring atoms.  For the sugar ring the atom sequence
    should be:

    [O4', C1', C2', C3', C4']

    """
    if len(angles)!=5:
        raise "pucker: ring size %d not yet supported." % len(angles)

    # Just check the first item in angles (assume the rest are the same type).
    # Also assume input angle values are either int of float.
    if type(angles[0]) not in (int, float):  # then it's an atom.Atom instance
        angles = [Dihedral(*[angles[i] for i in t]).value() for t in [(1,2,3,4),(2,3,4,0),(3,4,0,1),(4,0,1,2),(0,1,2,3)]]

    sin36 = sin(36*pi/180)
    sin72 = sin(72*pi/180)

    phase = atan2(angles[2]+angles[4] - (angles[1]+angles[3]),
                  2*angles[0]*(sin36+sin72))

    amplitude = angles[0] / cos(phase)

    return (amplitude, phase)



def pucker_all(atom_names):
    """Calculate ring pucker amplitude and phase for all rings in the
    structure.  

    This function returns a list with (segid, resid, ((amplitude,
    phase)) tuples, with the ring amplitude and phase in the specified
    segment and residue defined by [C. Altona and M. Sundaralingam,
    JACS 94, 8205 (1972)] (and calculated with the pucker() function).

    atom_names is a sequence with the names of the atoms in the ring (each a
    string).  For nucleotides and nucleosides this sequence should be:

    ["O4'", "C1'", "C2'", "C3'", "C4'"]
    
    """
    result = []
    tag = [(x.segmentName(), x.residueNum()) for x in atomSel.AtomSel('tag')]
    for (segid, resid) in tag:
        atoms = [atomSel.AtomSel('segid "%s" and resid %i and name %s' %
                                 (segid, resid, x)) for x in atom_names]
        if sum([len(x) for x in atoms]) == 5:
            result.append((segid, resid, pucker(atoms)))
    return result
    

def puckerRao(angles):
    """Calculate ring pucker amplitude and phase according to an
    alternate formula.

    The amplitude and phase are retuned in a (amplitude, phase) tuple,
    and are as defined by:

    S.T. Rao, E. Westhof and M. Sundaralingam, Acta Cryst. A37,
    421-425 (1981). 
    
    The angles argument is a sequence with the endocyclic torsion angles (in
    radians).  In the case of the sugar ring of nucleotides and nucleosides,
    angles should be:

    [nu2, nu3, nu4, nu0, nu1]

    where nu2 is the torsion around the C2'-C3' bond, nu3, around
    C3'-C4', nu4, around C4'-O4', etc.  Alternatively, angles can be a
    sequence with <m atom>.Atom instances representing the covalent
    sequence of ring atoms.  For the sugar ring example the atom
    sequence should be:

    [O4', C1', C2', C3', C4']

    """
    # Coded by Guillermo A. Bermejo.
    
    if len(angles) != 5:
        raise "pucker2: ring size %i not supported" % len(angles)

    # Just check the first item in angles (assume the rest are the same type).
    # Also assume input angle values are either int of float.
    if type(angles[0]) not in (int, float):  # then it's an atom.Atom instance
        angles = [Dihedral(*[angles[i] for i in t]).value() for t in [(1,2,3,4),(2,3,4,0),(3,4,0,1),(4,0,1,2),(0,1,2,3)]]

    alphas = [2*pi*angles.index(x)/5 for x in angles] # eq. 2 in ref

    A = 0.4 * sum([x*cos(2*alphas[angles.index(x)]) for x in angles]) # eq. 8 

    B = -0.4 * sum([x*sin(2*alphas[angles.index(x)]) for x in angles]) # eq. 9

    amplitude = (A**2 + B**2)**0.5  # from eq. 10

    phase = atan2(B, A)  # from eq. 11

    return (amplitude, phase)


def puckerRao_all(atom_names):
    """Calculate ring pucker amplitude and phase for all rings in the
    structure. 

    This function returns a list with (segid, resid, ((amplitude,
    phase)) tuples, with the ring amplitude and phase in the specified
    segment and residue defined by [S.T. Rao, E. Westhof and
    M. Sundaralingam, Acta Cryst. A37, 421-425 (1981)] (and calculated
    with the pucker2() function).

    atom_names is a sequence with the names of the atoms in the ring (each a
    string).  For nucleotides and nucleosides this sequence should be:

    ["O4'", "C1'", "C2'", "C3'", "C4'"]
    
    """
    result = []
    tag = [(x.segmentName(), x.residueNum()) for x in atomSel.AtomSel('tag')]
    for (segid, resid) in tag:
        atoms = [atomSel.AtomSel('segid "%s" and resid %i and name %s' %
                                 (segid, resid, x)) for x in atom_names]
        if sum([len(x) for x in atoms]) == 5:
            result.append((segid, resid, pucker2(atoms)))
    return result


def distortion(atoms):
    """calculate (in-plane) ring distortion coordinates. Returns
    amplitude, phase as a tuple."""
    if len(atoms)!=5:
        raise "distortion: ring size %d not yet supported." % len(atoms)

    pos = [a.pos() for a in atoms]
    rc = reduce(lambda x,y:x+y,pos) / len(pos)
    r = [q-rc for q in pos]

    A = reduce(lambda y,i:y+sin(2*i*pi/5)*r[i],list(range(5)),Vec3(0,0,0))
    B = reduce(lambda y,i:y+cos(2*i*pi/5)*r[i],list(range(5)),Vec3(0,0,0))

    zHat = unitVec(cross(A,B))
    z = [dot(q,zHat) for q in r]

    uHat = unitVec(A)
    vHat = cross(zHat,uHat)

    u = [dot(q,uHat) for q in r]
    v = [dot(q,vHat) for q in r]

    psi = -atan2(u[0],v[0])

    x = [u[j]*cos(psi)+v[j]*sin(psi) for j in range(5)]
    y = [-u[j]*sin(psi)+v[j]*cos(psi) for j in range(5)]

    prodSum = reduce(lambda s,j: s+x[j]*y[j], list(range(5)))
    diffSum = reduce(lambda s,j: s+x[j]**2-y[j]**2, list(range(5)))

    gamma = 0.5 * atan2(2*prodSum,diffSum)
    
    if cos(2*gamma)>0.9:
        s2 = 1./cos(2*gamma) * diffSum
    else:
        s2 = 2./sin(2*gamma) * prodSum
        pass

    return (sqrt(s2),gamma)
