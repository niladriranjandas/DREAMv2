"""tools to aid assist with the shape potential term

support tools for the <m shapePot>.ShapePot term.
"""

def contribAngle(shape,index):
    """computes an angle between the ensemble member with given index and the
    total aggregate shape tensor.

    Note: not tested.
    """
    from mat3 import Mat3, eigen
    from math import atan2, pi
    from vec3 import dot
    
    vec00 = shape.getContrib(index,0).vector
    vec01 = shape.getContrib(index,1).vector
    vec02 = shape.getContrib(index,2).vector
    vec10 = shape.getTarget(0).vector
    vec11 = shape.getTarget(1).vector
    vec12 = shape.getTarget(2).vector
    
    mat = Mat3(dot(vec00,vec10), dot(vec00,vec11), dot(vec00,vec12),
               dot(vec01,vec10), dot(vec01,vec11), dot(vec01,vec12),
               dot(vec02,vec10), dot(vec02,vec11), dot(vec02,vec12))

    e=eigen(mat)

    return atan2(e[1].value().imag,e[1].value().real)*180/pi


def analyze(potList):
    "perform analysis of ShapePot terms and return nicely formatted summary"

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'ShapePot')

    if not potList: return ret

    instanceNames = [x.instanceName() for x in potList]
    instanceNames.sort()

    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]
        print("ShapePot term:", name)
        print(term.info())
        pass

    orientTerms=[]
    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]
        if term.orientScale()>0.: orientTerms.append(term)
        pass
            

    #orientation analysis
    if orientTerms:
        import atomSelAction
        def trace(m):
            ret=0
            for i in range(m.rows()):
                ret += m[i,i]
                pass
            return ret
        
        ret += "orientational differences\n"

        ret += "%-11s  %7s  %8s  %8s  %9s  %8s\n" % \
               (" " , "size", "tens ang", "tens ave",
                "rigid ang", "rig. ave")

        for term in orientTerms:
            name = term.instanceName()
            esim = term.simulation()
            
            print(term.showVectors())
            aveTensAngle=0.
            if term.targetType()=='pairwise':
                tensAngle = term.rotation(0)
                #NOTE: this average is not a proper ensemble average-
                #  the correct weight is not used.
                #Also, this ave. is not the same as the pot. ave. in
                #  ShapePot::energyMaybeDerivs2
                for index in range(esim.size()):
                    aveTensAngle += term.rotation(index)
                    pass
                if esim.size()>1:
                    aveTensAngle /= esim.size()-1
                    pass
            else:
                tensAngle = term.rotation()
                for index in range(esim.size()):
                    aveTensAngle += esim.weight(index) * term.rotation(index)
                    pass
                pass

            
            #    rigid-body rmsd orientation diff to target
            #    ensemble ave rigid-body rmsd orientation diff to target
            from ensembleSharedObj import SharedObj
            sharedAngle = []
            sharedAxis = []
            for index in range(esim.size()):
                sharedAngle.append( SharedObj(esim) )
                sharedAxis.append( SharedObj(esim) )
                pass
            fit = atomSelAction.Fit( esim.members(0).atomPosArr(),
                                     term.atomSel() )
            fit.init_unused("all")
            angleTol=1e-8
            from math import acos, pi
            rotMat = fit.rotation()
            cosTheta=0.5*trace(rotMat)-0.5
            cosTheta=min(cosTheta,1-1e-8); cosTheta=max(cosTheta,-1+1e-8)
           
            rigidAngle = acos(cosTheta) *180/pi
            from vec3 import Vec3, unitVec
            from cdsVector import CDSVector_double as CDSVector
            from cdsMatrix import RMat as CDSMatrix
            from cdsMatrix import inverse
            if rigidAngle<angleTol:
                rigidVec = Vec3(0,0,0)
            #elif rigidAngle<(180-angleTol):
            #    rigidVec = unitVec(Vec3(rotMat[2,1]-rotMat[1,2],
            #                            rotMat[0,2]-rotMat[2,0],
            #                            rotMat[1,0]-rotMat[0,1]))
            else:  #theta=180
                if 1-rotMat[2,2]<angleTol:
                    rigidVec = Vec3(0,0,1)
                else:
                    m=CDSMatrix(2,2)
                    m.fromList([(rotMat[0,0]-1,rotMat[0,1]),
                                (rotMat[1,0]  ,rotMat[1,1]-1)])
                    v=CDSVector(2)
                    v.fromList([rotMat[0,2], rotMat[1,2]])
                    r = -inverse(m) * v
                    rigidVec = unitVec(Vec3(r[0],r[1],1))
                    pass
                pass
            sharedAngle[esim.member().memberIndex()].set(rigidAngle)
            sharedAxis[esim.member().memberIndex()].set(tuple(rigidVec))

            print("rotation matrix:\n",rotMat)

            print("term",name,"rigid body rotation:")
            print(" mem  angle  about vector: ")
            for index in range(esim.size()):
                print("  %2d %7.2f" % (index,sharedAngle[index].barrierGet()), end=' ')
                print(" (%7.4f, %7.4f, %7.4f)" % sharedAxis[index].barrierGet())
                pass

            #note that this average is not a proper ensemble average
            aveRigidAngle=0
            for index in range(esim.size()):
                aveRigidAngle += sharedAngle[index].barrierGet()
                pass
            if esim.size()>1:
                aveRigidAngle /= esim.size()-1
                pass

            ret += "%-11s  %7d  %8.3f  %8.3f  %9.3f  %8.3f\n" % \
                   (name,term.atomSel().size(),tensAngle,aveTensAngle,
                    rigidAngle,aveRigidAngle)
            
        pass
    
    sizeTerms=[]
    for name in instanceNames:
        term = [x for x in potList if x.instanceName()==name][0]
        if term.sizeScale()>0.: sizeTerms.append(term)
        pass
            

    #magnitude analysis
    from math import sqrt
    if sizeTerms:
        ret += "tensor magnitude terms\n"

        ret += "%-11s  %7s %7s %7s %7s %7s %7s %7s\n" % \
               (" " , "size", "X", "X diff", "Y", "Y diff", "Z", "Z diff")

        for term in sizeTerms:
            name = term.instanceName()
            print(term.showValues())

            ave=[0,0,0]
            rmsd=[0,0,0]
            for i in range(esim.size()):
                for a in range(3): ave[a] += term.getContrib(i,a).value
                pass
            for a in range(3): ave[a] /= esim.size()

            for i in range(esim.size()):
                for a in range(3):
                    rmsd[a] += (term.getContrib(i,a).value-ave[a])**2
                    pass
                pass
            for a in range(3):
                rmsd[a] = sqrt( rmsd[a]/esim.size() )
                pass
            
            index = term.simulation().member().memberIndex()
                
            values = [term.getContrib(index,a).value for a in (0,1,2)]
            ret += "%-11s  %7d %7.2f %7.3f %7.2f %7.3f %7.2f %7.3f\n" % \
                   (name,term.atomSel().size(),
                    values[0],rmsd[0],values[1],rmsd[1],values[2],rmsd[2])
            pass
        pass
    
    return ret


from simulationTools import registerTerm
registerTerm(analyze,"Shape Tenso Analysis","Shape",
r"""
The following quantities are printed for each <m shapePot>.ShapPot term:
  size       - number of atoms in the atom selection
  X/Y/Z      - 3 eigenvalues of the shape tensor
  X/Y/Z diff - the difference of the corresponding eigenvale between different
               ensemble members.
""")
