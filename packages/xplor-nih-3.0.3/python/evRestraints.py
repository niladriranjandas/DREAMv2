
"""
Generate restraints from evolutionary data
"""

def genDist(string,
            atomName="CB",
            glyAtomName="CA",
            dist=7.0,
            probThreshold=0,
            ):
    import os
    if os.path.exists(string):
        string = open(string).read()
        pass
    lines=string.split('\n')
    restraints=""
    numRestraints=0
    from atomSel import AtomSel
    for line in lines:
        try:
            i,j,prob = line.split()
            i=int(i)
            j=int(j)
            prob=float(prob)
        except:
            continue
        if prob<probThreshold:
            continue

        atomNamei=atomName
        atomNamej=atomName
        if AtomSel('resid %d'%i)[0].residueName()=='GLY':
            atomNamei=glyAtomName
        if AtomSel('resid %d'%j)[0].residueName()=='GLY':
            atomNamej=glyAtomName
        
        atomi = "resid %d and name %s" % (i,atomNamei)
        atomj = "resid %d and name %s" % (j,atomNamej)
        
        restraints += "probability %f\n" % prob
        restraints +=  "assign (%s) (%s) %f %f %f" % (atomi,atomj,dist,dist,0)
        restraints += "! probability: %f\n" % prob
        numRestraints+=1
        pass
    
    if numRestraints==0:
        raise Exception("no restraints read from " + string)
    return restraints

def genSS(string,
          phiH=-57,
          dphiH=7,
          psiH=-47,
          dpsiH=7,
          phiE=-127,
          dphiE=20,
          psiE=122,
          dpsiE=20,
          ):
    import os
    if os.path.exists(string):
        string = open(string).read()
        pass
    lines=string.split('\n')[1:]
    restraints=""
    numRestraints=0

    from atomSel import AtomSel
    from selectTools import oneToThree
    for index,line in enumerate(lines):
        line = line.replace('"','')
        try:
            oneChar,pred3,pred8 = line.split()
            resid = index + 1
            resType = AtomSel('resid %d'%resid)[0].residueName()
            if resType != oneToThree(oneChar):
                raise Exception("residue %d type mismatch: %s != %s" %(resid,
                                                                       oneChar,
                                                                       resType))
        except:
            continue
        if pred3==pred8:
            if pred3=="H":
                phi = phiH
                dphi = dphiH
                psi = psiH
                dpsi = dpsiH
            elif pred3=="E":
                phi = phiE
                dphi = dphiE
                psi = psiE
                dpsi = dpsiE
                pass


            cAtomm = "name C and resid %d" % (resid-1)
            nAtom  = "name N and resid %d" % resid
            caAtom = "name CA and resid %d" % resid
            cAtom  = "name C and resid %d" % resid
            nAtomp = "name N and resid %d" % (resid+1)
            restraint = "assign (%s) (%s) (%s) (%s) 1 %f %f 2" %(cAtomm,
                                                                 nAtom,
                                                                 caAtom,
                                                                 cAtom,
                                                                 phi,dphi)
            restraints += restraint + '\n'

            restraint = "assign (%s) (%s) (%s) (%s) 1 %f %f 2" %(nAtom,
                                                                 caAtom,
                                                                 cAtom,
                                                                 nAtomp,
                                                                 psi,dpsi)
            restraints += restraint + '\n'
            numRestraints+=2
            pass
        pass
    

    if numRestraints==0:
        raise Exception("no restraints read from " + string)
    return restraints

def genSS8(string,
           phiH=-57,
           dphiH=7,
           psiH=-47,
           dpsiH=7,
           phiE=-127,
           dphiE=20,
           psiE=122,
           dpsiE=20,
           hThreshold=0.5,
           sThreshold=0.5,
           ):
    #probabilities are in the order of H G I E B T S L(loops), the 8 secondary structure types used in DSSP 
    import os
    if os.path.exists(string):
        string = open(string).read()
        pass
    lines=string.split('\n')[1:]
    restraints=""
    numRestraints=0

    from atomSel import AtomSel
    from selectTools import oneToThree
    for index,line in enumerate(lines):
        line = line.replace('"','')
        try:
            resid,oneChar,pred = line.split()[:3]
            resid = int(resid)
            probs= line.split()[3:]
            print(probs)
            hProb = float(probs[0])
            sProb = float(probs[3])
            
#            resid = index + 1
            resType = AtomSel('resid %d'%resid)[0].residueName()
            if resType != oneToThree(oneChar):
                raise Exception("residue %d type mismatch: %s != %s" %(resid,
                                                                       oneChar,
                                                                       resType))
        except:
            continue
        havePred=False
        print(resid, hProb, sProb)
        if hProb>hThreshold:
            havePred=True
            phi = phiH
            dphi = dphiH
            psi = psiH
            dpsi = dpsiH
        elif sProb>sThreshold:
            havePred=True
            phi = phiE
            dphi = dphiE
            psi = psiE
            dpsi = dpsiE
            pass


        if havePred:
            cAtomm = "name C and resid %d" % (resid-1)
            nAtom  = "name N and resid %d" % resid
            caAtom = "name CA and resid %d" % resid
            cAtom  = "name C and resid %d" % resid
            nAtomp = "name N and resid %d" % (resid+1)
            restraint = "assign (%s) (%s) (%s) (%s) 1 %f %f 2" %(cAtomm,
                                                                 nAtom,
                                                                 caAtom,
                                                                 cAtom,
                                                                 phi,dphi)
            restraints += restraint + '\n'

            restraint = "assign (%s) (%s) (%s) (%s) 1 %f %f 2" %(nAtom,
                                                                 caAtom,
                                                                 cAtom,
                                                                 nAtomp,
                                                                 psi,dpsi)
            restraints += restraint + '\n'
            numRestraints+=2
            pass
        pass
    

    if numRestraints==0:
        raise Exception("no restraints read from " + string)
    return restraints
