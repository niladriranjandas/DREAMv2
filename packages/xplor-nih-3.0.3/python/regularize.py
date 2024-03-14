"""
Functions to modify atomic coordinates such that covalent restraints (bond,
angle, improper) are satisfied. 
"""
from functools import reduce

# decide how useful this code really is
# it is used by eginput/dna_refi/ensemble.py
def covalentMinimize(sel=0,
                     numSteps=100,
                     dEpred=1.   ):
    """Perform gradient optimization including only covalent terms (bonds,
    angles, impropers)

    This function resets the XPLOR constraints interaction settings.
    """
    
    from atomSel import AtomSel, notSelection
    if isinstance(sel,str): sel = AtomSel(sel)
    if not sel: sel = AtomSel("known")

    import simulation
    oldCurrentSimulation = simulation.currentSimulation()
    sim = sel.simulation()

    from xplorSimulation import getXplorSimulation
    xSim = getXplorSimulation(sim)
    simulation.makeCurrent( xSim )

    xplorCommand = xSim.fastCommand
    outputState=xSim.disableOutput()
    prevPrintFile = xplorCommand("set print=off end","prev_print_file")[0]

    savedMass=xSim.atomMassArr()
    indices=sel.indices()
    from atomSelAction import SetProperty
    sel.apply( SetProperty("mass",1.) )
    notSelection(sel).apply( SetProperty("mass",0.) )


    from xplorPot import XplorPot
    pots = [XplorPot(name) for name in ["bond","angl","impr"]]

    xplorCommand("""constraints
                      interaction (attr mass>0)
                                  (attr mass=0)
                      interaction (attr mass>0)
                                  (attr mass>0)
                    end""")

    from ivm import IVM
    ivm=IVM()
    ivm.setVerbose(0) #don't print info messages
    ivm.fix("attr mass = 0")
    from protocol import cartesianTopology
    cartesianTopology(ivm,sel=sel)

    
    from atomSelAction import SetProperty
    AtomSel("all",sel.simulation()).apply(SetProperty("mass",100.))

    from protocol import initMinimize
    initMinimize(ivm,potList=pots,numSteps=numSteps,
                 dEPred=dEpred,printInterval=10)

    ivm.run()

    xplorCommand("set print=%s end"%prevPrintFile)
    xSim.enableOutput(outputState)

    def restoreContext():
        xplorCommand("constraints inter (all) (all) end")
        xSim.setAtomMassArr(savedMass)
        oldCurrentSimulation.sync()
        simulation.makeCurrent( oldCurrentSimulation )
        pass

    restoreContext()
    return

def fixupCovalentGeom(sel="known",
                      useVDW=0,
                      useDynamics=1,
                      dynRatio=5,
                      maxIters=40,
                      verbose=0,
                      maxViols  =0,
                      bondTol   =0.01,
                      angleTol  =2.,
                      torsionTol=2.,
                      extraTerms=[],
                      suppressExceptions=False
                      ):
    """Given the atoms in selection sel, perform, minimization and (optionally)
    dynamics to remove bond, angle, and improper violations - so that their
    total number is less than or equal to maxViols.

    If useVDW is set, the nonbonded term will be used 1/4 of the time.
    if useDynamics is set, dynamics will be used 1/dynRatio of the time.
    maxIters is the total maximum number of iterations.

    Tolerances for bond, bond angle and improper torsion angles can also be
    specified using optional arguments.

    Additional <m pot>.Pot terms to satify may be specified with the
    extraTerms argument. Note that use of this argument probably will not work,
    due to internal use of old XPLOR dynamics and minimization engines.

    If this function is not successful in satisfying all covalent restraints,
    it will throw the exception regularize.CovalentViolation, if
    suppressExceptions is False. The structure is set to that with the
    minimal number of violations, and if more than one configuration
    has that number of violations, the structure with the minimal
    energy is used. 

    This function resets the XPLOR constraints interaction settings.

    If verbose>0, intermediate status is printed. Increase verbose for more
    verbosity up to a maximum value of 6.

    The return value is a dictionary with the following keys:

          'numViols' - number of violations
          'viols'    - list of violations of each term, bond, angle, improper
          'coords'   - coordinates corresoponding to min numViols
         
    """
    from selectTools import convertToAtomSel
    sel = convertToAtomSel(sel)
    
    #leave pseudoatoms alone
    from atomSel import AtomSel, intersection
    sel = intersection(sel, AtomSel("not pseudo",sel.simulation()))

    import simulation
    oldCurrentSimulation = simulation.currentSimulation()
    sim = sel.simulation()

    from xplorSimulation import getXplorSimulation
    xSim = getXplorSimulation(sim)
    simulation.makeCurrent( xSim )

    xplorCommand = xSim.fastCommand
    outputState=xSim.disableOutput()
    if verbose<=3:
        prevPrintFile = xplorCommand("set print=off end","prev_print_file")[0]
        pass
    if verbose>4:
        xplorCommand("set print=OUTPUT end")
        xplorCommand("set mess=on end")
        pass
    prevTorsionTol= float(xplorCommand("","dihed_tol")[0])
    xplorCommand("set toltorsion 0.9999 end")

    from atomSelAction import SetProperty
    savedMass=sim.atomMassArr()
    sel.apply( SetProperty("mass",-1.) )
    #wrong: this re-evaluates the selection - doesn't work in general...l
    
    AtomSel("attr mass>=0",sim).apply( SetProperty("mass",0.) )
    
    def covalentViols(xplorCmd):

        ret = [int(xplorCmd("print threshold %f %s end"%(name_thresh[1],name_thresh[0]),
                                   "violations")[0]) for name_thresh in (("bonds",bondTol),
                   ("angles",angleTol),
                   ("impropers",torsionTol))]
        for term in extraTerms:
            ret.append( term.violations() )
        return ret
    
    import random

    AtomSel("attr mass # 0",sim).apply(SetProperty("mass",100.))

    actions = ["min"]*(dynRatio-1)
    if useDynamics:
        sel.apply(SetProperty("fric",10.))
        actions.append("dyn")
        pass


    if useVDW:
        from protocol import initNBond
        initNBond(nbxmod=2,cutnb=6.5,tolerance=2.0,rcon=0.01)
        pass
    
    from ivm import IVM
    from protocol import initMinimize, initDynamics, cartesianTopology
    minState={}
    def fixupLoop(ivm,maxIters,verbose,maxViols,minState):
        oldPots=None
        lSim=getXplorSimulation(ivm.simulation)

        from xplorPot import XplorPot
        potCombinations = [[XplorPot(name) for name in list] for list in [["bond"],
                               ["bond","angl"],
                               ["bond","impr"],
                               ["bond","angl","impr"]]]

        if useVDW:
            potCombinations.append( potCombinations[-1]+[XplorPot('VDW')] )
            pass
        
        for term in extraTerms:
            potCombinations.append( potCombinations[-1]+[term] )
            if useVDW:
                potCombinations.append( potCombinations[-3]+[term] )
                pass
            pass

        pots = potCombinations[3] #first time try all terms
        action = actions[0]       # and minimization
        for iter in range(0,maxIters):
    
            viols=covalentViols(lSim.fastCommand)
            if verbose>0: print("iter", iter, "  violations:",viols)
            numViols=reduce(lambda x,y:x+y,viols)
    
    
            if numViols<minState['numViols']:
                minState['numViols']=numViols
                minState['viols']   =viols
                minState['coords']  =lSim.atomPosArr()
                pass
    
    
            if numViols<=maxViols: break

            if oldPots: action = random.choice(actions)
            
            #
            # randomly choose which minimization scheme to apply
            #
            while 1 and oldPots:
                pots = random.choice(potCombinations)
                if pots!=oldPots: break
                pass
            oldPots=pots
    
            import potList
            potList = potList.PotList()
            for pot in pots:
                potList.add( pot )
                pass
    
            if verbose>0:
                print("\t action: %s  " % action, end=' ')
                print("potential terms:", end=' ')
                print([p.instanceName() for p in potList])
                pass
    
            printInterval=0
            if verbose>2:
                printInterval=1
                pass
    
            if verbose>3:
                ivm.setVerbose(                 
                      ivm.printResetCM         |
                      ivm.printVelFromCartCost |
                      ivm.printTemperature     |
                      ivm.printEnergy          |
                      ivm.printCMVel           |
                      ivm.printStepDebug       |
                      ivm.printStepInfo        |
                      ivm.printNodeDef         |
                      ivm.printLoopDebug       |
                      ivm.printLoopInfo )
                pass
            if verbose>4:
                ivm.setVerbose(ivm.verbose() | ivm.printNodePos)
                pass
            if verbose>5:
                ivm.setVerbose(ivm.verbose() | ivm.printNodeTheta)
                pass
            if action=="min":
                energy=max(1.0,abs(potList.calcEnergy()))
                initMinimize(ivm,numSteps=500,printInterval=printInterval,
                             dEPred=energy/10.,
                             potList=potList)
            else:
                bathTemp=1000
                initDynamics(ivm,numSteps=500,stepsize=0.001,
                             printInterval=printInterval,
                             bathTemp=bathTemp,
                             potList=potList)
                from atomAction import randomizeVelocities
                randomizeVelocities( bathTemp, "attr mass>0")
                pass
            import _publicIVM
            try:
                ivm.run()
            except _publicIVM.IVMError:
                print('ivm error. Continuing...')
                pass
            if verbose>0:
                print("number of IVM dof: %d" % ivm.dof())
            pass
        return minState
    def fixupLoop_xplor(ivm,maxIters,verbose,maxViols,minState):
        """ fixupLoop using XPLOR's Cartesian dynamics/minimizer instead of
        the IVM"""
        oldPots=None
        lSim=getXplorSimulation(ivm.simulation)
        lXplorCommand=lSim.fastCommand

        from xplorPot import XplorPot
        potCombinations = [["bond"],
                           ["bond","angl"],
                           ["bond","impr"],
                           ["bond","angl","impr"]]
        
        if useVDW:
            potCombinations.append( potCombinations[-1]+['VDW'] )
            pass
        
        for term in extraTerms:
            potCombinations.append( potCombinations[-1]+[term] )
            if useVDW:
                potCombinations.append( potCombinations[-3]+[term] )
                pass
            pass

        pots = potCombinations[3] #first time try all terms
        action = actions[0]       # and minimization
        startCoords=lSim.atomPosArr()
        iter = 0
        while iter < maxIters:
    
            viols=covalentViols(lSim.fastCommand)
            if verbose>0: print("iter", iter, "  violations:",viols)
            numViols=reduce(lambda x,y:x+y,viols)
    
    
            if numViols<minState['numViols']:
                minState['numViols']=numViols
                minState['viols']   =viols
                minState['coords']  =lSim.atomPosArr()
                pass
    
    
            if numViols<=maxViols: break

            if oldPots: action = random.choice(actions)
            
            #
            # randomly choose which minimization scheme to apply
            #
            while 1 and oldPots:
                pots = random.choice(potCombinations)
                if pots!=oldPots: break
                pass
            oldPots=pots

            lXplorCommand("flags exclude * include %s end" % " ".join(pots))
    
            if verbose>0:
                print("\t action: %s  " % action, end=' ')
                print("potential terms:", end=' ')
                print(" ".join(pots))
                pass
    
            printInterval=0
            if verbose>2:
                printInterval=1
                pass
    
            if action=="min":
                energy=max(1.,float(lXplorCommand("energy end","ener")[0]))
                lXplorCommand("""
                              set abort=off end
                              mini powell
                                nstep=500
                                drop=%f
                                nprint=%d
                              end""" % (energy/10,printInterval))
            else:
                lXplorCommand("""
                              set abort=off end
                              dynamics verlet
                                nstep=500
                                timestep=0.001
                                nprint=%d
                                tbath=3000
                                iasvel=Maxwell
                                ntrfrq=10
                              end""" % printInterval)  
                pass
            iter += 1
            #need to check unknown atoms from XPLOR side because
            #syncing back is not happening here.
            unknown = lXplorCommand("vector identity ( store1 ) ( not known )",
                                    "select")[0]
            if int(unknown)>0:
                if verbose>0: print("\t Unknown coordinates. Resetting.")
                iter=0
                lSim.setAtomPosArr( startCoords )
                lSim.syncTo() # need this because we're using fastCommand above
                pass
            pass
        lSim.syncFrom() # need this because we're using fastCommand above
        return
    
    def restoreContext():
        if verbose<=3:
            xplorCommand("set print=%s end"%prevPrintFile)
            pass
        sim.setAtomMassArr(savedMass)
        oldCurrentSimulation.sync()
        simulation.makeCurrent( oldCurrentSimulation )
        xplorCommand("constr inter (all) (all) end")
        xplorCommand("set toltorsion %f end"%prevTorsionTol)
        xSim.enableOutput(outputState)
        pass

    def printViolsAndRaiseException(xSim,minState):
        if verbose<=3:
            xplorCommand("set print=%s end"%prevPrintFile)
            pass

        viols = covalentViols(xSim.fastCommand)
        mString= "  (violations:  bond: %d  angle: %d  improper: %d" % \
              tuple(viols[:3])
        for term in extraTerms:
            print("Violated restraints in term: ",term.instanceName())
            if hasattr(term,"restraints"):
                for r in term.restraints():
                    if r.violated():
                        print("%40s %10f"% (r.name(), r.diff()))
                        pass
                    pass
                pass
            mString += "  %s: %d" % (term.instanceName(), term.violations())
            pass
        mString += ")\n"
        
        print(mString)

        mess = "Covalent geometry still violated after fixupCovalentGeom\n"
        mess += mString
        
        if not useVDW:
            mess += "Try increasing maxIters or enabling the useVdw flag."
        else:
            mess += "Try increasing maxIters."
            pass

        restoreContext()
        if not suppressExceptions:
            raise CovalentViolation(mess)
        return

    # increment through all residues, calling fixupLoop on each.
    # other residues can be fixed
    # or grouped appropriately

    from atomSel import intersection, notSelection
    for atom in intersection(sel,AtomSel("tag")):
        segid = atom.segmentName()
        resid = atom.residueNum()

        import simulation
        mainSimulation = simulation.currentSimulation()
        mainActiveAtoms = intersection(sel,
                                       AtomSel('segid "%s" and resid %d' %
                                               (segid,resid)) )
        if len(mainActiveAtoms)==0:
            continue
        
        from xplorSimulation import XplorSimulation
        rSim = XplorSimulation(1)
        #rSim.fastCommand("set mess on echo on print OUTPUT end")
        rSim.fastCommand("constraints inter (all) (all) end")
        rSim.fastCommand("constraints fix (attr mass=0) end")
        # the next line to suppress SCRATC-warning messages
        rSim.fastCommand("restraints dihe reset end end")
        rSim.deleteAtoms( notSelection(mainActiveAtoms),
                          noSync=True                   )
        rSim.syncFrom()
        simulation.makeCurrent(rSim)
        activeAtoms = AtomSel("all")

        ivm=IVM()
        ivm.setVerbose(0) #don't print info messages
        cartesianTopology(ivm)

        if verbose:
            print("processing segid:%s resid:%.4d" % (segid, resid))
            pass
        
        minState['numViols']=1e30 #reset minstate update in fixupLoop
        fixupLoop_xplor(ivm,2*maxIters,verbose-1,0,minState)

        simulation.makeCurrent(xSim)

        if minState['numViols']>0 and verbose:
            print("Error regularizing resid: %d segid: %s" % (resid,segid))
            #   printViolsAndRaiseException(rSim,minState)
            pass
            

        for i in range(len(activeAtoms)):
            mainActiveAtoms[i].setPos( activeAtoms[i].pos() )
            pass
        
        pass

    ivm=IVM()
    ivm.setVerbose(0) #don't print info messages
    ivm.fix("attr mass = 0")
    cartesianTopology(ivm,sel=sel)
    
    #now fixup everything together
    xplorCommand("""constraints
    interaction (attr mass>0) (attr mass=0)
    interaction (attr mass>0) (attr mass>0)
    end""")

    minState['numViols']=1e30 #enable minstate update in fixupLoop
    minState = fixupLoop(ivm,maxIters,verbose,maxViols,minState)

    pass

    if minState['numViols']>maxViols:
        print("fixupCovalentGeom: Covalent geometry still violated at exit.")
        print("using best structure:")
        xSim.setAtomPosArr( minState['coords'] )
        printViolsAndRaiseException(xSim,minState)
        pass

    restoreContext()
    return minState

def fixupCovalentGeomIVM(sel="known",
                         rigidRegions=(),
                         translateRegions=(),
                         useVDW=0,
                         useDynamics=1,
                         dynRatio=5,
                         maxIters=40,
                         verbose=0,
                         maxViols  =0,
                         bondTol   =0.01,
                         angleTol  =2.,
                         torsionTol=2.,
                         extraTerms=[],
                         suppressExceptions=False
                         ):
    """given the atoms in selection sel, perform, minimization and (optionally)
    dynamics to remove bond, angle, and improper violations - so that there
    total number is less than or equal to maxViols.

    rigidRegions is a sequence of selections which specify those regions to
    move as a rigid unit.
    
    translateRegions is a sequence of selections which specify those regions to
    move as a rigid unit, allowing only translational motions (i.e. without
    rotation).
    
    If useVDW is set, the nonbonded term will be used 1/4 of the time.
    if useDynamics is set, dynamics will be used 1/dynRatio of the time.
    maxIters is the total maximum number of iterations.

    tolerances for bond, bond angle and improper torsion angles can also be
    specified using optional arguments.

    Additional terms to satify may be specified with the extraTerms
    argument. This should be a list of <m pot>.Pot terms.

    If this function is not successful in satisfying all covalent restraints,
    it will throw the exception regularize.CovalentViolation, if
    suppressExceptions is False. The structure is set to that with the
    minimal number of violations, and if more than one configuration
    has that number of violations, the structure with the minimal
    energy is used. 

    This function resets the XPLOR constraints interaction settings.

    If verbose>0, intermediate status is printed. Increase verbose for more
    verbosity.

    This function is not <m ensembleSimulation>-thread safe, meaning that
    that if it is called for one ensemble member, it must be called for all
    ensemble members. 
    """
    from selectTools import convertToAtomSel
    sel = convertToAtomSel(sel)

    #leave pseudoatoms alone
    from atomSel import AtomSel, intersection
    sel = intersection(sel, AtomSel("not pseudo",sel.simulation()))

    import simulation
    oldCurrentSimulation = simulation.currentSimulation()
    sim = sel.simulation()

    from xplorSimulation import getXplorSimulation    
    xSim = getXplorSimulation(sim)
    simulation.makeCurrent( xSim )

    xplorCommand = xSim.fastCommand
    outputState=xSim.disableOutput()
    prevPrintFile = xplorCommand("set print=off end","prev_print_file")[0]

    savedMass=[]
    indices=sel.indices()

    savedMass=sim.atomMassArr()

    from atomSel import notSelection
    
    massCnt=0
    for atom in notSelection(sel):
        atom.setMass(massCnt)
        pass

    #convert rigidRegions (and translateRegions) to list of AtomSels
    rigidRegions=list(rigidRegions)
    for i in range(len(rigidRegions)):
        rigidRegions[i]=convertToAtomSel(rigidRegions[i])

    translateRegions=list(translateRegions)
    for i in range(len(translateRegions)):
        translateRegions[i]=convertToAtomSel(translateRegions[i])
        pass

    rigidRegions += translateRegions

    for rsel in rigidRegions:
        massCnt -= 1
        for atom in rsel:
            atom.setMass(massCnt)
            pass
        pass
    
        

    def covalentViols():

        ret = [int(xplorCommand("print threshold %f %s end"%(name_thresh1[1],name_thresh1[0]),
                                   "violations")[0]) for name_thresh1 in (("bonds",bondTol),
                   ("angles",angleTol),
                   ("impropers",torsionTol))]
        for term in extraTerms:
            ret.append( max(0,term.violations()) )
        return ret

    
    

    from potList import PotList
    from xplorPot import XplorPot
    allTerms=PotList()
    for name in ["bond","angl","impr"]:
        allTerms.append(XplorPot(name))
        pass
    
    potCombinations = [[XplorPot(name) for name in list] for list in [["bond"],
                           ["bond","angl"],
                           ["bond","impr"],
                           ["bond","angl","impr"]]]

    if useVDW:
        from protocol import initNBond
        initNBond(nbxmod=2,cutnb=6.5,tolerance=2.0,rcon=0.01)
        potCombinations.append( potCombinations[-1]+[XplorPot('VDW')] )
        pass

    for i in range(len(potCombinations)):
        for term in extraTerms:
            potCombinations.append(potCombinations[i]+[term])
            pass
        pass

    consInterCmd="""constraints
                      interaction (attr mass>0) (attr mass>0)
                      interaction (attr mass=0) (attr mass>0)
                      """
    for i in range(len(rigidRegions)):
        mi = -1 - i
        consInterCmd += "interaction (attr mass>0) (attr mass=%d)\n"% mi
        consInterCmd += "interaction (attr mass=0) (attr mass=%d)\n"% mi
        for j in range(i+1,len(rigidRegions)):
            mj = -1 - j
            consInterCmd += "inter (attr mass=%d) (attr mass=%d)\n" % (mi,mj)
            pass
        pass
    consInterCmd+="end"

    numTranslate=0
    for aSel in translateRegions:
        numTranslate += len(aSel)
        pass

    print('num fixed:', len(AtomSel("attr mass=0")))
    print('num grouped:', len(AtomSel("attr mass<0")))
    print('  [of these, %d atoms are in groups which can only translate]' % \
          numTranslate)
    print('num free:', len(AtomSel("attr mass>0")))

    #print consInterCmd
    
    xplorCommand(consInterCmd)


    import random

    minNumViols=1e30
    minEnergy=1e30

    from ivm import IVM, IVMError
    ivm=IVM()
    ivm.setVerbose(0) #don't print info messages
    ivm.fix("attr mass = 0")
    for aSel in rigidRegions:
        ivm.group(aSel)
        pass
    for aSel in translateRegions:
        ivm.hinge("translate",aSel)
        pass
    from protocol import cartesianTopology, initMinimize, initDynamics
    cartesianTopology(ivm,sel=sel)


    from atomAction import SetProperty
    actions = ["min"]*(dynRatio-1)
    if useDynamics:
        AtomSel("all",sel.simulation()).apply(SetProperty("mass",100.))
        for region in rigidRegions:
            region.apply(SetProperty("mass",100./(1+len(region))))
            pass
        sel.apply(SetProperty("fric",10.))
        actions.append("dyn")
        pass


    oldPots=None
    for iter in range(0,maxIters):

        viols=covalentViols()
        energy=allTerms.calcEnergy()
        if verbose: print("iter", iter, "  violations:",viols)
        numViols=reduce(lambda x,y:x+y,viols)

        if numViols<=maxViols: break

        if numViols<minNumViols or (numViols==minNumViols and
                                    energy<minEnergy):
            minNumViols=numViols
            minEnergy=energy
            minViols=viols
            minCoords=xSim.atomPosArr()
            pass


	#
	# randomly choose which minimization scheme to apply
	#
        while 1:
            pots = random.choice(potCombinations)
            if pots!=oldPots: break
            pass
        oldPots=pots

        potList = PotList()
        for pot in pots:
            potList.add( pot )
            pass

        action = random.choice(actions)
        if verbose:
            print("\t action: %s  " % action, end=' ')
            print("potential terms:", end=' ')
            print([p.instanceName() for p in potList])
            pass

        printInterval=0
        if verbose>1:
            printInterval=1
            pass

        if verbose>2:
            ivm.setVerbose(                 
                  ivm.printResetCM         |
                  ivm.printVelFromCartCost |
                  ivm.printTemperature     |
                  ivm.printEnergy          |
                  ivm.printCMVel           |
                  ivm.printStepDebug       |
                  ivm.printStepInfo        |
                  ivm.printNodeDef         |
                  ivm.printLoopDebug       |
                  ivm.printLoopInfo )
            pass
        if verbose>3:
            ivm.setVerbose(ivm.verbose() | ivm.printNodePos)
            pass
        if verbose>4:
            ivm.setVerbose(ivm.verbose() | ivm.printNodeTheta)
            pass
        if action=="min":
            energy=abs(potList.calcEnergy())
            initMinimize(ivm,numSteps=500,printInterval=printInterval,
                         dEPred=energy/10.,
                         potList=potList)
        else:
            initDynamics(ivm,numSteps=500,stepsize=0.001,
                         printInterval=printInterval,
                         bathTemp=1000,initVelocities=True,
                         potList=potList)
            pass
        try:
            ivm.run()
        except IVMError as e:
            print("Encountered ivm error:", e.message)
            pass
        pass

    xplorCommand("set print=%s end"%prevPrintFile)
    xSim.enableOutput(outputState)

    def restoreContext():
        sim.setAtomMassArr(savedMass)
        oldCurrentSimulation.sync()
        xplorCommand("constr inter (all) (all) end")
        simulation.makeCurrent( oldCurrentSimulation )
        pass

    if numViols>maxViols:

        print("fixupCovalentGeomIVM: Covalent geometry still violated at exit.")
        print("using best structure:")

        xSim.setAtomPosArr( minCoords )
        viols = covalentViols()
        mString= "  (violations:  bond: %d  angle: %d  improper: %d" % \
              tuple(viols[:3])
        for term in extraTerms:
            print("Violated restraints in term: ",term.instanceName())
            if hasattr(term,"restraints"):
                for r in term.restraints():
                    if r.violated():
                        print("%40s %10f"% (r.name(), r.diff()))
                        pass
                    pass
                pass
            mString += "  %s: %d" % (term.instanceName(), term.violations())
            pass
        mString += ")\n"
        
        print(mString)

        mess = "Covalent geometry still violated after fixupCovalentGeom\n"
        mess += mString
        
        if not useVDW:
            mess += "Try increasing maxIters or enabling the useVdw flag."
        else:
            mess += "Try increasing maxIters."
            pass

        restoreContext()
        if suppressExceptions:
            return
        else:
            raise CovalentViolation(mess)
        pass

    restoreContext()
    return

class CovalentViolation(Exception):
    def __init__(s,mess):
        Exception.__init__(s,mess)
    pass

def addUnknownAtoms_fast(verbose=0,
                         maxFixupIters=500,
                         simulation=0):
    """add in unknown atoms so that covalent and vdw terms are satisfied by
    randomly placing atoms, and calling <m protocol>.fixupCovalentGeom() to
    correct the covalent geometry, using maxFixupIters as the maxIters
    argument to that function.

    This is a fast version which will only reliably work if the missing atoms
    are separated by a single bond from a known atom.

    if verbose=True, details of the minimization procedure are printed.

    """
    from xplorSimulation import getXplorSimulation
    xSim = getXplorSimulation(simulation)

    from atomSel import AtomSel
    from random import uniform
    from vec3 import Vec3, unitVec
    unknownSel = AtomSel("not known",xSim)
    for atom in unknownSel:
        try:
            boundAtom = AtomSel("known and bondedto id %d" %
                                (atom.index()+1),xSim)[0]
            boundAtomPos = boundAtom.pos()
        except IndexError:
            #isolated atoms (with no bonds) 
            min=-5
            max=5
            boundAtomPos = Vec3(uniform(min,max),uniform(min,max),
                                uniform(min,max))
            pass
        dir = Vec3(uniform(0,1),uniform(0,1),uniform(0,1))
        atom.setPos( boundAtomPos + unitVec(dir))
        pass

    if verbose:
        print('added %d atoms.' % len(unknownSel))
        print('now fixing covalent geometry...')
        pass
    try:
        fixupCovalentGeom(useVDW=1,maxIters=maxFixupIters,sel=unknownSel,
                          dynRatio=10, verbose=verbose)
    except CovalentViolation as x:
        print('addUnknownAtoms_fast: covalent exceptions:')
        print(str(x))
        pass
        
    pass

def addUnknownAtoms(dyn_stepsize=0.02,
                    dyn_numStepMul=1,
                    maxFixupIters=20,
                    maxFastIterations=100,
                    verbose=0,
                    simulation=0):
    """Iteratively delete unknown atoms bound to unknown atoms and call
    addUnknownAtoms_fast.

    This algorithm will fail for loops of unknown atoms- and the slow
    addUnknownAtoms_old routine will be called.
    """
    thisSimulation=simulation
    import simulation
    oldCurrentSimulation = simulation.currentSimulation()
    from xplorSimulation import XplorSimulation, getXplorSimulation

    xsim = getXplorSimulation(thisSimulation)

    cnt=0
    from atomSel import AtomSel
    unknownSel=AtomSel("not known",xsim)
    while len(unknownSel)>0:
        cnt += 1
        if verbose:
            print("iteration: %d, unknown atoms: %d"% (cnt, len(unknownSel)))
            pass
        simulation.makeCurrent(xsim)

        subSim=XplorSimulation(True) # a cloned simulation

        subSim.fastCommand("delete sele=(%s) end" %
                           "not bondedto known and not known")
        subSim.syncFrom()

        unknown=AtomSel("not known",subSim)

        if len(unknown)==0:
            if verbose:
                print("Encountered a loop.")
                print("  Falling back to the slow addUnknownAtoms_old path.")
                pass
            addUnknownAtoms_old(dyn_stepsize=dyn_stepsize,
                                dyn_numStepMul=dyn_numStepMul,
                                verbose=verbose,
                                maxFixupIters=maxFixupIters)
        else:
            simulation.makeCurrent(subSim)
            addUnknownAtoms_fast(verbose,
                                 maxFixupIters,
                                 subSim)

            for atom in unknown:
                AtomSel('atom "%s" %d %s' %(atom.segmentName(),
                                            atom.residueNum(),
                                            atom.atomName()),
                        xsim)[0].setPos(atom.pos())
                pass
            pass
        
        unknownSel.reevaluate()
        if cnt>=maxFastIterations:
            print("Warning: addUnknownAtoms: too many calls to addUnknownAtoms_fast.")
            break
        pass


    simulation.makeCurrent(oldCurrentSimulation)
    return
#alias for backwards compatibility
addUnknownAtoms_new=addUnknownAtoms
    


def addUnknownAtoms_old(dyn_stepsize=0.02,
                        dyn_numStepMul=1,
                        maxFixupIters=500,
                        verbose=0,
                        simulation=0):
    """Add unknown atoms, satisfying covalent and nonbonded terms.

    This routine is slow, but it is rather robust.

    dyn_stepsize specifies the timestep size during the MD phase.  Reduce this
    if you have convergence problems.

    dyn_numStepMul is a multiplier for the number of molecular dynamics steps
    taken.  Increase this to get better convergence.

    if verbose=True, details of the minimization procedure are printed.

    This function resets the XPLOR constraints interaction settings.

    """

    from xplorSimulation import getXplorSimulation
    xSim = getXplorSimulation(simulation)

    from atomSel import AtomSel
    unknownSel=AtomSel("not known")
    if len(unknownSel)==0:
        return

    if not len( AtomSel("bondedto not known and not known",xSim) ):
        if verbose:
            print("using the fast addUnknownAtoms protocol.")
            pass
        return addUnknownAtoms_fast(verbose,simulation=simulation,
                                    maxFixupIters=maxFixupIters)
    

    outputState=xSim.disableOutput()
    if not verbose: xSim.command("set print=off end")

    dyn_numStep = 500 * dyn_numStepMul
    dyn_ramp_numStep = 100 * dyn_numStepMul

    xSim.command("    evaluate ($timestep = %f)" % dyn_stepsize)
    xSim.command("    evaluate ($ramp_nstep = %f)" % dyn_ramp_numStep)
    xSim.command("    evaluate ($nstep = %f)" % dyn_numStep)


    xSim.command(r"""

vector do (fbeta=10) (not known)    {*Friction coefficient for MD heatbath.*}
vector do (q=mass)   (all)      {* Save the real masses for later *}
vector do (mass=100) (not known)            {*Uniform heavy masses to speed*}
vector do (fbeta=0) (known)    {*Friction coefficient for MD heatbath.*}
vector do (mass=0)  (known)              {*Uniform heavy masses to speed*}

    eval ($knoe=0.1)
constraints fix (known) end
vector do (vx = 0) (known)
vector do (vy = 0) (known)
vector do (vz = 0) (known)

evaluate ($init_t = 1000 )    {* Initial simulated annealing temperature.*}

vector do (vx = maxwell($init_t)) (attr mass>0)
vector do (vy = maxwell($init_t)) (attr mass>0)
vector do (vz = maxwell($init_t)) (attr mass>0)
vector do (x=(random()-0.5)*20) (attr mass>0)
vector do (y=(random()-0.5)*20) (attr mass>0)
vector do (z=(random()-0.5)*20) (attr mass>0)

!try bonds first
flags exclude * include bond  end

constraints 
 interaction (attr mass>0) (attr mass=0) 
 interaction (attr mass>0) (attr mass>0)
end

!mini powell
!  drop=1e5
!  nprint=1
!  tolg=1e-5
!  nstep=1000
!end

flags exclude * include bond angle dihe cdihe impr vdw  end

!energy end
!mini powell
!  drop=1e5
!  nprint=1
!  tolg=1e-5
!  nstep=1000
!end

    evaluate ($kbon = 0.00005  )                 {* Bonds.                 *}
    evaluate ($kang = 0.00005  )                 {* Angles.                *}

    !constraints 
    !	 interaction (all) (all) 
    !	 weights bond $kbon angl $kang impr $kimp vdw 0 elec 0 end 
    !end
    !dynamics  verlet
    !	   nstep=5000   timestep=$timestep   iasvel=current   
    !	   tcoupling=true tbath=$init_t  nprint=50  iprfrq=0 
    !end	    

    while ($kbon < 0.01) loop stage1
      
            evaluate ($kbon = min(0.25, $kbon * 1.25))
            evaluate ($kang = $kbon)
            evaluate ($kimp = $kbon/10)
	   
            noe scale * $knoe end
            !restraints dihed scale 0. end 
            constraints 
                interaction (attr mass>0) (attr mass=0) 
                weights bond $kbon angl $kang impr $kimp vdw 5e-4 elec 0 end 
                interaction (attr mass>0) (attr mass>0) 
                weights bond $kbon angl $kang impr $kimp vdw 5e-4 elec 0 end 
            end
	    	    
            dynamics  verlet
                  nstep=$ramp_nstep   timestep=$timestep   iasvel=current   
                  tcoupling=true tbath=$init_t  nprint=50  iprfrq=0 
            end	    
	    
    end loop stage1

    parameter                {* Parameters for the repulsive energy term.  *}
      nbonds
        repel=0.9            {* Initial value for repel - modified later.  *}
        nbxmod=-3            {* Initial value for nbxmod - modified later. *}
        wmin=0.01 
        cutnb=4.5 ctonnb=2.99 ctofnb=3. 
        tolerance=0.5 
      end
    end
    ! add vdw and slowly increase its weight
    parameter nbonds 
         atom cutnb 100 tolerance 45 repel=1.2 
         rexp=2 irexp=2 rcon=1.0 nbxmod 4
    end end
    flags exclude * include bond angle impr dihe cdihe vdw end

    constraints 
        interaction (attr mass>0) (attr mass=0) 
        weights bond $kbon angl $kang impr $kimp vdw 0.002 elec 0 end 
        interaction (attr mass>0) (attr mass>0) 
        weights bond $kbon angl $kang impr $kimp vdw 0.002 elec 0 end 
    end
    dynamics  verlet
          nstep=$nstep   timestep=$timestep   iasvel=current   
          tcoupling=true tbath=$init_t  nprint=50  iprfrq=0 
    end	    
    constraints 
        interaction (attr mass>0) (attr mass=0) 
        weights bond $kbon angl $kang impr $kimp vdw 0.005 elec 0 end 
        interaction (attr mass>0) (attr mass>0) 
        weights bond $kbon angl $kang impr $kimp vdw 0.005 elec 0 end 
    end
    dynamics  verlet
          nstep=$nstep   timestep=$timestep   iasvel=current   
          tcoupling=true tbath=$init_t  nprint=50  iprfrq=0 
    end	    
    constraints 
        interaction (attr mass>0) (attr mass=0) 
        weights bond $kbon angl $kang impr $kimp vdw 0.01 elec 0 end 
        interaction (attr mass>0) (attr mass>0) 
        weights bond $kbon angl $kang impr $kimp vdw 0.01 elec 0 end 
    end
    dynamics  verlet
          nstep=$nstep   timestep=$timestep   iasvel=current   
          tcoupling=true tbath=$init_t  nprint=50  iprfrq=0 
    end	    
    constraints 
        interaction (attr mass>0) (attr mass=0) 
        weights bond $kbon angl $kang impr $kimp vdw 0.5 elec 0 end 
        interaction (attr mass>0) (attr mass>0) 
        weights bond $kbon angl $kang impr $kimp vdw 0.5 elec 0 end 
    end
    dynamics  verlet
          nstep=$nstep   timestep=$timestep   iasvel=current   
          tcoupling=true tbath=$init_t  nprint=50  iprfrq=0 
    end	    

!flags exclude * include bond angle end

energy end
mini powell
  drop=100
  nprint=1
  tolg=1e-5
  nstep=10000
end

    constraints 
        interaction (attr mass>0) (attr mass=0) 
        weights bond $kbon angl $kang impr $kimp vdw 0.1 elec 0 end 
        interaction (attr mass>0) (attr mass>0) 
        weights bond $kbon angl $kang impr $kimp vdw 0.1 elec 0 end 
    end

mini powell
  drop=100
  nprint=1
  tolg=1e-5
  nstep=10000  
end

vector do (mass=q) (all) {* Return the masses to sane values *}

""")
    xSim.command("constraints inter (all) (all) end")
    xSim.command("set print=OUTPUT end")
    xSim.enableOutput(outputState)

    try:
        fixupCovalentGeom(useVDW=1,maxIters=maxFixupIters,sel=unknownSel,
                          dynRatio=10, verbose=verbose)
    except CovalentViolation as x:
        print('addUnknownAtoms: covalent exceptions:')
        print(str(x))
        pass
    return

