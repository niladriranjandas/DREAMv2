"""high-level refinement tools
"""

from sys import modules
from potList import PotList
from pdbTool import PDBTool
from functools import reduce


class StructureLoop:
    """ class which performs loop over structure calculations.
    Constructor: StructureLoop()
      arguments:
        numStructures       -    number of structures to calculate
        startStructure      -    id of the first structure to calculate
        structureNums       -    sequence of explicit structure numbers to
                                 calculate.
        pdbFilesIn          - A list of existing coordinate files, or a
                              string including a <s glob> wildcard. If
                              specified, structLoopAction will be called
                              after reading each specified file. In this
                              case, if the arguments numStructures or
                              structureNums are specified, the coordinates
                              will be initialized to the appropriate structure
                              number modulus the length of pdbFilesIn. If
                              numStructures and structureNums are not
                              specified, len(pdbFilesIn) structures will be
                              calculated. This argument can also be a template
                              string including the STRUCTURE literal which is
                              replaced by structure number. For ensemble
                              simulations, a MEMBER literal should also be
                              present, if separate coordinates for
                              are desired for each ensemble member. Coordinates
                              are initialized by calling <m protocol>.initCoords,
                              the return value of which is stored in
                              the initCoordsRet member of the instance object
                              (passed as the argument to the structLoopAction
                              function).
        pdbFilesInErase     - the value to pass to the erase argument of
                              <m protocol>.initCoords, when pdbFilesIn is
                              specified.
                              [default: False]
        inMemberModulo      - is used in conjunction with pdbFilesIn in
                              the case that an ensemble calculation is
                              using input structures from a smaller ensemble.
                              This number specifies the size of the smaller
                              input ensemble.
        structLoopAction    - a user-defined function which takes
                              one argument: an instance of this class.
                              If this argument is omitted, new structures
                              are not calculated. Rather, existing structures
                              in files specified by pdbFilesIn are read-in,
                              and analyzed. There is a special null class
                              member named sharedData created during run()
                              which can be written to in structLoopAction
                              and read back out after run() is called,
                              but at that point, structData is a
                              dictionary whose keys are structure numbers. The
                              argument to the function is an instance of this
                              class, such that structure-specific info can be
                              accessed. Useful members include:

                                structNum    - the current structure number .
                                filename()   - returns the output filename
                                               for the current structure.
                                extraRemarks - extra data to place in the
                                               REMARKs section of the output
                                               PDB.
                                               
        pdbTemplate         - template string used to create a filename
                              for the pdbFile method. The filename is
                              generated using the makeFilename method.
                              Note that the template string should always
                              include the STRUCTURE literal so that distinct
                              structure files are generated. For
                              <m ensembleSimulation> calculations a MEMBER
                              literal should also be present *after* the
                              STRUCTURE literal.
                              [default: "SCRIPT_STRUCTURE.pdb" for
                              non-EnsembleSimulations and
                              "SCRIPT_STRUCTURE_MEMBER.pdb" for
                              EnsembleSimulations]
        doWriteStructures   -    write the structures using writeStructure
                                 after each call to structLoopAction. The
                                 potList argument is used as the PotList for
                                 analysis, the member altPotList (which
                                 defaults to the value of averageCrossTerms)
                                 is used as the second argument, and the
                                 member extraRemarks is used as the third
                                 argument.
                                 [default: False]
        writeCrashStructures -   If structLoopAction (the structure calculation
                                 function) throws an exception, a crash
                                 coordinate file will be created if this is
                                 True.
                                 [default: True]
        writeSelection       -   Atom selection specifying which atoms to write
                                 to output files, whether those written from
                                 the pdbTemplate, crash structures, or the
                                 averaged coordinates.
                                 [default: all]
        potList              -   The potList argument to writeStructures if
                                 doWriteStructures is True. If it is not set,
                                 it defaults to averagePotList.
        ivm                  -   An <m ivm>.IVM object, used for consistency
                                 checking. If set, a check is made if fitting
                                 of atomic coordinates is attempted while
                                 (non-pseudo) atoms are fixed in the IVM. If so,
                                 an exception is raised.
        calcMissingStructs   -   If True, structures are only calculated if
                                 their associated filename does not exist.
                                 However, all structures specified by
                                 numStructures or structureNums are read in
                                 and used for analysis.
                                 [default: False]
                                 
      if numStructures<0,   existing files matching pdbTemplate are processed,
      starting at startStructure, and stopping when a file does not exist.
      This mode of operation does not work with structure parallelism, but
      it does work with ensemble parallelism.

     There are additional arguments if you would like an average structure
     to be calculated. If averaging is enabled, the output structure files
     will be fit by averageFitSel to that which has the lowest averageSortPots
     energy. Ensemble Structures from <m ensembleSimulation>.EnsembleSimulation
     calculations structures are fit to the respective ensemble member
     of the ensemble which has the lowest averageSortPots energy. This
     approach to ensemble fitting is often not appropriate.
     
       averageFilename      - name for output structure file. If not specified,
                              the average structure will not be calculated.
       averagePotList       - potential terms to use for average structure
                              calculation. These terms are reported on in the
                              .stats file.
       averageRegularize    - flag determining whether or not structure
                              regularization (gradient minimization) is
                              carried out on the average structure, by
                              minimizing against averagePotList.
                              [default: True, if averageFilename is set]
       averageFixedRegions  - sequence of regions held fixed in space during
                              structure regularization.
       averageRigidRegions  - sequence of regions held rigid during structure
                              regularization.
       averageSortPots      - potential terms used for sorting structures. The
                              top fraction or number of structures is reported
                              on in the .stats file.
                              [defaults to averagePotList]
       averageCrossTerms    - potential terms to report on, but not to use
                              in refinement of average structure.
       averageContext       - function to call to initialize force constants
                              and other parameters for the energy/observable
                              summary printed at the beginning of the output
                              structure files, .viols files, and statistical
                              restraint analysis. It is also used to setup the
                              average minimization run.  
       averageFitSel        - atom selection used to fit structures in
                              calculation of average structure and to fit
                              structures to the lowest energy
                              structure. It defaults to "name CA" for non-ensemble
                              calculations. Setting this to None will disable
                              this fitting.
       averageCompSel       - atom selection used for heavy atom rmsd
                              comparison metric [not name H* and not PSEUDO].
       averageTopFraction   - fraction of structures to use in calculation of
                              average structure. The structures are sorted
                              by energy, and the top fraction is retained.
                              [1 - all structures].
       averageTopNum        - number of structures to use in calculation of
                              average structure. Specify only one of 
                              averageTopFraction or averageTopNum.

       averageAccept        - function to call to assess whether it is
                              acceptable, i.e. meets violation, rmsd
                              requirements. The function has a single argument:
                              averagePotList [defaults to accepting all
                              structures.]

       averageRefineSteps   - number of minimization steps to take during
                              regularization refinement of the average
                              structure [50]

       genViolationStats    - flag controlling whether statics are gathered
                              (over all structures) on which restraints are
                              violated most often. The results are collected
                              in a file named pdbTemplate.stats. Statistics
                              will be gathered for all terms in
                              averagePotList, so this attribute must be
                              specified for this facility to work.
       storeCoordinates     - if True, atomic coordinates are saved in memory
                              for final analysis, but this may use excessive
                              memory if many structures are calculated. By
                              default, storeCoordinates=False.

       averageRestrain      - flag the control the inclusion of probDist energy
                              potential in calculating the average structure
                              from the ensemble. averageRestrainSel helps in
                              selecting atoms from which the density map is
                              created. inconsistentAveStruct is a flag that
                              is set to True if the sortEnergy (and associated
                              violation count) of the calculated average
                              structure is greater than all of the structures
                              used to compute it.
                              

    method: run():
       performs numStructures loop of action. In each pass, the coordinates
       of the current Simulation are reset to their initial values and the
       instance variable count is incremented. Also, the global random seed
       is incremented. If the current simulation is an EnsembleSimulation,
       the seed-setting takes this into account.

       After run() has completed, average structure coordinates are left in
       the current Simulation, if they have been calculated. If restraint
       statistics are generated, the StructureLoop instance will have the
       following members when run() returns:
         restraintStats 
         restraintStatsCross
       These are <m restraintStats>.RestraintStats objects corresponding to
       potential terms in averagePotList and AverageCrossTerms, respectively.
       The precision of the calculated structures is stored in the members
         fitRMSD
         compRMSD
       corresponding to averageFitSel and averageCompSel, respectively. The
       member
         structInfo
       contains a list of StructInfo objects, one for each structure reported
       on the the .stats file.
       Also, the cpu time spent within the run() method will be contained
       in the members
         cpuTime
         cpuTimes
         cpuTimeTot
       The first contains the local process's cpu time, the second is an array
       of times from each process, and the third is the sum.
    """
    def __init__(s,numStructures=-1,startStructure=0,
                 structureNums=[],
                 pdbFilesIn=[],
                 pdbFilesInErase=False,
                 inMemberModulo=None,
                 structLoopAction="",
                 pdbTemplate=None,
                 doWriteStructures=False,
                 calcMissingStructs=False,
                 writeCrashStructures=True,
                 writeSelection="all",
                 potList=None,
                 ivm=None,
                 genViolationStats=False,
                 storeCoordinates=False,
                 averageFilename="",
                 averagePotList=PotList(),
                 averageRegularize=-1,
                 averageFixedRegions=[],
                 averageRigidRegions=[],
                 averageSortPots=None,
                 averageCrossTerms=[],
                 averageContext=lambda : 1,
                 averageFitSel=-1,
                 averageCompSel="not name H* and not PSEUDO",
                 averageTopFraction=1,
                 averageTopNum=    -1,
                 averageAccept=lambda potList: 1,
                 averageRefineSteps=50,
                 averageRestrain=False,
                 averageRestrainSel="name CA or name C or name N or name O"):
        import sys
        from simulation import Simulation_currentSimulation
        import xplor

        #find ensembleSimulation, if present, and set skip member
        sim = Simulation_currentSimulation()
        s.esim = 0
        s.initCoords = sim.atomPosArr()
        s.skip=0
        memberIndex=0
        s.structNum=-1
        if ( sim.type() == "EnsembleSimulation"):
            from ensembleSimulation import EnsembleSimulation_currentSimulation
            s.esim = EnsembleSimulation_currentSimulation()
            memberIndex = s.esim.member().memberIndex()
            s.skip = memberIndex
            s.inMemberModulo = inMemberModulo if inMemberModulo else \
                               s.esim.size()+1
            memberIndex %= s.inMemberModulo
            pass
        


        if numStructures>0 and len(structureNums):
            raise Exception("specify only one of numStructures or " +
                            "structureNums")
        if len(structureNums):
            numStructures = len(structureNums)
#            if type(pdbFilesIn)==type("string"):
#                raise Exception("if structureNums is specified, then "
#                                "pdbFilesIn must be a sequence")
        elif len(pdbFilesIn):
            if type(pdbFilesIn)==type("string"):
                #assume pbdFilesIn is a template, otherwise its a list of filenames
                import glob
                template = pdbFilesIn.replace("STRUCTURE","*")
                template = template.replace("MEMBER",
                                            str(memberIndex))
                zeroFiles = glob.glob(template)
                if len(zeroFiles)==0:
                    raise Exception("template specified in pbdFilesIn matches"
                                    " no files")
                prefIndex=pdbFilesIn.find('STRUCTURE')
                # turn glob ordering into ordering based on structure number
                # -- if STRUCTURE is in template specification
                if prefIndex<0:
                    #globs are not sorted- sorting is required
                    # if a glob is specified with MEMBER
                    pdbFilesIn=sorted(zeroFiles)
                else:
                    suffIndex=prefIndex+9
                    sepString=None
                    if pdbFilesIn.find('MEMBER')>-1:
                        if pdbFilesIn.find('MEMBER')<suffIndex:
                            raise Exception("MEMBER must come after STRUCTURE "
                                            "in pdbFilesIn template")
                        suffIndex=pdbFilesIn.find('MEMBER')+6
                        sepString = pdbFilesIn[pdbFilesIn.find('STRUCTURE')+9:
                                               pdbFilesIn.find('MEMBER')]
                        pass
                    suffIndex = suffIndex - len(pdbFilesIn)
                    structsFiles=[]
                    for file in zeroFiles:
                        structNum= file[prefIndex:suffIndex]
                        if sepString:
                            structNum = structNum[:structNum.find(sepString)]
                            file=pdbFilesIn[:prefIndex]+structNum+\
                                  pdbFilesIn[prefIndex+9:]
                            pass
                        structsFiles.append( (int(structNum),file) )
                        pass
                    structsFiles.sort(key=lambda x: x[0])
                    pdbFilesIn = [t[1] for t in structsFiles]
                    pass
                pass
            #assume pdbFilesIn is a list of filenames
            if numStructures<0:
                numStructures = len(pdbFilesIn)
            pass
        s.pdbFilesInErase = pdbFilesInErase
        s.numStructures=numStructures
        s.pdbFilesIn = pdbFilesIn
        s.structLoopAction = structLoopAction
        if not pdbTemplate:
            pdbTemplate = "SCRIPT_STRUCTURE_MEMBER.pdb" if s.esim else \
                          "SCRIPT_STRUCTURE.pdb"
            pass
        s.pdbTemplate = pdbTemplate
        s.processID = xplor.p_processID
        proc = s.processID
        #if env.has_key("NUM_THREADS"):
        #    print "num_threads =", env["NUM_THREADS"]
        s.numProcs = xplor.p_numProcs

        s.genViolationStats=genViolationStats
        s.averageRestrain=averageRestrain
        s.averageRestrainSel=averageRestrainSel
        s.averageFilename =averageFilename
        s.averagePotList = convertToPotList( averagePotList )
        if not averageFilename:
            averageRegularize=False
            pass
        s.averageRegularize = averageRegularize 
        s.averageFixedRegions = averageFixedRegions
        s.averageRigidRegions = averageRigidRegions
        if not averageSortPots:
            averageSortPots = averagePotList
            pass
        s.averageSortPots = convertToPotList( averageSortPots )
        s.averageCrossTerms=convertToPotList( averageCrossTerms )
        s.averageContext   =averageContext
        from atomSel import intersection, AtomSel
        
        
        if averageFitSel==-1:
            if s.esim:
                averageFitSel=None
                pass
            else:
                averageFitSel = "name CA"
                print("averageFitSel set to:", averageFitSel)
                pass
            pass
        
        if averageFitSel and ivm and intersection(ivm.fixedAtoms(),
                                                  AtomSel("not pseudo",
                                                          ivm.simulation)):
            raise Exception("averageFitSel should be none if ivm has fixed " +
                            "non-pseudo atoms")
        s.averageFitSel    = '' if not averageFitSel else averageFitSel
        s.averageCompSel   =averageCompSel if averageCompSel else ""
        s.averageTopFraction=averageTopFraction
        s.averageTopNum     =averageTopNum
        s.averageAccept     =averageAccept
        s.averageRefineSteps=averageRefineSteps
        s.inconsistentAveStruct=0

        s.doWriteStructures=doWriteStructures
        s.writeCrashStructures = writeCrashStructures
        if storeCoordinates==None:
            if genViolationStats and not doWriteStructures:
                s.storeCoordinates = True
            else:
                s.storeCoordinates = False
                pass
            pass
        else:
            s.storeCoordinates = storeCoordinates
            pass
        if s.storeCoordinates: s.coords={}
        s.writeSelection = writeSelection
        s.potList = potList if potList!=None else s.averagePotList
        s.altPotList = s.averageCrossTerms
        s.extraRemarks = None

        s.calcMissingStructs=calcMissingStructs

        s.procStructMap=[]
        if numStructures>=0:
            # this logic used if number of jobs is greater than the desired 
            # number of structures
            s.numProcs = min(s.numProcs,numStructures)
            if proc >= s.numProcs:
                print('StructureLoop: this process has no work. Exiting...')
                xplor.p_numProcs=-1 # to avoid barrier error in xplorFinal:34
                sys.exit()
                pass

            if not structureNums:
                structureNums = list(range(startStructure,
                                      startStructure+numStructures))
                pass
            
            for i in range(s.numProcs):
                i_start = (i     * numStructures) // s.numProcs
                i_stop  = ((i+1) * numStructures) // s.numProcs
                s.procStructMap.append( structureNums[i_start:i_stop] )
#                (i     * numStructures) / s.numProcs
#                start = (i     * numStructures) / s.numProcs + startStructure
#                stop  = ((i+1) * numStructures) / s.numProcs + startStructure
#                s.procStructMap.append( (start,stop) )
                pass
            s.structNums = s.procStructMap[proc]
        else:
            # setup for processing existing files-
            if s.esim: s.skip=-2000
            #FIX: check this
            s.start=startStructure
            s.stop=-1
            pass
        
        # barrier across all parallel processes- this is required after
        # extra processes have shutdown so that further communication is
        # successful.
        from ensembleSimulation import commBarrier
        commBarrier(xplor.p_comm)
        return

    def run(s):
        import simulationWorld
        import protocol
        from simulation import Simulation_currentSimulation
        from inspect import currentframe, getouterframes
        simWorld = simulationWorld.SimulationWorld_world()
        sim = Simulation_currentSimulation()
        initCoords = sim.atomPosArr()
        s.initSeed = simWorld.random.seed()
        s.cpuTime = simWorld.cpuTime()

        s.doFitStructure = (s.genViolationStats or
                            s.averageFilename) and s.averageFitSel 
        
        cnt = 0
        s.sharedByStruct = {}
        while 1:
            s.sharedData=NullClass()
            if cnt==len(s.structNums): break
            if len(s.structNums):
                s.structNum = s.structNums[cnt]
            elif s.numStructures<0:
                s.structNum = s.start + cnt

                #logic for processing existing files- and the number
                # of structures is not specified.
                #stop if the file doesn't exist.
                try:
                    import os
                    os.stat( s.makeFilename(s.pdbTemplate) )
                except:
                    break
                pass
            # count is for backward compatibility
            s.count = s.structNum

            import os
            doCalcStructure= bool(s.structLoopAction)
            if len(initCoords) != sim.numAtoms():
                # This clause allows the number of atoms to change during
                # a structure calculation. Note that there will be bogus
                # coords for the mismatched atoms, so structLoopAction had
                # better take care of them. In particular, this is needed
                # for waterRefineTools.refine when keepWaters=True.
                print("StructureLoop.run: Warning: the number of atoms", end=' ')
                print("has changed.")
                initCoords.resize( sim.numAtoms() )
                pass
            sim.setAtomPosArr( initCoords )
            s.inputPDB = s.pdbFileIn()
            writeThisStructure = s.doWriteStructures

            structurePresent=True
            if s.calcMissingStructs: #check if structure already calculated
                if s.esim:
                    for i in range(s.esim.size()):
                        if not os.path.exists(s.filename(memberNum=i)):
                            structurePresent=False
                            pass
                        pass
                    pass
                else: #not an ensemble
                    if not os.path.exists(s.filename()):
                        structurePresent=False
                        pass
                    pass
                pass
            
            if s.calcMissingStructs and structurePresent:
                        
                if simWorld.logLevel()!='none':
                    print("StructureLoop:", end=' ')
                    print("structure %d already calculated" % s.structNum)
                    pass
                s.pdbFile().read()
                doCalcStructure=False
                pass
            elif s.inputPDB:
                if simWorld.logLevel()!='none':
                    print("StructureLoop:", end=' ')
                    print("initializing coordinates from %s" % s.inputPDB)
                    pass
                try:
                    s.initCoordsRet = \
                    protocol.initCoords(s.inputPDB , maxUnreadEntries=None,
                                        erase=s.pdbFilesInErase)
                except IOError:
                    print("Error reading structure: %s. Skipping..." % s.inputPDB)
                    doCalcStructure=False
                    writeThisStructure=False
                    pass
                pass
            
            if doCalcStructure:
                if simWorld.logLevel()!='none':
                    print("StructureLoop: calculating structure %d" % s.structNum)
                s.randomSeed = s.initSeed + s.structNum + s.skip*s.numStructures
                simWorld.setRandomSeed( s.randomSeed )
                
                try:
                    if type(s.structLoopAction) == type("string"):
                        #structLoopInfo = s
                        #exec( s.structLoopAction, vars( modules["__main__"] ), locals() )
                        global_dict = getouterframes( currentframe() )[1][0].f_globals
                        local_dict = getouterframes( currentframe() )[1][0].f_locals
                        local_dict["structLoopInfo"] = s
                        exec( s.structLoopAction, global_dict, local_dict )
                        del local_dict["structLoopInfo"]
                    else:
                        s.structLoopAction(s)
                        pass
                    pass
                except Exception as e:
                    import traceback
                    traceback.print_exc()
                    print("Error calculating structure %d: %s" % (s.structNum,
                                                                  e.args[0]))
                    if s.writeCrashStructures:
                        s.writeStructure(s.potList,s.altPotList,s.extraRemarks,
                                         nameSuffix=".crash")
                        pass
                    writeThisStructure=False
                    pass
                pass
            if writeThisStructure:
                if not s.doFitStructure or not s.averageAccept(s.averagePotList):
                    # in these cases, averageContext will not be called below,
                    # and structure will not be rewritten.
                    s.averageContext()
                    pass
                s.writeStructure(s.potList,s.altPotList,s.extraRemarks)
                s.sharedByStruct[s.structNum] = s.sharedData
                pass
            if s.storeCoordinates: s.coords[s.structNum] = sim.atomPosArr()
            if s.esim: s.esim.barrier()
            cnt += 1
            pass

        import xplor
        comm = xplor.p_comm
        if s.esim and s.esim.size()>1:
            from ensembleSimulation import Comm
            comm = Comm(comm,s.esim)
            pass
        s.sharedP = comm.collect(s.sharedByStruct)

        s.sharedData={}
        for data in s.sharedP:
            for structNum in list(data.keys()):
                s.sharedData[structNum] = data[structNum]
                pass
            pass

        if s.averageFilename or s.genViolationStats:
            s.genAveStats(comm,sim)
            pass

        s.cpuTime = simWorld.cpuTime() - s.cpuTime
        s.cpuTimeTot=0
        s.cpuTimes=[]
        s.cpuTimes=comm.collect( s.cpuTime )
        for time in s.cpuTimes:
            s.cpuTimeTot += time
            pass
        return s

    def genAveStats(s,comm,sim):
        """ generate averages, statistics for calculated structures
        """

        #1) accumulate RestraintStats
        #2)  and accumulate energy/rms/violation statistics
        from restraintStats import RestraintStats
        if s.genViolationStats:
            rStats = RestraintStats()
            rStatsCross = RestraintStats()
            pass
        structInfo = []
        #need to have the zero structure for the fit -
        # this can be a long wait - disable the timeout (wait forever).
        comm.barrier(timeout=None)  

        from atomSelAction import Fit
        from atomSel import AtomSel
        if s.numStructures<0:
            raise "genAveStats does not support numStructures<0"
        
        # now read each structure, evaluate energies, and
        # collect restraint statistics
        for s.structNum in s.structNums:
            try:
                #read files
                if s.storeCoordinates:
                    sim.setAtomPosArr( s.coords[s.structNum] )
                else:
                    s.pdbFile().read()
                    pass
                sim.sync()
                s.averageContext()
                if s.genViolationStats:
                    rStats.accumulate(s.averagePotList,s.structNum)
                    rStatsCross.accumulate(s.averageCrossTerms,s.structNum)
                    pass
                #calc energy
                try:
                    sortEnergy   = s.averageSortPots.calcEnergy()
                    sortViol=0 
                    if s.averageRestrain: 
                        sortViol     = s.averageSortPots.violations()
                        pass
                        
                except AttributeError:
                    sortEnergy = 0
                    pass

                # does the structure pass the accept() function
                accept = s.averageAccept(s.averagePotList)

                structInfo.append( StructInfo(s.pdbFile().filename(),
                                              sortEnergy,
                                              sortViol,
                                              accept,
                                              s.structNum,
                                              sim.atomPosArr() if s.storeCoordinates else None) )
            except IOError as e:
                #file is missing, and can't be read
                print("genAveStats: error reading pdb:", e.args[0])
                pass
            pass
        
        # collect all information in proc 0
        structInfo = comm.collect(structInfo)
        #convert from array of arrays to single array
        structInfo = reduce(lambda x,y:x+y, structInfo, [])
        structMap={}
        for struct in structInfo:
            structMap[struct.id] = struct
            pass

        if s.processID==0: 
            #sort structures by energy
            structInfo.sort( key=lambda x: x.sortEnergy)

            numCalcdStructs=len(structInfo)

            # choose only structures which pass accept()
            structInfo = [s for s in structInfo if s.accept]

            #take top goodFraction structures
            from math import ceil
            if s.averageTopNum>=0:
                goodNum = s.averageTopNum
            else:
                goodNum = int(ceil(numCalcdStructs *
                                   s.averageTopFraction))
                pass

            structInfo = structInfo[0:min(goodNum,len(structInfo))]

            s.structInfo = structInfo
            pass

        
        s.lowEnergyStructure = comm.distribute(structInfo[0].id
                                               if len(structInfo) else -1)
        
        averageFitTo = None
        if s.lowEnergyStructure>=0 and s.averageFitSel:
            if s.esim and s.esim.size()>1:
                print("Warning: fitting EnsembleSimulation structures")
                pass
            from xplorSimulation import getXplorSimulation
            import simulation
            curSim=simulation.currentSimulation()
            curCoords=curSim.atomPosArr()
            if s.storeCoordinates:
                averageFitTo = comm.distribute(structInfo[0].atomPosArr)
            else:
                from pdbTool import PDBTool
                xpdb=PDBTool( s.filename(s.lowEnergyStructure))
                xpdb.read()
                pass
            averageFitTo=curSim.atomPosArr()
            curSim.setAtomPosArr(curCoords)

            def readFitWrite(structNum):
                """read structure, fit to lowest energy structure, updated
                to current context, and write out
                """
                s.structNum = structNum
                try:
                    #read files
                    pdbFile = s.pdbFile()
                    pdbFile.read()
                    sim.sync()
                    try:
                        AtomSel("known").apply( Fit(averageFitTo,
                                                    s.averageFitSel) )
                    except SystemError as e:
                        print("Error fitting structures:", e.args[0])
                        pass
                    s.averageContext()
                    s.inputPDB = s.pdbFileIn() # so the reported input file is ok
                    s.writeStructure(extraRemarks=s.getExtraRemarks(pdbFile),
                                     makeBackup=False)

                except IOError:
                    #file is missing, and can't be read
                    pass
                return
            
            procs = comm.barrier()

            # do the low energy structure first to avoid a race
            if s.lowEnergyStructure in s.structNums:
                readFitWrite( s.lowEnergyStructure )
                pass

            procs = comm.barrier()

            # now, fit structures to lowest energy structure and rewrite
            for structNum in s.structNums:
                if structNum!=s.lowEnergyStructure:
                    readFitWrite(structNum)
                    pass
                pass
            pass
        
        #3) barrier
        #  communicate these to proc 0
        #  proc 0 will then read in structures, calc ave and Fit
        from ensembleSimulation import singleThread, multiThread
        procs = comm.barrier()

        if s.processID!=0:
            #wait for collection of restraint info
            # for process 0, this is called in rStats.summarizeTerms() below
            if s.genViolationStats:
                rStats.collect(comm)
                if s.averageCrossTerms:
                    rStatsCross.collect(comm)
                    pass
                pass
            return #final regularization, reporting are done in proc 0
        
        structIDs=[]
        for proc in range(s.numProcs):
            pStructIDs=s.procStructMap[proc]
            if proc in procs:
                structIDs += pStructIDs
            else:
                host = comm.info(proc).remoteHost
                print("StructureLoop: no results obtained from ", \
                      "process %d running on %s" % (proc,host))
                print("\tskipping structures %s" % str(pStructIDs))
                pass
            pass
        if not s.storeCoordinates:
            for s.structNum in structIDs:
                try:
                    #reread files
                    s.pdbFile().read()
                    sim.sync()
                
                    structMap[s.structNum].atomPosArr = sim.atomPosArr()
                except IOError:
                    print('error reading', s.pdbFile().filename())
                    #file is missing, and can't be read
                    pass
                pass
            pass

        if s.averageFilename or s.genViolationStats:


            if not structInfo:
                print("no acceptable structures:")
                print("   averages and statistics not collected!")
                return

            #  - perform average
            ( aveCoords,fRMSD_array,cRMSD_array,bFactors, remarks ) = \
              calcAverageStruct( [s.atomPosArr for s in structInfo] ,
                                 fitSel=s.averageFitSel,
                                 compSel=s.averageCompSel,
                                 potList= s.averagePotList if
                                            s.averageRegularize else [],
                                 regularizeSteps=s.averageRefineSteps,
                                 averageRestrainSel=s.averageRestrainSel,
                                 averageRestrain=s.averageRestrain,
                                 fixedRegions=s.averageFixedRegions,
                                 rigidRegions=s.averageRigidRegions,
                                 )
                                  
               
        
            structureDetails = "%-30s  %10s  %9s  %9s\n" % ("","sort ",
                                                            "fit ",
                                                          "   comparison")
            structureDetails += "%-30s  %10s  %9s  %9s\n" % ("Filename:",
                                                             "energy",
                                                             "RMSD","RMSD")
            
            fRMSD_ave=0.
            cRMSD_ave=0.
            for i in range(len(structInfo)):
                fRMSD = fRMSD_array[i]
                cRMSD = cRMSD_array[i]
                fRMSD_ave += fRMSD
                cRMSD_ave += cRMSD
                
                structureDetails += "%-30s  %10.2f  %9.3f  %9.3f\n" % \
                                    (structInfo[i].filename,
                                     structInfo[i].sortEnergy,
                                     fRMSD,cRMSD)
                pass

            if len(structInfo):
                fRMSD_ave /= len(structInfo)
                cRMSD_ave /= len(structInfo)
                pass
            
            structureDetails += "\n%-30s  %10s  %9.3f  %9.3f\n" % \
                                ("average:","",fRMSD_ave,cRMSD_ave)
            if s.esim and s.esim.size()>1:
                structureDetails += "\nWARNING: ensemble simulation!\n"
                structureDetails += "   Structures are not fit thus"
                structureDetails += "   RMSD values are likely meaningless\n"
                pass

            s.fitRMSD  = fRMSD_ave
            s.compRMSD = cRMSD_ave
            ###################################################
            if s.averageRestrain:
                sim.setAtomPosArr(aveCoords)
                sortEnergy= s.averageSortPots.calcEnergy()
                sortViols = s.averageSortPots.violations()
                s.inconsistentAveStruct=True
                for struct in structInfo:
                    if (sortEnergy<=struct.sortEnergy and
                        sortViols <=struct.sortViols):
                        s.inconsistentAveStruct=False
                        pass
                    pass

                if s.inconsistentAveStruct:
                    structureDetails+= "\n"
                    structureDetails += "Warning: Calculated Average have energy or violation \n"
                    structureDetails += "greater than the calculated structures.\n"
                    structureDetails += "%-30s  %9.5s  %9.1s  \n" % ("Filename:",
                                                                     "Energy","Violations")
                    structureDetails += "%-30s  %10.2f  %9.3f  \n" % \
                                        ("average:",sortEnergy,sortViol)
                    structureDetails+="\n"
                    pass
                pass
            ####################################################################
            structureDetails += "\n  fit selection: "
            structureDetails += s.averageFitSel + '\n'
            structureDetails += "\n  comparison selection: "
            structureDetails += s.averageCompSel + '\n'
            pass
                
        if s.genViolationStats:

            out = "\n  Results for the top %d (of %d) structures\n\n" % \
                   (len(structInfo),numCalcdStructs)
            structNums=[s.id for s in structInfo]
            out+=rStats.summarizeTerms(comm,structNums)
            if s.averageCrossTerms:
                out += "\n Cross Validated Terms\n\n"
                out+=rStatsCross.summarizeTerms(comm,structNums)
                pass
            extraQuantityStats = rStats.summarizeExtraQuantities(structNums)
            if extraQuantityStats!='':
                out += "\n Statistics for Additional Quantities\n\n"
                out+=extraQuantityStats
                pass
            extraQuantityStats = rStatsCross.summarizeExtraQuantities(structNums)
            if extraQuantityStats!='':
                out += "\n Statistics for Additional Cross-Validated Quantities\n\n"
                out+=extraQuantityStats
                pass
            out += rStats.summarizeViolations(structNums)
            if s.averageCrossTerms:
                out += "\n Cross Validated Terms\n\n"
                out+=rStatsCross.summarizeViolations(structNums)
                pass

            out+="\n Energy terms used for sorting:\n "
            for pot in s.averageSortPots:
                out += " " + pot.instanceName()
                pass
            out += "\n"

            out += "\n" + structureDetails
            if singleThread():
                open(genFilename(s.pdbTemplate,'##')+".stats",
                     'w'                                     ).write(out)
                pass
            multiThread()
            s.restraintStats      = rStats
            s.restraintStatsCross = rStatsCross
            pass

        if s.averageFilename:
            sim.setAtomPosArr( aveCoords )

            remarks += structureDetails

            from atomSel import AtomSel
            from pdbTool import PDBTool
            avePDB = PDBTool(s.makeFilename(s.averageFilename),
                             s.writeSelection                 )
            for atom in AtomSel('known'):
                avePDB.setAux2(atom,bFactors[atom.index()])
            avePDB.addRemarks(remarks)
            s.averageContext()
            pl = s.averagePotList
            avePDB.addRemarks(analyze(pl,s.averageCrossTerms,
                 outFilename=s.makeFilename(s.averageFilename+'.viols')))
            if singleThread():
                avePDB.write()
                pass
            multiThread()
            pass
            
        return

    def pdbFileIn(s):
        """ return the input pdb filename corresponding to structNum. If
        running in the context of an EnsembleSimulation, the MEMBER literal
        in the filename is replaced with the appropriate ensemble member
        identifier.        
        """
        if not s.pdbFilesIn:
            return None
        if type(s.pdbFilesIn)==type("string"):
            filename=s.pdbFilesIn.replace("STRUCTURE",str(s.structNum))
        else:
            filename = s.pdbFilesIn[ s.structNum%len(s.pdbFilesIn) ]
            pass
        if s.esim:
            memberIndex=s.esim.member().memberIndex() % s.inMemberModulo
        else:
            memberIndex=0
            pass
        filename = filename.replace('MEMBER',str(memberIndex))
        return filename
    
    def pdbFile(s,nameSuffix=""):
        """ return a PDBTool object whose filename is generated by
            makeFilename() from the pdbTemplate argument of the class
            constructor. nameSuffix specifies an optional additional suffix
            for the generated filename.
            """
        #a reference to the PDBTool objet is saved here to allow
        #access like info.pdbFile().filename()
        # if this is not done, the PDBFile object is reaped before the
        # filename() string is returned: garbage out
        if (not hasattr(s,'pdbFileObj') or
            s.pdbFileObj.filename() != s.filename()):
            s.pdbFileObj = PDBTool( s.filename(suffix=nameSuffix),
                                    s.writeSelection             )
        return s.pdbFileObj
    def filename(s,
                 structNum=None,
                 memberNum=None,
                 suffix=""):
        """ return filename generated by makeFilename() from the pdbTemplate
        argument of the class constructor. If structNum is not specified, it
        defaults to s.structNum.
        """
        template = s.pdbTemplate
        if not template:
            template="SCRIPT_STRUCTURE_MEMBER.struct"
            pass
        return s.makeFilename(template,structNum,memberNum) + suffix
    def makeFilename(s,template,
                     structNum=None,
                     memberNum=None):
        """ create a filename given a template. In the template:
             the following substitutions are made:
                 SCRIPT    -> name of input script (without .py suffix)
                 STRUCTURE -> structNum if specified, else s.structNum
                 MEMBER    -> memberNum if specified, else ensemble member
                              index
        """
        from simulation import Simulation_currentSimulation
        sim = Simulation_currentSimulation()
        if structNum==None: structNum = s.structNum
        if sim.type() == "EnsembleSimulation":
            from ensembleSimulation import EnsembleSimulation_currentSimulation
            esim=EnsembleSimulation_currentSimulation()
            if memberNum==None: memberNum = esim.member().memberIndex()
        else:
            if memberNum==None: memberNum = 0
            pass
        return genFilename(template,structNum,memberNum)
    def getExtraRemarks(s,pdbFile):
        """
        Return those REMARKs records previously written by writeStructure
        using the extraRemarks argument.
        """
        remarks = pdbFile.remarks()
        foundHeader=False
        ret = []
        for remark in remarks:
            if foundHeader: ret.append(remark)
            if remark.startswith("StructureLoop Extra Remarks"):
                foundHeader=True
                pass
            pass
        return ret
                
            
    def writeStructure(s,potList=None, 
                       altPotList=None,
                       extraRemarks="",
                       nameSuffix="",
                       makeBackup=True):
        """perform analysis using analyze(), then write a structure with
        the analysis information included as remarks.

        The filename is generated from the pdbTemplate argument of the
        StructureLoop constructor.

        A summary is written out for each term in potList (and in altPotList,
        if specified), and more detailed violation information is output
        to a file named filename + '.viols' .

        If potList is not specified, it defaults to s.averagePotList.
        If altPotList is not specified, it defaults to s.averageCrossTerms.

        extraRemarks, if specified, is extra (string or a list of
        strings) information in the REMARKS section of the PDB file,
        printed after the usual summary.

        nameSuffix specifies an optional additional suffix for the
        generated filename.

        """
        if not potList:    potList = s.averagePotList
        if not altPotList: altPotList = s.averageCrossTerms
        import protocol
        remarks = "  generated by"
        remarks += " simulationTools.StructureLoop.writeStructure\n"
        if s.inputPDB:
            remarks += "  coordinates initialized from %s\n" % s.inputPDB
            pass
        remarks += "  seed info: initial: %d structure id: %d\n\n" % \
                   (protocol.initialRandomSeed(), s.structNum)

        if extraRemarks:
            if type(extraRemarks) in (type([]),type(tuple())):
                try:
                    extraRemarks='\n'.join(extraRemarks)
                except TypeError:
                    import traceback
                    traceback.print_exc()
                    print("continuing...")
                    pass
                pass
            remarks += "\nStructureLoop Extra Remarks\n" + extraRemarks
            pass

        from simulationTools import writeStructure
        return writeStructure(s.pdbFile(nameSuffix=nameSuffix),
                              potList,
                              altPotList,
                              s.filename()+ ".viols",
                              remarks,
                              makeBackup)
        return
    def structureNum(s):
        " return the current structure number, or -1."
        return s.structNum
    pass

class NullClass:
    """dummy empty class"""
    pass



class StructInfo:
    """
    helper class for generating structure statistics in StructureLoop.
    """
    def __init__(s,filename,
                 energy,
                 viols,
                 accept,
                 id,
                 atomPosArr=None):
        s.filename=filename
        s.sortEnergy=energy
        s.sortViols=viols
        s.accept=accept
        s.id=id
        if atomPosArr:
            s.atomPosArr = [tuple(v) for v in atomPosArr]
            pass
        return
    pass
        
def writeStructure(pdbFile,
                   potList,
                   altPotList=None,
                   analyzeFilename=None,
                   remarks="",
                   makeBackup=True,
                   clearRemarks=True):
    """perform analysis using analyze(), then write a structure with
    the analysis information included as remarks.
    
    The filename is generated from the pdbTemplate argument of the
    StructureLoop constructor.
    
    A summary is written out for each term in potList (and in altPotList,
    if specified), and more detailed violation information is output
    to a file named filename + '.viols' .

    If potList is not specified, it defaults to s.averagePotList.
    If altPotList is not specified, it defaults to s.averageCrossTerms.

    extraRemarks, if specified, is extra (string or a list of
    strings) information in the REMARKS section of the PDB file,
    printed after the usual summary.

    nameSuffix specifies an optional additional suffix for the
    generated filename.
    
    """
    import protocol
    extraRemarks = remarks
    remarks = analyze(potList,altPotList,
                      analyzeFilename)

    if remarks:
        remarks += "For help on header terms, try:\n"
        remarks += "   headerHelp TERM\n"
        try:
            from os import environ
            headerHelpURL = environ["XPLORNIH_HEADERHELP_URL"]
            remarks += "   or " +  headerHelpURL + '\n'
        except KeyError:
            pass
        pass

    if extraRemarks:
        if type(extraRemarks) in (type([]),type(tuple())):
            extraRemarks='\n'.join(extraRemarks)
            pass
        remarks += "\n" + extraRemarks + "\n"
        pass

    import os, pwd, time
    user="unknown"
    try:
        user = pwd.getpwuid(os.geteuid())[0]
    except (IndexError,KeyError):
        from os import environ
        try:
            user = environ["LOGNAME"]
        except KeyError:
            pass
        pass
    remarks += "-"*65 + '\n'
    remarks += "user: %-15s            date:" % user
    remarks += " " + time.asctime() + '\n'
    remarks += "-"*65 + '\n'
        

    if clearRemarks: pdbFile.clearRemarks()
    pdbFile.addRemarks( remarks )
    pdbFile.setMakeBackup(makeBackup)
    pdbFile.write()
    return

def convertToPotList(obj):
    """ convert a single potential term or a sequence to a
    <m potList>.PotList, if necessary. 
    """
    from potList import PotList
    if obj==None: obj=[]
    try:
        len(obj)
    except TypeError:
        obj=[obj]
        pass
    if not "PotList" in obj.__class__.__name__:
        pl = PotList()
        for p in obj:
            pl.append(p)
            pass
        obj = pl
        pass
    return obj
        

def calcAverageStruct(structs,fitSel,compSel,
                      potList=[],regularizeSteps=50,
                      averageRestrainSel="",averageRestrain=False,
                      fixedRegions=[],
                      rigidRegions=[]): 
    """compute unregularized average structure given structs.
    The structures are first fit using fitSel, and analysis is performed
    using compSel.

    For homogeneous <m ensembleSimulation>.EnsembleSimulations with Ne>1, this
    routine calculate the average of the ensemble averages.

    if potList is set, the sum of the specified terms 
    will be will be minimized with respect to the average coordinates, after
    the straight average structure has been calculated.

    fixedRegions and rigidRegions arguments are passed to regularizeRefine.
    """

    from cdsVector import vec_norm, sqrt
    from atomSelAction import Fit, RMSD
    from atomSel import AtomSel
    #Local copy of potlist
    avepotList=list(potList)
    from simulation import Simulation_currentSimulation
    sim = Simulation_currentSimulation()
                                      
    sim.setAtomPosArr( structs[0] )
    if averageRestrain:
       atmCoordList=[]
       atmCoordList.append(sim.atomPosArr())
    fitTo = ensembleAvePos()
    sim.setAtomPosArr( fitTo )
    aveCoords = sim.atomPosArr() # must be separate copy
    var = vec_norm(aveCoords)**2
    for struct in structs[1:]:
        sim.setAtomPosArr(struct)
        sim.setAtomPosArr( ensembleAvePos() )
        if fitSel:
            try:
                AtomSel("known").apply( Fit(fitTo,fitSel) )
            except SystemError as e:
                print("Error fitting structures:", e.args[0])
                pass
            pass
        if averageRestrain: 
           atmCoordList.append(sim.atomPosArr()) 
        aveCoords += sim.atomPosArr()
        var += vec_norm(sim.atomPosArr())**2
        pass
    if averageRestrain: 
       atmSel=averageRestrainSel
       #Calculating the target Grid
       from atomProb import AtomProb
       from selectTools import convertToAtomSel
       map = AtomProb(convertToAtomSel(atmSel),atmCoordList)      
       map.calc()
       targetMap=map.getGrid()
       sim.setAtomPosArr(aveCoords)
       from probDistPotTools import create_probDistPot
       prob=create_probDistPot("prob",targetMap,AtomSel(atmSel),\
                                     potType="cross_correlation",scale=100)
       avepotList.append(prob)
    aveCoords /= len(structs)
    var /= len(structs)
    bFactors = var - vec_norm(aveCoords)**2
    from math import pi
    bFactors.scale( 8*pi**2 )
    sim.setAtomPosArr( aveCoords )

    #get the tensor atoms in reasonable shape-
    #   averaging will have scrambled them
    from varTensorTools import calcTensor
    for t in getRegisteredTopoTerms("VarTensor",sim):
        calcTensor(t,suppressExceptions=True)
#        t.setFreedom("varyDa, varyRh")       #allow all tensor parameters float
        pass
    
    if len(avepotList)>0:
        minimizeRefine(avepotList,
                       fixedRegions=fixedRegions,
                       rigidRegions=rigidRegions,
                       refineSteps=regularizeSteps)
        pass

    # make sure final tensors are consistent with structure
    for t in getRegisteredTopoTerms("VarTensor",sim):
        calcTensor(t,suppressExceptions=True)
        pass
    
    aveCoords = sim.atomPosArr()
        

    fitTo=aveCoords
    aveRMSD=0.
    aveRMSDcomp=0.
    rmsdArray=[]
    rmsdCompArray=[]
    for struct in structs:
        sim.setAtomPosArr(struct)
        sim.setAtomPosArr( ensembleAvePos() )
        if fitSel:
            try:
                AtomSel("known").apply( Fit(fitTo,fitSel) )
            except SystemError as e:
                print("Error fitting structures:", e.args[0])
                pass
            pass
                    
        comparer=RMSD(fitTo)
        if fitSel:
            AtomSel(fitSel).apply(comparer)
            pass
        rmsd = comparer.rmsd()
        aveRMSD += rmsd
        rmsdArray.append( rmsd )
                    
        AtomSel(compSel).apply(comparer)
        rmsdComp = comparer.rmsd()
        rmsdCompArray.append( rmsdComp )
        aveRMSDcomp += rmsdComp
        pass
    
    aveRMSD      /= len(structs)
    aveRMSDcomp /= len(structs)

    remarks  = "average structure over %d files\n" % len(structs)
    remarks += "fitted using atoms: %s\n" % fitSel
    remarks += "RMSD diff. for fitted atoms: %f\n" % aveRMSD
    remarks += "comparison atoms: %s\n" % compSel
    remarks += "RMSD diff. for comparison atoms: %f\n" % aveRMSDcomp
    remarks += "B array (last column) is rms diff from mean\n"
    return (aveCoords,rmsdArray,rmsdCompArray,bFactors, remarks)

def ensembleAvePos():
    """return the (unregularized) average coordinates of the current ensemble
    For heterogeneous ensembles, this average is not meaningful, and the
    full ensemble is returned.
    """
    from simulation import Simulation_currentSimulation
    sim = Simulation_currentSimulation()
    esim=0
    if ( sim.type() == "EnsembleSimulation"):
        from ensembleSimulation import EnsembleSimulation_currentSimulation
        esim = EnsembleSimulation_currentSimulation()
        pass
    if not esim: return sim.atomPosArr()

    ave = esim.members(0).weight() * esim.members(0).atomPosArr()
    for i in range(1,esim.size()):
        if esim.members(i).numAtoms() != esim.members(i-1).numAtoms():
            return sim.atomPosArr()
        ave += esim.members(i).weight() * esim.members(i).atomPosArr()
        pass
    return ave
    

def minimizeRefine(potList,
                   refineSteps=50,
                   minimizeSteps=100,
                   xplorPots=['BOND','ANGL','IMPR'],
                   scaleMultiplier=0.001,
                   rigidRegions=(),
                   fixedRegions=(),
                   translateRegions=(),
                   simulation=None,
                   ):
    """ simple refinement using gradient minimization in Cartesian coordinates.
    Some potential terms are held fixed during minimization, while others are
    scaled up from a smaller value progressively, during each round of
    minimization.

    refineSteps specifies how many rounds of minimization to perform, while
    minimizeSteps are how many minimization steps are taken at each energy
    scale value.
    
    xplorPots are XPLOR terms which are to always be present, and for which
              the scale constant is held constant.
    scaleMultiplier specifies the initial value for the scale constant.

    rigidRegions specifies selections of atoms which do not move relative to
    each other. This argument can also be a single string.

    fixedRegions specifies selections of atoms which do not move at
    all. This argument can also be a single string.

    translateRegions specifies selections of atoms which can translate
    as a rigid body, but not rotate.
    """
    from simulation import currentSimulation
    currentSim = currentSimulation()
    newSim = simulation if simulation else currentSim
    import simulation
    simulation.makeCurrent(newSim)
    potList = convertToPotList(potList)
    flattenedPots = flattenPotList(potList)
    if type(fixedRegions)==type("string"): fixedRegions = [ fixedRegions ]
    if type(rigidRegions)==type("string"): rigidRegions = [ rigidRegions ]

    # refine here
    #  remove bond, angle, impr, vdw terms- if they exist- use them as is.
    # if they don't exist, add them in with default scale values.
    #
    # for rest of terms, loop over them, with MultRamp running from .01 .. 1
    # times the nominal scale values.
    hasReqdTerms = {}
    for p in xplorPots: hasReqdTerms[p] = False
    for pot in flattenedPots:
        if potType(pot) == 'XplorPot' and pot.instanceName() in xplorPots:
            hasReqdTerms[pot.instanceName()] = True
            pass
        pass
    addedPots=[]
    from xplorPot import XplorPot
    from avePot import AvePot
    for pType in xplorPots:
        # FIX: if this is required for EnsembleSimulations, it will break
        if not hasReqdTerms[pType]:
            addedPots.append(pType)
            potList.append( AvePot(XplorPot,pType) if
                            newSim.type()=="EnsembleSimulation" else
                            XplorPot(pType) )
        pass
    rampedParams = []
    for pot in potList:
        if potType(pot) == 'XplorPot' and pot.instanceName() in xplorPots:
            continue
        rampedParams.append( MultRamp( scaleMultiplier*pot.scale(),
                                       pot.scale(),
                                       "potList['%s'].setScale(VALUE)"%
                                       pot.instanceName() ) )
        pass
    
    from ivm import IVM
    import protocol

    minc = IVM(newSim)
    for aSel in rigidRegions: minc.group(aSel)
    for aSel in fixedRegions: minc.fix(aSel)
    for aSel in translateRegions:
        minc.group(aSel)
        minc.hinge('translate',aSel)
        pass
    protocol.cartesianTopology(minc)

    if refineSteps<=0:
        simulation.makeCurrent( currentSim )
        return
    
    #run initial powell minimization with Bond & Angles only - 
    # so impropers are well-defined.
    protocol.initMinimize(minc,
                           potList=[potList['BOND'],potList['ANGL']],
                           numSteps=200)

    minc.run()

    protocol.initMinimize(minc,
                          numSteps=minimizeSteps,
                          potList=potList,
                          dEPred=10)


    from simulationTools import AnnealIVM
    AnnealIVM(initTemp =0, finalTemp=0,
              numSteps=refineSteps,
              ivm=minc,
              rampedParams = rampedParams).run()

    #remove added pots
    for pot in addedPots:
        potList.remove(pot)
        pass
                             
    simulation.makeCurrent( currentSim )
    return

#used when creating temporary files to avoid races in EnsembleSimulation
#calculations
#threadSuffix=""
def mktemp(suffix='', prefix='tmp', dir=None):
    """ return a temporary file name which is unique for each process and
    thread within an ensemble simulation.
    """
    import tempfile
    from ensembleSimulation import threadSuffix
    suffix += '.' + threadSuffix
    return tempfile.mktemp(suffix,prefix,dir)

def potType(pot):
    """ return the potential's type.
    For most potential terms, the type is given by the potName() accessor,
    exception:
      for <m avePot>AvePot, potType(subPot()) is returned.
    """
    if pot.potName().startswith('ave_'):   return potType(pot.subPot())
    return pot.potName()
    
def flattenPotList(potList,
                   includePotLists=False):
    """
    Given a <m potList>.PotList, or other sequence of potential terms,
    return a list of potentials obtained from recursively the sequence,
    and any contained sequences.

    If includePotLists is True, include the non-toplevel PotLists themselves
    in the returned list.
    """
    from potList import PotListPtr
    ret = []
    for pot in potList:
        try:
            ret += flattenPotList(pot)
            if pot.potName()!='PotList' or includePotLists:
                ret.append(pot)
                pass
            pass
        except AttributeError:
            ret.append(pot)
            pass
        pass
    return ret

# this variable is a dictionary of dictionaries of lists of
# <m pot>.Pot objects. These are *not* PotProxy objects usually used
# and destroyed when their refcnt goes to zero.
registeredTopoTerms={}
def registerTopoTerm(term,
                     simulation=None):
    """
    add the given <m pot>.Pot to a list associated with its simulation.
    These objects will be automatically processed by topologySetup and
    massSetup.
    """

    #first, clean up any old mess
    gcRegisteredTopoTerms()
    
    potType = term.potName()
    
    if simulation==None: simulation = term.simulation()

    #only register PotProxy objects
    if hasattr(term,'__topoRegistered') and term.__topoRegistered==True:
        return
    
    term.__topoRegistered=True
        
    id = simulation.lookupID()
    if not id in registeredTopoTerms: registeredTopoTerms[id] = {}
    if not potType in registeredTopoTerms[id]:
        registeredTopoTerms[id][potType] = set()
        pass

    registeredTopoTerms[id][potType].add(term)

    from xplorSimulation import getXplorSimulation
    xsim = getXplorSimulation(simulation)

    #register to sub-XplorSimulation, if it's distinct
    if id!=xsim.lookupID(): registerTopoTerm(term,xsim)
    return

def getRegisteredTopoTerms(potType,
                           simulation=None):
    """
    return a list of all <m pot>.Pot objects with potName potType registered to
    the specified simulation.
    """
    ret=[]
    from simulation import currentSimulation
    if simulation==None: simulation = currentSimulation()
    try:
        ret += registeredTopoTerms[simulation.lookupID()][potType]
    except:
        pass
    try:
        if sim.type()=="EnsembleSimulation":
            from ensembleSimulation import EnsembleSimulation_currentSimulation
            ret += registeredTopoTerms[
       EnsembleSimulation_currentSimulation().member().lookupID()][potType]
            pass
        pass
    except:
        pass
    return ret

def gcRegisteredTopoTerms(verbose=False):
    """Remove any registered Topo terms without external references.
    Set verbose to True to debug garbage collection.
    """
    termCnt={}
    
    for simID in list(registeredTopoTerms.keys()):
        typeDict = registeredTopoTerms[simID]
        for termSet in list(typeDict.values()):
            terms=list(termSet)
            
            for term in terms:
                if term in list(termCnt.keys()):
                    termCnt[term] += 1
                else:
                    termCnt[term] = 1
                    pass
                pass
            pass
        pass

    #there may be multiple references, so need 2 passes to get correct
    # deletion reference count

    import sys
    for term in list(termCnt.keys()):
        # termCnt for sets, 1 to for function temporary
        # 1 for pot.instanceData, 2 for local
        
        if verbose:
            print('gcRegisteredTopoTerms: term: %s cnt: %d' %(
                term.instanceName(), sys.getrefcount(term)-termCnt[term]))
            pass
        if sys.getrefcount(term)-termCnt[term] >6:
            del termCnt[term]
        elif verbose:
            print('    removing ', term.instanceName())
            pass
        pass
    

    for simID in list(registeredTopoTerms.keys()):
        typeDict = registeredTopoTerms[simID]
        for termSet in list(typeDict.values()):
            terms=list(termSet)

            for term in terms:
                if term in list(termCnt.keys()):
                    termSet.remove(term)
                    pass
                pass
            pass
        pass
    return



def genFilename(template,
                structure=0,
                member=0):
        """ create a filename given a template. In the template:
             the following substitutions are made:
                 SCRIPT    -> name of input script (without .py suffix)
                 STRUCTURE -> the structure argument
                 MEMBER    -> the member argument
        """
        import re

        from sys import argv

        #convert to string
        if isinstance(structure,int): structure = "%d" % structure

        import xplor
        scriptName =  re.sub(r'\.[^.]*$','',xplor.scriptName)
        if not scriptName:
            scriptName='stdin'
            pass
        if scriptName.startswith("SCRIPT"):
            filename = template.replace("SCRIPT","%s" %scriptName)
        else:
            #not first, just use basename:
            import os
            baseName=os.path.basename(scriptName)
            filename = template.replace("SCRIPT","%s" %baseName)
            

        filename = filename.replace("STRUCTURE",structure)
        filename = filename.replace("MEMBER","%d" % member)
        return filename
    

def sortFilesByEnsMember(ensembleFiles):
    """
    given a list of files corresponding to an ensemble, extract the
    ensemble member numbers and use these to sort ensembleFiles and return
    the result.
    """

    ret = ensembleFiles

    #find the differing field in the files - presume this is
    #  the ensemble member index
    minBegField=len(ensembleFiles[0])
    for file in ensembleFiles[1:]:
        begField=0
        while begField<len(ensembleFiles[0]):
            begField+=1
            if not ensembleFiles[0].startswith(file[:begField]):
                break
            pass
        begField-=1
        minBegField=min(minBegField,begField)
        pass

    maxEndField=-len(ensembleFiles[0])
    for file in ensembleFiles[1:]:
        endField=-1
        while endField>-len(ensembleFiles[0]):
            endField-=1
            if not ensembleFiles[0].endswith(file[endField:]):
                break
            pass
        endField+=1
        maxEndField=max(maxEndField,endField)
        pass


    try:
        keyedEnsembleFiles={}
        for file in ensembleFiles:
            key = int(file[minBegField:maxEndField])
            keyedEnsembleFiles[key] = file
            pass
        ret = [keyedEnsembleFiles[key] for key in sorted(keyedEnsembleFiles)]
        pass
    except ValueError:
        pass
            
    return ret


class RunAction:
    """Base class for function or command string to be run in the defining
       frame at a later time. There is a value associated with this object which
       is always zero.
       
      update() - increments value. It will not change the value beyond
                 that specified by stopValue
      value()  - return the current value
      init(ns)        - set number of steps and initialize val to startValue
      finalize()      - set value to the final value.
    """
    def __init__(s,action):
        from inspect import currentframe, getouterframes
        s.action = action
        s.callingFrame = getouterframes( currentframe() )[2][0]
        s.val = 0
        s.startValue = 0
        s.stopValue  = 0
        s.numSteps=1
        s.curStep=0
        return
    def init(s,numSteps,caller=0):
        s.setNumSteps(numSteps)
        s.curStep=0
        s.val = s.startValue
        s.runAction(caller)
        return
    def finalize(s,caller=0):
        s.val = s.stopValue
        s.runAction(caller)
        return
    def runAction(s,caller=0):
        if not s.action: return
        s.fractionDone = float(s.curStep) / s.numSteps if s.numSteps else 1
        if type(s.action) == type("string"):
            global_dict = s.callingFrame.f_globals
            local_dict = s.callingFrame.f_locals
            local_dict["rampInfo"] = s
            local_dict["caller"] = caller
            exec( s.action.replace("VALUE","%e" %s.val),
                  global_dict, local_dict )
            del local_dict["rampInfo"]
        else:
            s.action(s.val)
            pass
        return
    def setNumSteps(s,numSteps):
        s.numSteps = numSteps
        return
    def value(s):
        return s.val
    def update(s,caller=0):
        s.runAction(caller)
        return 0.
    pass


class MultRamp(RunAction):
    """convenience class for multiplicatively (geometrically)
    ramping a value from startValue to stopValue over numberSteps
    constructor: MultRamp(startValue, stopValue, action)
    methods:
      update() - increments value. It will not change the value beyond
                 that specified by stopValue
      value()  - return the current value
      setNumSteps(ns) - set number of steps
      init(ns)        - set number of steps and initialize val to startValue
      finalize()      - set value to the final value.
    """
    def __init__(s,startValue,stopValue,action=None):
        RunAction.__init__(s,action)

        s.val = startValue
        s.startValue = max(float(startValue),1e-30)
        s.stopValue  = max(float(stopValue),1e-30)
        s.dirIncreasing=0
        if stopValue>startValue: s.dirIncreasing=1
        return
    def update(s,caller=0):
        s.val *= s.factor
        s.curStep += 1
        if (s.dirIncreasing and s.val>s.stopValue) or \
           (not s.dirIncreasing and s.val<s.stopValue):
            s.val = s.stopValue
            pass
        s.runAction(caller)
            
        return s.val
    def setNumSteps(s,numSteps):
        s.factor=1
        if numSteps: s.factor = (s.stopValue/s.startValue)**(1./numSteps)
        s.numSteps = numSteps
        return
    pass

class LinRamp(RunAction):
    """convenience class for linearly
    ramping a value from startValue to stopValue over numberSteps
    constructor: MultRamp(startValue, stopValue, action)
    methods:
      update() - increments value. It will not change the value beyond
                 that specified by stopValue
      value()  - return the current value
      setNumSteps(ns) - set number of steps
      init(ns)        - set number of steps and initialize val to startValue
    """
    def __init__(s,startValue,stopValue,action=None):
        RunAction.__init__(s,action)

        s.val = startValue
        s.startValue = float(startValue)
        s.stopValue  = float(stopValue)
        s.dirIncreasing=0
        if stopValue>startValue: s.dirIncreasing=1
        return
    def update(s,caller=0):
        s.val += s.factor
        s.curStep += 1
        if (s.dirIncreasing and s.val>s.stopValue) or \
           (not s.dirIncreasing and s.val<s.stopValue):
            s.val = s.stopValue
            pass
        s.runAction(caller)
            
        return s.val
    def setNumSteps(s,numSteps):
        s.factor=0
        if numSteps>0: s.factor = (s.stopValue-s.startValue)/numSteps
        s.numSteps = numSteps
        return
    pass

class StaticRamp(RunAction):
    """convenience class for static parameter setup.
      update() - increments value. It will not change the value beyond
                 that specified by stopValue
      value()  - return the current value
      setNumSteps(ns) - set number of steps
      init(ns)        - set number of steps and initialize val to startValue
    """
    def __init__(s,action,stride=1):
        """
        action is a function or string to be executed.
        stride specifies how often the function is called. A stride value of
        1 specifies that the function is called every time update is called.
        Larger values of stride specify that action is called by update()
        if caller.step%stride=0.
        """
        RunAction.__init__(s,action)
        s.stride=stride
        return
    def update(s,caller=0):
        s.curStep += 1
        if caller and 'step' in dir(caller) and caller.step%s.stride!=0:
            return
        s.runAction(caller)
            
        return s.val
    pass

class IVMAction(RunAction):
    """convenience class for static parameter setup.
      update() - increments value. It will not change the value beyond
                 that specified by stopValue
      value()  - return the current value
      setNumSteps(ns) - set number of steps
      init(ns)        - set number of steps and initialize val to startValue
    """
    def __init__(s,action):
        RunAction.__init__(s,action)
        return
    pass

class InitialParams:
    """Constructor takes a list of ramped parameters. The constructor invokes
    each parameter such that it set to its initial value, unless the optional
    call argument is specified as False.
    Also, this object can be called as a function with zero arguments, to
    set parameters to initial values. 
    """
    def __init__(s,pList,
                 call=True):
        s.fractionDone=0
        s.numSteps=1
        s.pList = pList
        if call: s.__call__()
        return
    def __call__(s):
        for p in s.pList:
            p.init(1,caller=s)
            pass
        return
    pass

class FinalParams:
    """constructor takes a list of ramped parameters. The constructor invokes
    each parameter such that it set to its final value, unless the optional
    call argument is specified as False.
    Also, this object can be called as a function with zero arguments, to
    set parameters to final values.
    """
    def __init__(s,pList,
                 call=True):
        s.fractionDone=1
        s.numSteps=1
        s.pList = pList
        if call: s.__call__()
        return
    def __call__(s):
        for p in s.pList:
            p.finalize(caller=s)
            pass
        return
    pass
        
#from potList import PotList
#class PotListWithContext(PotList):
#    """a <m potList>.PotList with an initializing function called before
#    calcEnergy/calcEnergyAndDerivs
#    """
#    def __init__(s,name="",
#                 potList=None,
#                 context=None):
#        if not potList: potList = PotList()
#        PotList.__init__(s,name)
#        if potList: s.copy(potList)
#        s.context = context
#        from inspect import currentframe, getouterframes
#        s.callingFrame = getouterframes( currentframe() )[1][0]
#        return
#    def calcEnergy(s):
#        if s.context: s.context()
#        return PotList.calcEnergy(s)
#    def calcEnergyAndDerivs(s,derivs):
#        if s.context: s.context()
#        return PotList.calcEnergyAndDerivs(s,derivs)
    



class AnnealIVM:
    """class to perform simulated annealing using molecular dynamics. """
    def __init__(s,  initTemp,  finalTemp,
                 ivm=0,
                 numSteps=None,   tempStep=None,
                 rampedParams   ={},
                 extraCommands  =0,
                 toleranceFactor=1000):
        """construct by specifying the intial and final annealing temperatures,
        and an <m ivm>.IVM object.

        if tempStep is specified, it will be used to determine the number of
          dynamics runs at different temperatures. If if is omitted, numSteps
          will be used (and tempStep implicitly determined).

        rampedParams is a list of MultRamp and LinRamp objects which specify
        refinement parameters to adjust during the simulated annealing run.
        extraCommands is a function or string which is run before dynamics at
        each temperature. If it is a function, is is passed the current
        AnnealIVM instance as the argument.

        During each temperature step the following instance variables are set:

          fractionDone - fraction (from 0 to 1) of annealing completed.
          bathTemp     - current temperature.
          step         - current step number.

        toleranceFactor is used to adjust the <m ivm>.IVM's energy tolerance
        as the temperature changes. The energy tolerance is calculated as
        
            eTolerance = temp / toleranceFactor

        The ivm argument doesn't actually have to be an <m ivm>.IVM, but it
        must be an object which has the following methods defined:
          setBathTemp
          setETolerance
          run

        """
        from inspect import currentframe, getouterframes
        s.initTemp = initTemp
        s.finalTemp = finalTemp
        s.extraCommands = extraCommands
        if tempStep:
            s.tempStep=tempStep
            s.numSteps = int( (initTemp - finalTemp)/float(tempStep) )
        elif numSteps:
            s.numSteps = numSteps -1 
            s.tempStep = (initTemp - finalTemp)/float(s.numSteps)
        else:
            raise("AnnealIVM: neither numSteps nor tempStep is defined")
        s.ivm = ivm
        s.params = rampedParams
        s.toleranceFactor = toleranceFactor
        s.callingFrame = getouterframes( currentframe() )[1][0]
        return
    def run(s):
        s.bathTemp = s.initTemp
        s.printTemp()

        s.step = 0
        s.initParameters()
        
        s.runExtraCommands()
        s.runIVM(s.bathTemp) #run at initial temperature
        
        for s.step in range(1,s.numSteps+1):
            s.bathTemp -= s.tempStep
            s.fractionDone = float(s.step) / (s.numSteps+1)
            s.printTemp()

            s.updateParameters()
            s.runExtraCommands()

            s.runIVM(s.bathTemp)

            pass
        return

    def printTemp(s):
        import simulationWorld
        simWorld = simulationWorld.SimulationWorld_world()
        if simWorld.logLevel() != 'none':
            print("AnnealLoop: current temperature: %.2f" % s.bathTemp)
            pass
        return

    def runExtraCommands(s):
        if s.extraCommands:
            if type(s.extraCommands) == type("string"):
                global_dict = s.callingFrame.f_globals
                local_dict = s.callingFrame.f_locals
                local_dict["annealLoopInfo"] = s
                exec( s.extraCommands, global_dict, local_dict )
                del local_dict["annealLoopInfo"]
            else:
                s.extraCommands(s)
                pass
            pass
        return
            
    def initParameters(s):
        "initialize ramped parameters"
        for param in s.params:
            param.init( s.numSteps, s ) #calls runAction
            pass
    def finalParameters(s):
        "sets parameters to final values"
        for param in s.params:
            param.val = param.stopValue
            param.runAction(s)
            pass
    def updateParameters(s):
        for param in s.params:
            param.update(s)
            pass
        return

    def runIVM(s,temp):

        if not s.ivm:
            return

        s.ivm.setBathTemp( temp )
        if s.toleranceFactor>0:
            s.ivm.setETolerance( temp / float(s.toleranceFactor) )
            pass
        s.ivm.run()
        
        return
    pass



def verticalFormat(name,
                   nameValuePairs,
                   minWidth=6,
                   numericFormat=".2f"):
    line2 = " " + name + ":"
    line1 = " " * len(line2)
    for (name,value) in nameValuePairs:
        strLen = max(minWidth,len(name))
        format = " %" + "%ds" % strLen
        line1 += format % name
        format = " %" + ("%d" % strLen) + numericFormat
        line2 += format % value
        pass
    return (line1, line2)

    
                   


if __name__ == "__main__":
    #
    # tests
    #

    
    from ensembleSimulation import EnsembleSimulation
    import sys
    import simulationWorld
    esim = EnsembleSimulation("esim",2)

    initSeed=751
    simWorld = simulationWorld.SimulationWorld_world()
    simWorld.setRandomSeed(751)

    result = ""
    expectedResult="""running structure 0 in process 0  thread: 0 seed: 751
running structure 1 in process 0  thread: 0 seed: 752
running structure 2 in process 0  thread: 0 seed: 753
running structure 3 in process 0  thread: 0 seed: 754
running structure 4 in process 0  thread: 0 seed: 755
running structure 5 in process 0  thread: 0 seed: 756
running structure 6 in process 0  thread: 0 seed: 757
running structure 7 in process 0  thread: 0 seed: 758
running structure 8 in process 0  thread: 0 seed: 759
running structure 9 in process 0  thread: 0 seed: 760
"""
    def concat(str):
        global result
        result += str
        return 0

    
    sys.stdout.write("StructureLoop: action as function...")
    StructureLoop(numStructures=10,
                  structLoopAction=lambda loopInfo: \
                  concat("running structure %d" % loopInfo.count +
                         " in process %d " % loopInfo.processID +
                         " thread: %d" % esim.member().memberIndex() +
                         " seed: %d\n" % loopInfo.randomSeed)
                  ).run()

    if result == expectedResult:
        print("ok")
    else:
        print("FAILED")
        pass


    simWorld.setRandomSeed(751)
    result=""
        
    sys.stdout.write("StructureLoop: action as string...")
        
    StructureLoop(numStructures=10,
                  structLoopAction=r'''global result
result += "running structure %d" % structLoopInfo.count + \
" in process %d " % structLoopInfo.processID + \
" thread: %d" % esim.member().memberIndex() + \
" seed: %d\n" % structLoopInfo.randomSeed
'''
                  ).run()
        
    if result == expectedResult:
        print("ok")
    else:
        print("FAILED")
        pass
            
    del esim

    sys.stdout.write("MultRamp: action as function...")
    expectedResult = "0.10 0.16 0.25 0.40 0.63 1.00 1.58 2.51 3.98 6.31 10.00 "
    result = ""
    def multProc(v):
        global result
        result +=  "%.2f " % v
        return
    
    param = MultRamp(0.1,10,multProc)
    numSteps=10
    param.init(numSteps)
    for i in range(numSteps):
        param.update()
        pass

    if result == expectedResult:
        print("ok")
    else:
        print("FAILED")
        pass
    

    sys.stdout.write("MultRamp: action as string...")
    expectedResult = "0.10 0.16 0.25 0.40 0.63 1.00 1.58 2.51 3.98 6.31 10.00 "
    result = ""
    def multProc(v):
        global result
        result +=  "%.2f " % v
        return
    
    param = MultRamp(0.1,10,"global result; result += '%.2f ' % VALUE")
    numSteps=10
    param.init(numSteps)
    for i in range(numSteps):
        param.update()
        pass

    if result == expectedResult:
        print("ok")
    else:
        print("FAILED")
        pass

    expectedResult="""0.10 extra: 1000.000000
0.32 extra: 900.000000
1.00 extra: 800.000000
3.16 extra: 700.000000
10.00 extra: 600.000000
"""

    sys.stdout.write("AnnealIVM: extraCmd as function...")
    result=""
        
    def coolProc(s):
        global result
        result += 'extra: %f\n' % s.bathTemp
        return

    AnnealIVM(initTemp=1000,
              finalTemp=600,
              tempStep=100,
              rampedParams=(param,),
              extraCommands=coolProc).run()

    if result == expectedResult:
        print("ok")
    else:
        print("FAILED")
        pass


    sys.stdout.write("AnnealIVM: extraCmd as string...")
    result=""
        
    AnnealIVM(initTemp=1000,
              finalTemp=600,
              tempStep=100,
              rampedParams=(param,),
              extraCommands=
              r"result += 'extra: %f\n' % annealLoopInfo.bathTemp").run()

    if result == expectedResult:
        print("ok")
    else:
        print("FAILED")
        pass

    pass

def testGradient(pots,
                 eachTerm=0,
                 alwaysPrint=0,
                 components=[0],
                 tolerance=1e-3,
                 eTolerance=1e-8,
                 epsilon=1e-7,
                 selString="all",
                 sim=0):
    """
    check the analytic gradient in the given potential terms against finite
    difference values.

    If the eachTerm argument is set, each term in the potList is tested
    individually.

    If alwaysPrint is set, all gradient terms will always be printed.

    components specificies which of the three gradient components x, y, z 
    (specified by the numbers 0,1,2) to test. By default this is just the
    x (0) component. A value of "random" specifies that a random component is
    tested for each atom. Random component specification is not yet supported
    for EnsembleSimulations.

    tolerance specifies the agreement expected between numerical and
    analytic gradient

    eTolerance specifies the ageement relative to the magnitude of the energy

    epsilon specifies the stepsize used in the finite different gradient
    calculation:

            dE/dqi(xyz) = E(qi(xyz)+epsilon) - E(qi(xyz)) / epsilon

    where E(0) is the energy evaluated at the nominal coordinate value
    qi(xyz).

    The Simulation is specified with the optional sim argument. The selString
    argument can be used to test a subrange of atoms' gradients.

    For each atom, the following is printed:
    atom identifying string  numerical gradient  analytical gradient
    """

    from potList import PotList
    potList = PotList()

    try:
        len(pots)
        for term in pots:
            potList.append(term)
            pass
        pass
    except TypeError:
        potList.append(pots)
        pass


    if not sim:
        from simulation import currentSimulation
        sim = currentSimulation()
        pass
    
    #is this an EnsembleSimulation?
    from ensembleSimulation import EnsembleSimulation_currentSimulation
    esim = EnsembleSimulation_currentSimulation()
    if esim and sim.name() == esim.name():
        sim = esim
    else:
        esim=0
        pass
    

    if eachTerm:
        ret=1
        for term in potList:
            print("testing gradient of potential term:", term.instanceName())
            if not testGradient(term):
                ret=0
            pass
        return ret
    
    from vec3 import Vec3

    # total gradient
    
    import random
    from derivList import DerivList
    dlist = DerivList()
    dlist.init(sim)
    energy = potList.calcEnergyAndDerivs(dlist)
    dlist.sync()
    derivs = dlist.get(sim)
    randComponents = True if components=="random" else False
    
    ret=1
    from atomSel import AtomSel
    if esim:
        header = "\n   Gradient Error Report\n\n"
        header += "%19s  ens %3s  %8s     %8s" %("Atom       ","xyz",
                                               "Numerical","From pot")
        isHeader=False
        sharedObj = esim.sharedObj()
        for j in range( esim.size() ):
            member = esim.members(j)
            sel = AtomSel(selString,member)
            for atom in sel:
                i = atom.index()
                for xyz in components:
                    if j == esim.member().memberIndex():
                        initPos = Vec3( member.atomPos(i) )
                        pos = Vec3( initPos )
                        pos[xyz] += epsilon
                        member.setAtomPos(i,pos)
                        pass
                    denergy = potList.calcEnergy()-energy
                    bad=False
                    if j == esim.member().memberIndex():
                        sharedObj.set( None )
                        grad = denergy/epsilon
                        bad=(abs(grad-derivs[i][xyz])>
                             tolerance*(1+
                                        max(abs(grad),abs(derivs[i][xyz]))) and
                             abs(grad-derivs[i][xyz])>eTolerance*energy)
                        if (alwaysPrint or bad):
                            sharedObj.set( (grad, derivs[i][xyz],bad) )
                            pass
                        member.setAtomPos(i,initPos)
                        pass
                    result=sharedObj.barrierGet()
                    if result!=None:
                        ret=0
                        if not isHeader:
                            print(header)
                            isHeader=True
                            pass
                        lgrad,lderivs,lbad=result
                        print("%-25s"%member.atomByID(i).string(),\
                              j,xyz,lgrad,lderivs, end=' ')
                        if alwaysPrint and lbad:
                            print("  **ERROR**")
                        else:
                            print()
                            pass
                        
                        pass
                    pass
                pass
            pass
        pass
    else:
        #not an EnsembleSimulation
        header = "\n   Gradient Error Report\n\n"
        header += "%19s  %3s %8s     %8s" %("Atom       ","xyz",
                                        "Numerical","From pot")
        isHeader=False
        for atom in AtomSel(selString,sim):
            i = atom.index()
            if randComponents:
                components=[random.choice([0,1,2])]
                pass
            for xyz in components:
                initPos = Vec3( sim.atomPos(i) )
                pos = Vec3( initPos )
            
                pos[xyz] += epsilon
                sim.setAtomPos(i,pos)
                denergy = potList.calcEnergy()-energy
                grad = denergy/epsilon
                bad = (abs(grad-derivs[i][xyz]) >
                        tolerance*(1+
                                   max(abs(grad),abs(derivs[i][xyz]))) and
                        abs(grad-derivs[i][xyz]) > tolerance and
                        abs(grad-derivs[i][xyz])>eTolerance*abs(energy))
                if (alwaysPrint or bad):
                    if not isHeader:
                        print(header)
                        isHeader=True
                        pass
                    ret=0
                    print("%s %3d %10.5f %10.5f" %(sim.atomByID(i).string(),
                                                   xyz,
                                                   grad,
                                                   derivs[i][xyz]), end=' ')
                    if alwaysPrint and bad:
                        print("  **ERROR**")
                    else:
                        print()
                        pass
                    pass
                sim.setAtomPos(i,initPos)
                pass
            pass
        pass
    return ret


def saRefine(potList,
             refineSteps=50,
             xplorPots=['BOND','ANGL','IMPR','RAMA'],
             initTemp=10,
             finalTime=0.2,
             numSteps=100,
             htFinalTime=10,
             htNumSteps=1000,
             initVel=1,
             scaleMultiplier=0.001,
             rigidRegions=(),
             fixedRegions=(),
                   ):
    """ 
    Added by Robin A Thottungal 
    Add explanation 07/16/09
    refineSteps specifies how many rounds of minimization to perform.
    xplorPots are XPLOR terms which are to always be present, and for which
              the scale constant is held constant.
    scaleMultiplier specifies the initial value for the scale constant.

    initTemp specifies the initial temperature for high-temperature dynamics,
    and at the start of simulated annealing. The parameters finalTime, and
    numSteps specify dynamics duration and number of steps at each step of
    simulated annealing, while htFinalTime and htNumSteps specify these
    parameters for initial dynamics.

    rigidRegions specifies selections of atoms which do not move relative to
    each other.

    fixedRegion specifies selections of atoms which do not move at all.
    """
    pots = flattenPotList(potList)

    #first get the tensor atoms in reasonable shape-
    # averaging will have scrambled them
    from varTensorTools import getVarTensors, calcTensor
    varTensors = getVarTensors(pots)

    for t in varTensors:
        calcTensor(t,suppressExceptions=True)
        pass


    # refine here
    #  remove bond, angle, impr, vdw terms- if they exist- use them as is.
    # if they don't exist, add them in with default scale values.
    #
    # for rest of terms, loop over them, with MultRamp running from .01 .. 1
    # times the nominal scale values.
    from potList import PotList
    minPots = PotList()
    hasReqdTerms = {}
    for p in xplorPots: hasReqdTerms[p] = 0
    rampedParams = []
    for pot in pots:
        reqdTerm=0
        for pType in xplorPots:
            if potType(pot) == 'XplorPot' and pot.instanceName() == pType:
                minPots.append(pot)
                hasReqdTerms[pType] = 1
                reqdTerm=1
                continue
            pass
        if reqdTerm: continue
        minPots.append( pot )
        rampedParams.append( MultRamp( scaleMultiplier*pot.scale(),
                                       pot.scale(),
                                       "minPots['%s'].setScale(VALUE)"%
                                       pot.instanceName() ) )
        pass
    from xplorPot import XplorPot
    for pType in xplorPots:
        if not hasReqdTerms[pType]: minPots.append( XplorPot(pType) )
        pass
    """
    #Added by Robin A Thottungal on 05/11/09 for
    #running a powell minimization with Bond & Angles
    # to fix up the improper
    #begin:
    from ivm import IVM
    import protocol
    min=IVM()
    protocol.cartesianTopology(min)
    protocol.initMinimize(min,
                           potList=[XplorPot("BOND"),XplorPot("ANGL")],
                           numSteps=200)

    if refineSteps>0:
        min.run()
    #end
    """    
    #Cartesian topology
    from ivm import IVM
    minc = IVM()
    for aSel in rigidRegions: minc.group(aSel)
    for aSel in fixedRegions: minc.fix(aSel)
    import protocol
    protocol.cartesianTopology(minc)

    #high-temp dynamics
    protocol.initDynamics(minc,
                          potList=minPots,
                          bathTemp=initTemp,
                          initVelocities=initVel,
                          finalTime=htFinalTime,# stops at 800ps or 8000 steps
                          numSteps=htNumSteps,  # whichever comes first
                          printInterval=100)

    InitialParams(rampedParams)
    minc.run()
    
    #simulated annealing
    protocol.initDynamics(minc,
                          bathTemp=initTemp,
                          finalTime=finalTime,
                          numSteps=numSteps,
                          initVelocities=initVel,
                          potList=minPots)


    import varTensorTools
    for m in varTensors:
        m.setFreedom("varyDa, varyRh")       #allow all tensor parameters float
        varTensorTools.topologySetup(minc,m) #setup tensor topology
        pass

    if refineSteps>0:
        from simulationTools import AnnealIVM
        AnnealIVM(initTemp, finalTemp=0,
                  numSteps=refineSteps,
                  ivm=minc,
                  rampedParams = rampedParams).run()
        pass
                             
    # make sure final tensor is consistent with structure
    for t in varTensors:
        calcTensor(t,suppressExceptions=True)
        pass

    return

#backwards compatibility
from termAnalysis import registerTerm, registerExtraStats, analyze
from termAnalysis import getPotTerms, summarizePotentials
