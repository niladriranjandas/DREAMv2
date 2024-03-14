"""High-level initialization routines.

These routines provide high level interfaces to initialize potential terms
and to set up minimization and dynamics calculations.

The following functions are imported from the <m regularize> module:
fixupCovalentGeom, covalentMinimize, addUnknownAtoms.

Global variable pdbLocation: location from which PDB entries are loaded.
"""
from atomSel import AtomSel

#default topology/parameter files
parameters={}
parametersInitialized={}
topology={}
linkages={}
topologyInitialized={}
parVersion={}
topVersion={}

systems=('protein','nucleic','water','metal','axis')

parameters['protein'] = "protein.par"
parametersInitialized['protein']=[]
topology['protein'] = "protein.top"
linkages['protein'] = []
topologyInitialized['protein']=[]

parameters['nucleic'] = "nucleic.par"
parametersInitialized['nucleic']=[]
topology['nucleic'] = "nucleic.top"
linkages['nucleic'] = []
topologyInitialized['nucleic']=[]

parameters['water'] = "tip3p.parameter"
parametersInitialized['water']=[]
topology['water'] = "tip3p.topology"
linkages['water'] = []
topologyInitialized['water']=[]

parameters['metal'] = "ion.par"
parametersInitialized['metal']=[]
topology['metal'] = "ion.top"
linkages['metal'] = []
topologyInitialized['metal']=[]


#pseudo-atoms for rdc, pre, etc pots.
parameters['axis'] = "axes.par"
parametersInitialized['axis']=[]

from regularize import fixupCovalentGeom, covalentMinimize, addUnknownAtoms
from regularize import CovalentViolation

def initRandomSeed(seed=None):
    """set the initial random seed. If the seed argument is omitted, it is set
    to seconds from the epoch mod 10^7
    """
    if 'initialRandomSeed_' in globals():
        raise Exception("random seed already initialized")
    global initialRandomSeed_
    import time
    if seed==None:
        seed = int(time.time() % 1e7)
        import xplor
        seed = xplor.p_comm.distribute(seed)
        pass
    initialRandomSeed_ = seed
    import simulationWorld
    simWorld = simulationWorld.SimulationWorld_world()
    simWorld.setRandomSeed( initialRandomSeed_ )
    import ensembleSimulation
    if ensembleSimulation.currentSimulation():
        print("Warning: initRandomSeed should be called before creating an", end=' ')
        print(" EnsembleSimulation")
        pass
    if simWorld.logLevel() != 'none':
        print("random seed initialized to ", initialRandomSeed_)
        pass
    return

def initialRandomSeed():
    """return the initial random seed value, or that from
    <m simulationWorld>.random if initRandomSeed has not been called.
    """
    if 'initialRandomSeed_' in globals():
        return initialRandomSeed_
    else:
        import simulationWorld
        simWorld = simulationWorld.SimulationWorld_world()
        return simWorld.random.seed()
    pass

def genTopParFilename(filename,
                      useTopParVar=False,
                      suffix=""):
    """
    If the given filename exists, it is returned, else the file is
    looked for the in the TOPPAR directory, and a full path is returned
    if it is found there.

    If useTopParVar=True, the XPLOR TOPPAR variable is used for the
    filename (instead of expanding its value). This results in
    shorter paths, but the path is only valid within the XPLOR interpreter.

    If suffix is specified, the filename will be tried with "." + suffix
    appended.
    """
    import os
    from os import environ as env
    ret=filename
    found=False
    if os.access(filename,os.F_OK):
        found=True
        ret=filename
        pass
    elif os.access(env['TOPPAR'] + '/' + filename,os.F_OK):
        found=True
        if useTopParVar:
            ret = 'TOPPAR:' + filename
        else:
            ret = env['TOPPAR'] + '/' + filename
            pass
        pass
    elif suffix:
        ret = filename + "." + suffix
        if os.access(ret,os.F_OK):
            found=True
            pass
        elif os.access(env['TOPPAR'] + '/' + ret,os.F_OK):
            found=True
            if useTopParVar:
                ret = 'TOPPAR:' + ret
            else:
                ret = env['TOPPAR'] + '/' + ret
                pass
            pass
        pass
    if not found:
        mess = """genTopParFilename: could not find name %s as a structure type
        or as file in the current dir, or in TOPPAR""" % filename
        raise Exception(mess)
    return ret
            

from xplorSimulation import getXplorSimulation    

def initParams(files=[],
               string="",
               reset=0,
               use_dihe=False,
               weak_omega=0,
               simulation=0,
               system=None,
               silent=False):
    """files is a structure type or filename or a list of a combination of
    these two.

    valid structure types are: protein, nucleic, and water

    string is an optional string of XPLOR-formatted parameter statements.

    The default protein, nucleic and water parameter files are specified in
    the protocol.parameters dictionary. Default values can be modified from
    a script using e.g.
       import protocol
       protocol.parameters['protein'] = '/path/to/other/file.par' 

    If the argument is a filename, the current directory is first searched
    for the file, and then the TOPPAR directory is searched.
    XPLOR parameters are documented <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/node49.html
    here>.
    If reset is true, all previous parameter settings are discarded
    If weak_omega is true, improper force constants associated with the
    peptide bond are reduced by 1/2 to allow a small amount of flexibility.

    If use_dihe is True, the DIHE XPLOR term is enabled.

    If system is specified, specify that the parameters are forthe specified
    system. Valid values are 'protein', 'nucleic', 'water', and 'metal'.

    If silent is True, XPLOR messages are suppressed.
    """
    import os
    from os import environ as env

    xSim = getXplorSimulation(simulation)
    command = xSim.fastCommand
    simName = xSim.name()
    if silent: outputState = xSim.disableOutput()

    if reset:
        command("param reset end")
        command("eval ($edtaParamsInit=FALSE)")
        command("eval ($ionParamsInit=FALSE)")
        for key in list(parametersInitialized.keys()):
            if simName in parametersInitialized[key]:
                parametersInitialized[key].remove(simName)
                pass
            pass
        pass

    if weak_omega:
        command('eval ($weak_omega=1)')

    if use_dihe:
        command('eval ($use_dihe=1)')

    def addParams(file,system):
        if file in parameters:
            if simName in parametersInitialized[file]: return
            parametersInitialized[file].append(simName)
            system=file
            file = parameters[file]
            pass

        file = genTopParFilename(file,True,"par")
            
        varName=''
        if system:
            varName="%s_par_vers"%system
            pass
        (version) =command("param @%s end" %file,varName)
        loutputState = xSim.disableOutput()
        for system in systems:
            varName=''
            varName="%s_par_vers"%system
            ans = command('',(varName,))
            if ans[0]!='':
                parVersion[system] = ans[0]
                pass
            pass
        xSim.enableOutput(loutputState)
        return


    if type(files)==type("string"):
        files=[files]
        pass
    tmpFilename=None
    if string!="":
        from simulationTools import mktemp
        tmpFilename=mktemp(prefix="params")
        open(tmpFilename,"w").write(string)
        files.append(tmpFilename)
    for file in files:
        addParams(file,system)
        pass
    if tmpFilename:
        import os
        os.unlink(tmpFilename)
        pass

    if silent: xSim.enableOutput(outputState)
    return

def initTopology(files=[],
                 reset=0,
                 simulation=0,
                 system=None):
    """file is a structure type or filename or a list of a combination of
    these two.

    valid structure types are: protein, nucleic, water, and metal.

    The default protein, nucleic and water topology files are specified in
    the protocol.topology dictionary. Default values can be modified from
    a script using e.g.
       import protocol
       protocol.topology['protein'] = '/path/to/other/file.top' 

    First the current directory is searched for the file, and then TOPPAR
    is searched.
    XPLOR topology is documented <l http://nmr.cit.nih.gov/xplor-nih/doc/current/xplor/node46.html
    here>.
    If reset is true, all previous topology settings are discarded.

    A non-current <m simulation>.Simulation can be specified with the simulation
    argument.

    If system is specified, specify that the topology are for the specified
    system. Valid values are 'protein', 'nucleic', 'water', and 'metal'.

    """
    import os
    from os import environ as env

    xSim = getXplorSimulation(simulation)
    command = xSim.fastCommand
    simName = xSim.name()

    if reset:
        command("rtf reset end")
        command("eval ($ionTopoInit=FALSE)")
        for key in list(topologyInitialized.keys()):
            if simName in topologyInitialized[key]:
                topologyInitialized[key].remove(simName)
                linkages[key]=[]
                pass
            pass
        pass


    def addTopology(file,system):
        system=""
        if file in topology:
            if simName in topologyInitialized[file]: return
            topologyInitialized[file].append(simName)
            system=file
            file = topology[file]
            pass

        xfile = genTopParFilename(file,True,"top")

        varName=''
        if system:
            varName="%s_top_vers"%system
            pass
        (version) = command("rtf @%s end" %xfile,varName)
        if system:
            file = genTopParFilename(file,False,"top")
            #FIX: check for mismatch
            topVersion[system] = version[0]
            linkages[system] += [line[1:] for line in open(file).readlines()
                                 if line.startswith("!LINK")]
            pass
        return


    if type(files)==type("string"):
        files=[files]
        pass
    for file in files:
        addTopology(file,system)
        pass
    return

def commandNoEcho(xSim,cmd,verbose=None):
    """run the given XPLOR command in the specified
    <m xplorSimulation>.XplorSimulation.
    """
    import simulationWorld
    simWorld = simulationWorld.SimulationWorld_world()
    if verbose==None and simWorld.logLevel() != 'none':
        verbose=1
        pass
    if not verbose:
        outputState=xSim.disableOutput()
        pass
    xSim.command(cmd)
    if not verbose:
        xSim.enableOutput(outputState)
        pass
    return
    

def initStruct(files=0,
               string=None,
               erase=1,
               simulation=0):
    """Read XPLOR PSF files.

    The files argument is a filename (a string) or sequence of filenames to read.
    First the current directory is searched for the file, and then TOPPAR
    is searched.  Alternatively, a psf entry (starting with PSF...) can be
    directly passed as the files argument. In addition, a string value of the
    PSF can be passed via the string argument.

    Any pre-existing structure information is erased unless the erase argument
    is cleared. 
    """
    xSim = getXplorSimulation(simulation)

    if erase: commandNoEcho(xSim,"struct reset end",verbose=False)


    if string:
        commandNoEcho(xSim,"struct %s end" % string)
        pass

    if not files: return

    import re, os
    done=False
    if type(files)==type("string"):
        if re.search("PSF\s*\n",files):
            commandNoEcho(xSim,"struct %s end" % files)
            done=True
        else:
            files=[files]
            pass
        pass
    if not done:
        for file in files:
            try:
                os.stat(file)
                commandNoEcho(xSim,"struct @%s end" %file)
            except OSError:
                from os import environ as env
                tfile = env['TOPPAR'] + '/' + file
                try:
                    os.stat(tfile)
                    commandNoEcho(xSim,"struct @TOPPAR:%s end" %file)
                except OSError:
                    raise Exception("initStruct: could not find file " + file +\
                                    " in current dir, or in TOPPAR")
                pass
            pass
        pass
    updatePseudoAtoms(simulation)
    return

pdbLocation="ftp://ftp.wwpdb.org/pub/pdb/data/structures/divided/pdb/"

downloadPDB_connectAttempts=3
def downloadPDB(entry):
    """download a compressed PDB file from the RCSB database
    return the PDB record as a string. If the download fails, the function
    will sleep for three seconds and then retry. The total number of connection
    attempts is set via the module scope variable downloadPDB_connectAttempts
    whose default value is 3.
    """

    entry = entry.lower()

    if len(entry)!=4:
        raise Exception("entry must be a 4 character string:" + repr(entry))

    middle = entry[1:3]

    url = pdbLocation + middle + '/pdb' + entry + '.ent.gz'

    import urllib.request, urllib.error, urllib.parse
    import time
    for attempt in range(downloadPDB_connectAttempts):
        try:
            u=urllib.request.urlopen(url)
            contents=u.read()
            u.close()
            break
        except IOError as err:
            print('sleeping...')
            time.sleep( 3 ) # secs -
            pass
        pass
    else:
        raise Exception("error fetching %s: %s"%(url,str(err))+
                        repr(err.args))
    
    from io import BytesIO
    import gzip
    ret=gzip.GzipFile(fileobj=BytesIO(contents)).read()
    ret=str(ret,'utf8')

    return ret

def splitModel(filename,
               explicitModel):
    """
    Process a filename with a trailing model number separated by a ``:''
    character from the preceeding filename. The value of explicitModel
    overrides the colon model if it is not negative.
    """
    fields = filename.split(':')
    model=fields[-1]
    if len(fields)==1:
        return filename,explicitModel

    filename=":".join(fields[:-1])
    if explicitModel>=0:
        return filename,explicitModel
    try:
        model=int(model)
    except ValueError:
        return filename,explicitModel

    return filename,model
            

def loadPDB(file='',
            string='',
            entry='',
            model=-1,
            simulation=None,
            correctSymmetricSidechains=True,
            deleteUnknownAtoms=False,
            maxUnreadEntries=(20,0.2),
            processSSBonds=True,
            processBiomt=False,
            loadOccupanciesBfactors=False,
            verbose=-1,
            **kwargs):
    """Load a PDB entry from file, string, or from the RCSB database.
    autogenerate the PSF information and initialize coordinates.

    deleteUnknownAtoms and maxUnreadEntries are passed to initCoords().

    processSSBonds specifies whether the SSBOND pdb record is used in PDB
    generation, and is passed to <m psfGen>.pdbToPSF.

    processBiomt specifies whether the REMARK 350 BIOMT record is used. By
    default, Biomolecule 1 is read, specifying a different integer to the
    processBiomt argument will read that entry.

    loadOccupanciesBfactors is passed to initCoords.

    The following optional arguments are passed through to other functions:

      swapProtons and assumeSegid are handled by the function
      matchInexactAtomEntry.

      includeHETATM, fixMethylImpropers and fixProtonImpropers are handled by
      the function initCoords.

      suppressExceptions and failIfInsertion are passed through to
      psfGen.pdbToSeq.

    A side-effect is that protein and/or nucleic force field parameters are
    loaded by initParams().

    This function returns the InitCoordsResult from initCoords (described
    below). 
    """
    sim = simulation
    import simulation
    if sim:
        oldCurrentSimulation = simulation.currentSimulation()
        simulation.makeCurrent(sim)
        pass

    
    def finalize():
        updatePseudoAtoms()
        if sim:
            simulation.makeCurrent(oldCurrentSimulation)
            pass
        return

    file,model = splitModel(file,model)

    logfile=open("/dev/null","w")
    from simulationWorld import world
    import sys
    if verbose<0 and world().logLevel()!="none": verbose=True
    if verbose==True: logfile=sys.stdout

    if file:
        print('loading pdb file: %s ' % file, end=' ', file=logfile) ; logfile.flush()
        pdbEntry = open(file).read()
    elif string:
        pdbEntry = string
    elif entry:
        print('loading pdb entry: %s ' % entry, end=' ', file=logfile) ; logfile.flush()
        
        pdbEntry = downloadPDB(entry)

        print(' [downloaded]', end=' ', file=logfile) ; logfile.flush()
    else:
        raise Exception("entry or file argument must be specified.")
        
    
    import psfGen
    biomt=psfGen.pdbToPSF(pdbEntry,
                          processSSBonds=processSSBonds,
                          processBiomt=processBiomt,
                          **kwargs)
    kwargs.pop('failIfInsertion',None) # used in pdbToPSF
    kwargs.pop('suppressExceptions',None) # used in pdbToPSF
    kwargs.pop('deprotonateHIS', None) # used in pdbToPSF

    print(' [psf]', end=' ', file=logfile) ; logfile.flush()

    includeHETATM=kwargs.pop('includeHETATM', False)
    try:
        ret=initCoords(string=pdbEntry,model=model,verbose=verbose,
                       correctSymmetricSidechains=correctSymmetricSidechains,
                       includeHETATM=includeHETATM,
                       deleteUnknownAtoms=deleteUnknownAtoms,
                       maxUnreadEntries=maxUnreadEntries,
                       loadOccupanciesBfactors=loadOccupanciesBfactors,
                       biomt=biomt,**kwargs)
    except:
        finalize()
        raise

    print(' [coords]', file=logfile)

    finalize()

    return ret

def initBiomtCoords(biomt,
                    simulation=None):
    """Fill in coordinates from a <m psfGen>.Biomt object.

    Note that occupancy and bfactor values for duplicated coordinates will
    not be filled.
    """
    sim = simulation
    if not sim:
        import simulation
        sim=simulation.currentSimulation()
        pass
    
    from atomSelAction import Rotate, Translate
    for c in biomt.chains:
        fromAtoms=AtomSel('segid "%s"' % c,sim)
        startIndex=fromAtoms[0].index()
        stopIndex=fromAtoms[-1].index()
        for l in biomt.labels():
            atomPosArr=sim.atomPosArr()
            toAtoms=AtomSel('segid "%s%s"' % (c,l),sim)
            if len(fromAtoms) != len(toAtoms):
                raise Exception("selections have different numbers of atoms")
            newStartIndex=toAtoms[0].index()
            
            for index in range(startIndex,stopIndex+1):
                coord = atomPosArr[index]
                atomPosArr[index+newStartIndex-startIndex] = coord
                pass
            sim.setAtomPosArr(atomPosArr)
            toAtoms.apply( Rotate(biomt.entry(l).rot) )
            toAtoms.apply( Translate(biomt.entry(l).trans) )
            pass
        pass
    
    return



class InitCoordsResult:
    """ A class to simplify returning B factors, occupancies and remarks read
    during a call to initCoords. It contains two members, bfactors and
    occupancies, which are arrays which contain the respective properties
    from the PDB record for each atom references by its index."""
    def __init__(self):
        self.bfactors = []
        self.occupancies = []
        self.remarks = []
        pass
    pass


from trace import notrace_decorate
@notrace_decorate
def initCoords(files=[],string="",
               model=-1,
               verbose=-1,
               erase=False,
               useChainID=True,
               includeHETATM=False,
               correctSymmetricSidechains=True,
               strictResNames=False,
               maxUnreadEntries=(20,0.2),
               deleteUnknownAtoms=False,
               loadOccupanciesBfactors=False,
               biomt=None,
               selection="all",
               **kwargs):
    """
    Initialize coordinates from one of more pdb file,
      or from a string containing a PDB entry.

    If model is specified, the specified PDB MODEL record will be read.

    If erase is set to True, then atom positions in selection are cleared
    before reading the pdb record.

    If useChainID is True (the default), then the chain ID PDB field is
    used as the segment name (if set), If the chain ID field is blank it is
    ignored regardless.

    By default, HETATM PDB records are ignored. If includeHETATM=True, they
    will be read and treated as ATOM records.

    If deleteUnknownAtoms is True, then atoms whose coordinates are still
    unknown after the PDB file is read will be deleted.
    
    If correctSymmetricSidechains is True (the default), then the sidechains
    of phe, tyr, asp, and glu residues are checked to ensure that their
    symmetric torsion angles (chi3 for glu, chi2 for the rest) are in
    the correct (0..180 degree) range. If they are incorrect, then atomic
    positions are swapped using <m selectTools>.correctSymmetricSidechains.

    If strictResNames is True, then residue names in the PDB must match
    those in the PSF, otherwise, the residue names need not match
    and backbone atoms for mutants can be read. Note that sidechain atoms (or
    base atoms in nucleic acids) may be strangely placed if this option is
    False.

    PDB ATOM records which do not exactly match the expected input
    (the atom name is not exactly the same) are attempted to be
    matched using the function matchInexactAtomEntry. The argument
    maxUnreadEntries specifies the behavior if many unmatachable
    records are encountered. If maxUnreadEntries is an integer, an
    exception is thrown if more than that number of unmatchable
    records are encountered. If maxUnreadEntries is a two-membered
    tuple, the exception is not thrown if number of unmatchable
    records is less than the first entry or the number of unmatchable
    records is less than the second entry (a fraction) multiplied by
    the number of atoms. This check is disabled (no exception will be
    thrown) if maxUnreadEntries is set to None.

    If loadOccupanciesBfactors is True, occupancies and bfactors will be
    passed to the appropriate <m xplorSimulation>.XplorSimulation. If
    specified as True, the <m atomSelLang> facility will support searching
    the b and q attritbutes.

    biomt is an optional <m psfGen>.Biomt object used to generate identical
    copies of coordinates in ATOM records.
    
    If fixMethylImpropers is True, <m selectTools>.fixMethylImpropers will be
    called after coordinates are loaded.

    If fixProtonImpropers is True, <m selectTools>.fixProtonImpropers will be
    called after coordinates are loaded.

    swapProtons and assumeSegid are handled by the function
    matchInexactAtomEntry. 

    This function results an InitCoordsResult structure containing
    bfactors, occupancies, and remarks read from the PDB record.

    """
    fixMethylImpropers=kwargs.pop('fixMethylImpropers', False)
    fixProtonImpropers=kwargs.pop('fixProtonImpropers', False)

    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)
    
    xSim = getXplorSimulation(selection.simulation())
    command = xSim.command
    
    #don't erase restraints when adding atoms
    commandNoEcho(xSim, 'struct doScratch=false end',verbose=False)

    if verbose<0:
        verbose=0
        import simulationWorld
        simWorld = simulationWorld.SimulationWorld_world()
        if simWorld.logLevel() != 'none':
            verbose=1
            pass
        pass

    if erase:
        # FIX: this should be an option to PDBTool::read
        import vec3
        clearCoord = vec3.Vec3(9999.,9999.,9999.)
        for atom in selection:
            atom.setPos(clearCoord)
            pass
        pass
        
    

    if type(files)==type("string"):
        files=[files]
        pass

    if not files and not string:
        print("initCoords: warning no file or string was specified.")
    
    from pdbTool import PDBTool
    pdb = PDBTool("",selection,useChainID,strictResNames)
    pdb.setIncludeHETATM(includeHETATM)
    pdb.setVerbose(verbose)

    unknownEntries=[]
    for file in files:
        file,model = splitModel(file,model)
        pdb.setFilename(file)
        unknownEntries += pdb.read(model)
        pass

    if string:
        pdb.setFilename('')
        pdb.setContents(string)
        unknownEntries += pdb.read(model)
        pass

    if verbose>1:
        print("pdb file has format:", pdb.format())

    unreadEntries=0
    swapProtons=kwargs.pop('swapProtons', False)
    assumeSegid=kwargs.pop('assumeSegid', None)
    if len(kwargs):
        raise Exception("unexpected keyword argument(s): " + str(kwargs))
    for entry in unknownEntries:
        if not matchInexactAtomEntry(entry,pdb,selection,
                                     strictResNames,
                                     verbose=verbose,
                                     assumeSegid=assumeSegid,
                                     swapProtons=swapProtons):
            unreadEntries += 1
            pass
        pass

    if maxUnreadEntries != None:
        doException=False
        if type(maxUnreadEntries)==type(1):
            if unreadEntries>maxUnreadEntries:
                doException=True
                pass
            pass
        elif (unreadEntries>maxUnreadEntries[0] and
              unreadEntries>maxUnreadEntries[1]*len(selection)):
            doException=True
            pass
        if doException:
            raise Exception("too many unreadable ATOM entries: %d" %
                            unreadEntries)
        pass
    
    if verbose:
        if unreadEntries:
            print("initCoords: Warning: unable to read %d pdb ATOM entries" % \
                  unreadEntries, end=' ')
            print(" (nonpseudoatom)")
            pass
        unknown = AtomSel("not known")
        if len(unknown):
            print("initCoords: still %d unknown atomic coordinates" % \
                  len(unknown))
            pass
        pass

    if correctSymmetricSidechains:
        from selectTools import correctSymmetricSidechains
        correctSymmetricSidechains(sim=xSim)
        pass

    retVal = InitCoordsResult()
    retVal.bfactors = pdb.bfactors()
    retVal.occupancies = pdb.occupancies()
    retVal.remarks = pdb.remarks()

    if loadOccupanciesBfactors:
        b=retVal.bfactors
        q=retVal.occupancies

        import atomSelLang
        atomSelLang.addAttribute(selection.simulation(),'b',b)
        atomSelLang.addAttribute(selection.simulation(),'q',q)
        
        outputState=xSim.disableOutput()
        for atom in AtomSel("known"):
            index = atom.index()
            #can mess up if atom ordering is changed
            xSim.fastCommand("vector do (b=%f)  (id %d)" % (b[index],index+1))
            xSim.fastCommand("vector do (q=%f)  (id %d)" % (q[index],index+1))
            pass
        xSim.enableOutput(outputState)
        pass
            
    if biomt:
        initBiomtCoords(biomt,xSim)

    if deleteUnknownAtoms:
        unknown = AtomSel("not known",xSim)
        if verbose and len(unknown):
            print("initCoords:  Deleting %d atoms with unknown coordinates."  % len(unknown))
            pass

        retVal.bfactors    = syncArrayBeforeDelete(retVal.bfactors, unknown)
        retVal.occupancies = syncArrayBeforeDelete(retVal.occupancies, unknown)
        selection.simulation().deleteAtoms("not known")
        updatePseudoAtoms( selection.simulation() )
        pass

    if fixMethylImpropers:
        from selectTools import fixMethylImpropers
        fixMethylImpropers(verbose=True)
        pass
    if fixProtonImpropers:
        from selectTools import fixProtonImpropers
        fixProtonImpropers(verbose=True)
        pass


    return retVal


def syncArrayBeforeDelete(anArray, selectionToDelete):
    """ Given an array of values indexed by atom number, and a selection of
    atoms to be deleted, remove the entries in the array so that they will
    match up with the correct atoms after the atoms are actually deleted."""

    if len(anArray) != selectionToDelete.simulation().numAtoms():
        raise 'array size (%d) != number of atoms in simulation (%d)' % \
              (len(anArray), selectionToDelete.simulation().numAtoms())
    
    import cdsVector
    newV = cdsVector.CDSVector_double(len(anArray) - len(selectionToDelete))

    newc=0
    for c in range(len(anArray)):
        if not selectionToDelete.containsIndex(c):
            newV[newc] = anArray[c]
            newc += 1
            pass
        pass
    
    return newV

def writePDB(filename,
             writeChainID=False,
             selection="all",
             occupancies=None,
             bFactors=None,
             saveOccupanciesBfactors=False,
             remarks=[],
             ):
    """
    write atoms from specified selection to a PDB file with the specified
    filename

    If writeChainID is specified, the segid is stored in the chainID column.

    The selection argument should be an <m atomSel>.AtomSel object specifying
    atoms to write out.

    The arguments occupancies and bFactors, if specified, should be numerical
    sequences of length len(selection). These are placed in the occupancies and
    bFactor PDB fields.
    
    If saveOccupancyBfactors, occupancies and bfactors are read from the
    associated <m xplorSimulation>.XplorSimulation XPLOR process and placed in
    the appropriate pdb columns.

    remarks should be a sequence of strings without newlines, which are added
    as REMARKS statement entries to the PDB file.
    """
    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)
    

    from pdbTool import PDBTool
    pdb=PDBTool(filename,selection)
    pdb.setWriteChainID(writeChainID)

    if occupancies:
        if len(occupancies) != len(selection):
            raise Exception("occupancies should have length len(selection)")
            pass
        for i,atom in enumerate(selection):
            pdb.setAux1(atom,occupancies[i])
            pass
        pass

    if bFactors:
        if len(bFactors) != len(selection):
            raise Exception("bFactors should have length len(selection)")
            pass
        for i,atom in enumerate(selection):
            pdb.setAux2(atom,bFactors[i])
            pass
        pass


    if saveOccupanciesBfactors:
        xSim = getXplorSimulation(selection.simulation())
        outputState=xSim.disableOutput()
        for atom in selection:
            index=atom.index()
            (q,) = xSim.fastCommand("vector show elem (q) (id %d)" %
                                    (index+1), "result")
            (b,) = xSim.fastCommand("vector show elem (b) (id %d)" %
                                    (index+1), "result")
            pdb.setAux1(atom,q)
            pdb.setAux2(atom,b)
            pass
        xSim.enableOutput(outputState)
        pass
    for remark in remarks:
        pdb.addRemark(remark)
        pass
        
    pdb.write()
    return pdb.contents()

def writeCIF(filename,
             selection="all",
             entryName="XplorNIH-Structure",
             #writeChainID=False,
             #occupancies=None,
             #bFactors=None,
             remarks=[],
             ):
    """
    Write atoms from specified selection to an mmCIF file with the specified
    filename.

    The selection argument should be an <m atomSel>.AtomSel object specifying
    atoms to write out.

    No backup file is created in this version.
    """
    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)
    
    ofile=open(filename,"w")
    from cif import CifDatablock
    dataBlock = CifDatablock()
    
    from cifTools import genAtomSite
    from cif import Cif
    mmcif = Cif()
    mmcif.setUseTrailingPound( True )
    if remarks:
        if not type(remarks)==type(""):
            remarks = "\n".join(remarks)
            pass
        from cif import CifCategory
        cat = CifCategory()
        cat.addKey("text")
        cat.addValue("text", remarks)
        dataBlock["XplorNIH_remarks"] = cat
        pass
    dataBlock = genAtomSite(selection,dataBlock=dataBlock)
    mmcif[entryName] = dataBlock
    ofile.write( mmcif.asString() )
    return


def addPseudoResName(resname):
    """
    Add a residue name to the set which correspond to pseudo atoms.

    This function only updates the list of residue names which correspond
    to pseudoatoms. The function updatePseudoAtoms must be called after atoms
    are created to actually update the PSEUdo selection.
    """
    pseudoResNames.add(resname)
    return
    
def updatePseudoAtoms(sim=None):
    """
    Update the set of pseudo atoms selected by the PSEUdo atom selection
    keyword.
    """

    import simulation
    if not sim: sim=simulation.currentSimulation()
    xSim = getXplorSimulation(sim)

    outputState=xSim.disableOutput()
   
    xSim.syncTo()
    for resname in pseudoResNames:
        xSim.fastCommand("vector iden (pseudo) (pseudo or (resname %s))" %
                         resname)
        pass
    xSim.enableOutput(outputState)

    from atomSel import AtomSel, union

    sel=AtomSel("not all",sim)
    for resname in pseudoResNames:
        sel = union(sel,AtomSel("resname %s" % resname,sim))
        pass
    import atomSelLang
    atomSelLang.setNamedSelection(sim,"pseudo",sel.indices())

    if ( sim.type() == "EnsembleSimulation"):
        sim.sync()
        from ensembleSimulation import EnsembleSimulation_currentSimulation
        esim=EnsembleSimulation_currentSimulation()
        if esim==None:
            return
        for i in range(esim.size()):
            sim = esim.members(i)
            sel=AtomSel("not all",sim)
            for resname in pseudoResNames:
                sel = union(sel,AtomSel("resname %s" % resname,sim))
                pass
            atomSelLang.setNamedSelection(sim,"pseudo",sel.indices())
            pass
        pass
    
                           
    return

    
# residue names of pseudo atoms which may or may not be present
# in a coordinate record.
pseudoResNames=set()
import varTensorTools    #import modules to populate pseudoResNames
import diffPotTools
import prePotTools

def matchInexactAtomEntry(entry,pdbTool,selection,strictResNames,
                          assumeSegid=None,
                          swapProtons=False,
                          verbose=0):
    """ try to match a PDB record to an atom if the name is not exactly the
    same. The coordinates will only be overwritten if they are unknown.

    The current algorithm for matching is as follows:
      match H to HT1 
      else match if there's only a single undefined atom whose name starts
      with the same character.
      else
      match if atom names start with the same character and the remainder
      of the characters are simple tranpositions of each other.
      else
      match if there's only a single undefined atom whose name starts
      with and ends with the same character.
      else
      match if the first name of the name is a digit present in the correct
      atom name and the second character matches the first character of the
      correct name.
      else
      match if the two names have two matching characters.
      else
      match up   O --> OT1
               OXT --> OT2

      If the match is not unique, it is not made.

    The residue number and segment id must match.

    Entries with nonnumeric residue numbers are ignored.

    For atoms with residue names corresponding to those of pseudo atoms
    (e.g. ANI), this function prints a warning on missing atoms, but returns
    a True value.

    If strictResNames is True, then residue names in the PDB must match
    those in the PSF (default), otherwise, the residue names need not match
    and backbone atoms for mutants can be read. Note that sidechain atoms (or
    base atoms in nucleic acids) may be strangely placed if this option is
    False.

    If the argument swapProtons==True, this function has the following
    side effect: for gylcine, which should have HB1 and HB2, but
    the input PDB names them HB2 and HB3, stereo-assignment is maintained by
    swapping the positions. swapProtons defauls to False.

    The argument assumeSegid affects matchInexactAtomEntry's behavior when
    the segid field is blank. If the PSF has a blank segid, and the PDB has
    assuemSegid, the atom will match. Likewise, if the PDB has a blank segid
    (and chainID), the entry will match if the PSF has assumeSegid.
    """
    

    from pdbFormat import getAtomEntry
    name       = getAtomEntry('name',entry)
    altLoc     = getAtomEntry('altLoc',entry)
    resName    = getAtomEntry('resName',entry)
    chainID    = getAtomEntry('chainID',entry)
    resSeq     = getAtomEntry('resSeq',entry)
    iCode      = getAtomEntry('iCode',entry)
    try:
        x  	       = float(getAtomEntry('x',entry))
        y  	       = float(getAtomEntry('y',entry))
        z  	       = float(getAtomEntry('z',entry))
    except ValueError:
        if verbose: print('error reading coordinate values:' ,entry)
        return 1
        
    occupancy  = entry[ 54:60]
    tempFactor = entry[ 60:66]
    segID      = entry[ 72:76].strip() #lead/trail whitespace not significant
    element    = entry[ 76:78]
    charge     = entry[ 78:80]

    if pdbTool.useChainID():
        if chainID!=' ': segID = chainID+'   '
        pass

    readAtoms = pdbTool.readAtoms()

    if assumeSegid!=None and assumeSegid==segID: segID=""
    if assumeSegid!=None and segID.strip()=="": segID=assumeSegid
    
    name = name.strip()
    from atomSel import AtomSel, intersection
    resAtoms = intersection(selection,
                            AtomSel('not known and (resid %s and segid "%s")' %
                                    (resSeq,segID),
                                    selection.simulation()))
    resAtoms = [a for a in resAtoms if not readAtoms[a.index()]]
    matchAtom=None
    swapAtom=None
    rule=-1

    from psfGen import renameResidues
    resName = renameResidues([resName])[0]
    if strictResNames and resAtoms:
        recResName = resAtoms[0].residueName()
        if recResName != resName:
            from psfGen import variantResidues
            found=False
            for variant in sum(list(variantResidues.values()),[]):
                if variant.name==resName and variant.baseName==recResName:
                    found=True
                    break
                pass
            if not found:
                raise Exception("residue name mismatch %s != %s\n" % (resName,
                                                                      recResName)+
                                "  entry:" + entry)
            pass
        pass

    if not matchAtom:
        hn = [a for a in resAtoms if a.atomName()=="HN"]
        if hn and name=='H':
            matchAtom = hn[0]
            rule=1
            pass
        pass    
    if not matchAtom:
        ht1 = [a for a in resAtoms if a.atomName()=="HT1"]
        if ht1 and (name=='H1' or name=='H2' or name=='H'):
            matchAtom = ht1[0]
            rule=2
            pass
        pass    
    if not matchAtom:
        ht2 = [a for a in resAtoms if a.atomName()=="HT2"]
        if ht2 and name=='H3':
            matchAtom = ht2[0]
            rule=3
            pass
        pass    

    #phosphate oxygens 
    if not matchAtom:
        op1 = [a for a in resAtoms if a.atomName()=="O1P"]
        if name=='OP1' and op1 and not readAtoms[op1[0].index()]:
            matchAtom = op1[0]
            rule=4
            pass
        pass
    if not matchAtom:
        op2 = [a for a in resAtoms if a.atomName()=="O2P"]
        if name=='OP2' and op2 and not readAtoms[op2[0].index()]:
            matchAtom = op2[0]
            rule=5
            pass
        pass

    if (not matchAtom and resName in
       "ARG ASN ASP CYS GLN GLU GLY ILE LEU LYS MET PHE PRO SER TR{".split()):
        if name[0]=="H" and name[-1]=="3":
            for atom in [atom for atom in resAtoms if atom.atomName()[0]=="H"]:
                if len(name)!=len(atom.atomName()):
                    continue
                if name[:2] != atom.atomName()[:2]:
                    continue
                if atom.atomName()[-1]=="1":
                    matchAtom=atom
                    rule=6
                    swapAtomName=name[:-1]+'2'
                    try:
                        swapAtom = AtomSel('''known and name %s and resid %d
                                           and segid "%s"''' %
                                           (swapAtomName,atom.residueNum(),
                                            atom.segmentName()),
                                           selection.simulation())[0]
                    except:
                        pass
                    break
                pass
            pass
        pass
        
    if not matchAtom:
        # match atom name that begin with same character
        for atom in resAtoms:
            if len(name)>1 and atom.atomName().startswith(name[:-1]):
                #check for names with transposed characters
                if set(atom.atomName()) == set(name):
                    matchAtom = atom
                    rule=7
                    break
                elif matchAtom:   #if two match, then none match
                    matchAtom = None
                    rule=8
                    break
                else:
                    matchAtom = atom
                    rule=9
                    pass
                pass
            pass
        pass
    if not matchAtom:
        # match atom name that begin and end with same character
        atom = [a for a in resAtoms if (a.atomName()[0] ==name[0] and
                                 a.atomName()[-1]==name[0-1] )]
        if len(atom)==1:
            matchAtom = atom[0]
            rule=10
            pass
        pass
    if not matchAtom:
        for atom in resAtoms:
            if ((name[0].isdigit() and atom.atomName().find(name[0])>=0) and
                name[1:] == atom.atomName()[:-1]):
                if matchAtom:
                    matchAtom = None
                    break
                else:
                    matchAtom = atom
                    rule=11
                    pass
                pass
            pass
        pass
    if not matchAtom:
        # match if 2 or more characters are common
        for atom in resAtoms:
            count=0
            atomName = atom.atomName()
            tname=str(name)
            for c in atomName:
                tindex=tname.find(c)
                if tindex>=0:
                    count += 1
                    tname=tname[:tindex]+tname[tindex+1:]
                    pass
                pass
            if count>=2:
                if matchAtom:
                    matchAtom = None
                    break
                else:
                    matchAtom = atom
                    rule=12
                    pass
                pass
            pass
        pass
    #rules for strange nucleic base HN atom naming
    if not matchAtom and name=="HN'":
        h21 = [a for a in resAtoms if a.atomName()=="H21"]
        h41 = [a for a in resAtoms if a.atomName()=="H41"]
        h61 = [a for a in resAtoms if a.atomName()=="H61"]
        if h21 and not readAtoms[h21[0].index()]:
            matchAtom = h21[0]
            rule=13
        elif h41 and not readAtoms[h41[0].index()]:
            matchAtom = h41[0]
            rule=14
        elif h61 and not readAtoms[h61[0].index()]:
            matchAtom = h61[0]
            rule=15
            pass
        pass
    if not matchAtom and name=="HN''":
        h22 = [a for a in resAtoms if a.atomName()=="H22"]
        h42 = [a for a in resAtoms if a.atomName()=="H42"]
        h62 = [a for a in resAtoms if a.atomName()=="H62"]
        if h22 and not readAtoms[h22[0].index()]:
            matchAtom = h22[0]
            rule=16
        elif h42 and not readAtoms[h42[0].index()]:
            matchAtom = h42[0]
            rule=17
        elif h62 and not readAtoms[h62[0].index()]:
            matchAtom = h62[0]
            rule=18
            pass
        pass
        

    #special rules for the O-terminus
    if not matchAtom:
        ot1 = [a for a in resAtoms if a.atomName()=="OT1"]
        if name=='O' and ot1 and not readAtoms[ot1[0].index()]:
            matchAtom = ot1[0]
            rule=19
            pass
        pass
    if not matchAtom:
        ot2 = [a for a in resAtoms if a.atomName()=="OT2"]
        if name=='OXT' and ot2 and not readAtoms[ot2[0].index()]:
            matchAtom = ot2[0]
            rule=20
            pass
        pass

    if not matchAtom and resName=='WAT':
        oAtom = [a for a in resAtoms if a.atomName()=="OH2"]
        if name=="O" and oAtom and not readAtoms[oAtom[0].index()]:
            matchAtom=oAtom[0]
            rule=21
            pass
        pass

    if not altLoc in pdbTool.allowedAltLoc():
        matchAtom=False

    if matchAtom:
        if verbose:
            print("matchInexactAtomEntry:", \
                  "matching entry %s %s %s to atom %s %d %s" % \
                  (segID,resSeq,name,
                   matchAtom.segmentName(),matchAtom.residueNum(),
                   matchAtom.atomName()), end=' ')
            print(" [rule %d]" % rule)
            pass
                                                                
        matchAtom.setPos( [x,y,z] )
        try:
            occupancy = float(occupancy)
        except ValueError:
            occupancy = 0
            pass
        try:
            tempFactor = float(tempFactor)
        except ValueError:
            tempFactor = 0
            pass
        pdbTool.setAux1(matchAtom, occupancy)
        pdbTool.setAux2(matchAtom, tempFactor)

        if swapProtons and swapAtom:
            if verbose:
                print("matchInexactAtomEntry:", \
                      "swapping entries %s %s %s and " \
                      "%s %d %s" % \
                      (segID,resSeq,name,
                       swapAtom.segmentName(),swapAtom.residueNum(),
                       swapAtom.atomName()))
                
            tmp=swapAtom.pos()
            swapAtom.setPos(matchAtom.pos())
            matchAtom.setPos(tmp)
            tmp = pdbTool.aux1(swapAtom)
            pdbTool.setAux1(swapAtom, pdbTool.aux1(matchAtom) )
            pdbTool.setAux1(matchAtom, tmp)
            tmp = pdbTool.aux2(swapAtom)
            pdbTool.setAux2(swapAtom, pdbTool.aux2(matchAtom) )
            pdbTool.setAux2(matchAtom, tmp)
            pass            
        pass
    else:
        if verbose:
            print("matchInexactAtomEntry:", \
                  "found no match for entry %s %s %s" % (segID,resSeq,name))
            pass
#            print "matchInexactAtomEntry:", \
#                  "found no match for entry %s %s %s %s" % (segID,resSeq,resName,
#                                                            name)
        pass

    if resName in pseudoResNames:
        if verbose: print(" pseudo atom")
        matchAtom=True
        pass        
    
    return matchAtom != None

def addDisulfideBond(sel1,sel2):
    """
    deprecated. Calls psfGen.addDisulfideBond
    """
    import psfGen
    return psfGen.addDisulfideBond(sel1,sel2)
                  
def initNBond(cutnb=4.5,
              rcon=4.0,
              nbxmod=3,
              selStr="all",
              tolerance=0.5,
              repel=None,
              onlyCA=0,
              simulation=0,
              suppressException=False,
              ):
    """standard initialization of the XPLOR VDW term, configured to use the
    repel potential. The XPLOR nonbonded potential term is described
    <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/node117.html here>.

    Note that cutnb should be greater than rmax + 2*tolerance, where rmax
    is the largest vdw radius.

    selStr specifies which atoms to use in nonbonded energy calculations.
    
    If onlyCA is True, this string gets set to "name CA" - but this option is
    deprecated.

    The default value of repel is 0.8 for pre-3.1 version Xplor-NIH protein
    and nucleic acid parameters, and 0.9 for more later versions. If a repel
    value less than 0.9 specified to be used with post-3.0 parameter sets, an
    exception will be raised, unless suppressExceptions is set to True.
    """

    from repelPotTools import checkGetRepel
    repel = checkGetRepel(repel,suppressException)

    noSelStr="pseudo" #exclude these nonbonded interactions

    if onlyCA: selStr="name CA"

    xSim = getXplorSimulation(simulation)

    xSim.fastCommand("""     
          constraints
            interaction (%s and (not (%s))) (%s and (not (%s)))
            weights * 1 vdw 1 end 
            interaction  (not (%s) and (not (%s))) (not (%s))
            weights * 1 vdw 0 end 
          end""" % (selStr,noSelStr,selStr,noSelStr,selStr,noSelStr,noSelStr) )

        
    if repel==None: repel=0.8
    xSim.fastCommand("""
    parameters
    nbonds
    atom
    nbxmod %d
    wmin  =   0.01  ! warning off
    cutnb =   %f   ! nonbonded cutoff
    tolerance %f
    repel=    %f   ! scale factor for vdW radii = 1 ( L-J radii)
    rexp   =  2     ! exponents in (r^irex - R0^irex)^rexp
    irex   =  2
    rcon=%f      ! actually set the vdW weight
    end
    end 
    """ % (nbxmod,cutnb,tolerance,repel,rcon) )

    return
    

def initRamaDatabase(type='protein',
                     simulation=0,
                     selection="all"):
    """Standard initialization of the <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/rama-db.html
    RAMA torsion-angle database potential>.

    The argument type (a string) is either 'protein' or 'nucleic', specifying
    the system type (protein or nucleic acid, respectively).

    The selection argument (a string) can be used to specify a
    subregion of the molecular system in the case that one does not
    wish the RAMA database to apply to all residues. This option is
    only implemented for nucleic acids.
    """
    xSim = getXplorSimulation(simulation)
    outputState=xSim.disableOutput()
    from selectTools import numResidues
    updatePseudoAtoms(xSim)
    global nRama
    nRama=45*numResidues( AtomSel('not PSEUDO',xSim) )
    print(' allocating space for %d RAMA restraints' % nRama)
    xSim.command("""
    eval ($krama=1.)
    rama
    nres=%d
    end
    """% nRama)
    if type=='protein':
        xSim.command("""
        rama
        @QUARTS:2D_quarts_new.tbl
        @QUARTS:3D_quarts_new.tbl
        @QUARTS:forces_torsion_prot_quarts_intra.tbl
        end
        @QUARTS:setup_quarts_torsions_intra_2D3D.tbl
        """)
        xSim.enableOutput(outputState)
    elif type=='nucleic':
        xSim.command("""
        evaluate ($knuc=1.0)
        rama
        @QUARTS:nucleic_deltor_quarts2d.tbl
        @QUARTS:nucleic_deltor_quarts3d.tbl
        @QUARTS:nucleic_deltor_quarts4d.tbl
        @QUARTS:force_nucleic_quarts2d.tbl
        @QUARTS:force_nucleic_quarts3d.tbl
        @QUARTS:force_nucleic_quarts4d.tbl
        end
        vector identify (store6) (%s)
        @QUARTS:setup_nucleic_2d3d.tbl
        @QUARTS:setup_nucleic_4d.tbl
        """% selection)
        xSim.enableOutput(outputState)
    else:
        xSim.enableOutput(outputState)
        raise Exception("initRamaDatabase: unknown database type: "+ type)
    return


def initOrie(system='dna', pairs=[], selection='all'):
    """Set up the <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/posdb.html
    nucleobase-nucleobase positional knowledge-based potential>.

    The argument system is either 'dna' or 'rna'.  pairs is a sequence
    of (x, y) tuples, where x and y are the residue numbers of two
    bases involved in a Watson-Crick (WC) interaction.  Alternatively,
    pairs may be the name of a file specifying the residue numbers in
    WC pairs; the format of the file should be, e.g.,

    1  24
    2  23

    where residue 1 is paired with 24, and 2 with 23. The use of pairs requires
    that the specified residue numbers be unique. i.e. not replicated in
    separate segids.

    The selection argument can be used to specify a subset of the system to
    apply this term to, but its selection string is not consulted if the
    pairs argument is specified.

    """
    # Written by Guillermo A. Bermejo.

    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)
    
    def setup(system):
        from selectTools import numResidues
        numResidues = numResidues()
        numRestraints = numResidues * 500 #525  a guess (450 is too small)
        numResidues += 1 #seems to be an off-by-one bug in the Fortran code
        cmd = '''
            orient
               nres %d
               residues %d
               maxgauss 1000
               scale 1
               shape quart
               '''  % (numRestraints,numResidues)

        if system == 'dna':
            cmd += '''
               !based on structures 2.3 A or less
               @DNA_DNA_PAIRS:dna_dna_quarts_combo_res23.tbl  
            end

            @DNA_DNA_PAIRS:dna_gaussians_bases_setup.tbl
            '''
        elif system == 'rna':
            cmd += '''
               @RNA_RNA_PAIRS:rna_quarts_combo_res3.tbl 
               @RNA_RNA_PAIRS:rna_nonseq_delpos_quarts_res3.tbl
            end
            '''
            cmd += '@RNA_RNA_PAIRS:rna_quarts_setup.tbl\n'
            cmd += '@RNA_RNA_PAIRS:rna_nonconsec_setup.tbl\n'
        else:
            raise Exception("system must be either 'dna' or 'rna'")
        xSim.command(cmd)
        return

    if type(pairs) is str:  # assume pairs is a filename
        infile = open(pairs, 'r')
        contents = infile.read()
        infile.close()
        lines = contents.splitlines()
        pairs = []
        for line in lines:
            pair = line.split()
            pairs.append((int(pair[0]), int(pair[1])))

    xSim = getXplorSimulation(selection.simulation())
    outputState=xSim.disableOutput()
    updatePseudoAtoms(xSim)  # DO I NEED THIS?

    import atomSel
    segids = {}

    xSim.command('vector do (store4 = 0) (all)')        
    
    if pairs:

        # Assume first item in each pair belongs to one strand, 
        # the second to the other.
        for (x, y) in pairs:

            # Save original segids to restore them later on.
            # (This assumes resid is enough to determine the residue (i.e., 
            # there are no overlapping residue numbers between the chains).)
            segids[x] = atomSel.AtomSel('resid %i' % x)[0].segmentName()
            segids[y] = atomSel.AtomSel('resid %i' % y)[0].segmentName() 
            
            xSim.command('''
            vector do (segid = "NUC1") (resid %i)
            vector do (segid = "NUC2") (resid %i)
            ''' % (x, y))

            xSim.command('vector do (store4 = %i) (resid %i)' % (x, y))
            xSim.command('vector do (store4 = %i) (resid %i)' % (y, x))

        xSim.command('''       
        vector do (ustring = "") (all)
        vector do (ustring = "NUC2") (segid NUC1)
        vector do (ustring = "NUC1") (segid NUC2)
        ''')

        xSim.command('vector identify (store9) (segid NUC1 or segid NUC2)')
        
        setup(system) # final setup

        # Restore original segids.
        for (x, y) in pairs:
            xSim.command('''
            vector do (segid = "%s") (resid %i)
            vector do (segid = "%s") (resid %i)
            ''' % (segids[x], x, segids[y], y))

    elif not pairs:

        xSim.command('vector identify (store9) (%s)' % selection.string())

        setup(system)  # final setup

    xSim.enableOutput(outputState)

    return
            


def initDihedrals(filenames=[],
                  string="",
                  scale=1,
                  useDefaults=True,
                  reload=False,
                  simulation=0):
    """Initialize the XPLOR <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/cdih.html
    dihedral restraints (CDIH) potential term>.

    Parameters are:
       filenames   - either a single filename, or a sequence of filenames of
                     dihedral restraint assignment tables.
       string      - assignment table as a plain string.
       scale       - scale factor (defaults to 1).
       useDefaults - use the default sidechain restraints (default: True)
                     these force chi2 angles of PHE, TYR, ASP, and chi3 of GLU
                     to obey IUPAC naming conventions.
       reload      - if True, reload the previously-loaded dihedral restraints.
                     This is useful if the restraints are lost do to the
                     number of atoms changing (SCRATCHing).
    """

    xSim = getXplorSimulation(simulation)

    global prev_filenames, prev_string, prev_useDefaults
    if reload and 'prev_filenames' in globals():
        filenames=prev_filenames
        string=prev_string
        useDefaults=prev_useDefaults
    else:
        prev_filenames=  filenames
        prev_string=     string
        prev_useDefaults=useDefaults
#        reload=False
        pass
        

    outputState=xSim.disableOutput()
    if not reload:
        #FIX: this should be determined properly
        xSim.fastCommand("set display=none end")
        xSim.fastCommand("restraints dihed nass = 10000 end")  
        xSim.fastCommand("set display=OUTPUT end")
        pass
    xSim.command("""
    restraints dihed 
    reset
    scale %f
    end""" % scale)
    if type(filenames)==type("string"): filenames = [filenames]
    from xplorPot import XplorPot
    pot=XplorPot('CDIH',xSim)
    restraints0=pot.numRestraints()
    pot.modified.set()
    for file in filenames:
        if not file: continue
        xSim.fastCommand("restraints dihed @%s end" % file)
        pass
    if string:
        xSim.fastCommand("restraints dihed %s end" % string)
    if useDefaults:
        xSim.fastCommand("""
!
! phe_angles.tbl
!
! constrains the chi2 angles of phe, tyr, asp, and glu
! in a protein to -90..90. This so that the dihedral angles follow
! IUPAC naming convention Biochemistry 9, 3471 (1970).
!
! JJK 3/10/97
!

for $res in id (name cd1 and (resn phe or resn tyr)) loop ang

   restraints dihedral 
      assign 
         (byresidue id $res and name ca)
         (byresidue id $res and name cb)
         (byresidue id $res and name cg)
         (byresidue id $res and name cd1)
         1.0 90.0 90.0 2
   end

end loop ang

for $res in id (name od1 and resn asp) loop ang

   restraints dihedral 
      assign 
         (byresidue id $res and name ca)
         (byresidue id $res and name cb)
         (byresidue id $res and name cg)
         (byresidue id $res and name od1)
         1.0 0.0 90.0 2
   end

end loop ang

for $res in id (name oe1 and resn glu) loop ang

   restraints dihedral 
      assign 
         (byresidue id $res and name cb)
         (byresidue id $res and name cg)
         (byresidue id $res and name cd)
         (byresidue id $res and name oe1)
         1.0 0.0 90.0 2
   end

end loop ang
""")
        pass
    xSim.enableOutput(outputState)
    numRestraints=pot.numRestraints() - restraints0
    from simulationWorld import world
    xSim.sync()
    if not reload and world().logLevel()!="none":
        print("  read %d dihedral restraints." % numRestraints)
    return


def initPlanarity(filenames=[],
                  string="",
                  basePairs=[],
                  simulation=0):
    """Initialize the XPLOR <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/plan.html
    planarity restraint (PLAN) potential term>.

    Parameters are:
        filenames - either a single filename (a string), or a sequence of such
                    filenames of planarity restraint tables.
        string    - table as a plain string.
        bairPairs - a sequence of pairs of residues for which to
                    generate planarity restraints. Each base pair can be
                    represented as a pair of residue numbers, or as a pair of
                    tuples specifying (segid,resid).

    The restraint table(s) should only contain the <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/plan_syntax.html
    restraints-planar-statement>.
    Each statement is of the type:

        group
           selection=<selection>
           weight=<real>
        end

    which adds a new atom group to the planar restraints database, where
    <selection> is an XPLOR selection involving the atoms (at least four)
    restrained to form a plane, and <real> is the restraint's weight.

    Note that this function will clear any previously set up planarity
    restraints.

    """
    # Written by Guillermo A. Bermejo.
    
    xSim = getXplorSimulation(simulation)

    outputState = xSim.disableOutput()

    xSim.command("restraints plan initialize end")  # clear previous restraints

    if type(filenames) == str: filenames = [filenames]
    from xplorPot import XplorPot
    pot = XplorPot('PLAN', xSim)
    pot.modified.set()

    if basePairs:
        from selectTools import genPlanarityRestraints
        string += genPlanarityRestraints(basePairs,xSim)
        pass

    for filename in filenames:
        xSim.command("restraints plan @%s end" % filename)
    if string:
        xSim.command("restraints plan %s end" % string)    

    xSim.enableOutput(outputState)

##    from simulationWorld import world
##    if world().logLevel() != "none":
##        print "  read %i groups for planarity restraints." % groups

    return


# We might want to have the following function:
##def genPlanarityGroupStatement(pairs):
##    """Return a string with planarity restraint group statements from pairs.
##
##    This function helps set up the XPLOR <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/plan.html planarity restraint (PLAN) potential
##    term> by generating <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/plan_syntax.html restraints-planar-statement> from input residue pairs.
##
##    """



def initCollapse(sel    ="all",
                 Rtarget=None,
                 scale  =1):
    """Initialize the XPLOR <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/node412.html
    radius of gyration potential term>.

    Parameters are:
    
       sel - Selection string or atomSel.AtomSel instance specifying atoms to
             include in the calculation of the radius of gyration (Rgyr)
             (default: all atoms).
             
       Rtarget - Target Rgyr.  In practice the target Rgyr should be smaller
                 than the expected value to force compaction.  If Rtarget is
                 not specified, it is predicted from the number of residues, N,
                 by
                    Rtarget = 2.2 * N^0.38 - 1   (in Angstroms)
                 (i.e., the expected Rgyr minus one).
                 
       scale - Scale factor (defaults to 1).  Note that the per-assignment scale
               is always 100, so the energy is actually scaled by 100*scale.
    """

    if type(sel)==type("string"): sel = AtomSel(sel)

    xSim = getXplorSimulation( sel.simulation() )

    from selectTools import numResidues
    numResidues = numResidues(sel)
    if not Rtarget:
        Rtarget = (2.2 * numResidues**0.38 -1)
        pass
    
    xSim.command("""
    collapse
    assign (%s) 100.0 %f
    scale %f
    end""" % (sel.string(), Rtarget, scale))
    return

def initCarb13(filenames=[],
               scale=0.5):
    """Initialize the XPLOR <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/node384.html
    carbon chemical shift potential term>.

    arguments are:
       filenames   - either a single filename, or a sequence of filenames of
                     chemical shift assignment tables.

    In addition to these files C13SHIFTS:rcoil_c13.tbl and
    C13SHIFTS:expected_edited.tbl are always added.
    """

    xSim = getXplorSimulation()

    if type(filenames)==type("string"): filenames = [filenames]
    restraints=""
    nres=0
    for file in filenames:
        for line in open(file).readlines():
            if line.strip().lower().startswith('assi'): nres += 1
            restraints += line
            pass
        pass

    outputState=xSim.disableOutput()
    xSim.command("""
    carbon
      phistep=180
      psistep=180
      nres=%d
      class all
      force %f
      potential harmonic
      set echo=off mess=off end 
      @C13SHIFTS:rcoil_c13.tbl             !rcoil shifts
      @C13SHIFTS:expected_edited_c13.tbl   !13C shift database
      
      %s
      
      set echo=$prev_echo mess=$prev_messages end 
    end""" % (nres,scale,restraints))
    xSim.enableOutput(outputState)
    return

def initHBDA(filename,
             scale=500,
             simulation=0):
    """Initialize the XPLOR hbda potential term
    """
    xSim = getXplorSimulation(simulation)
    lines=open(filename).readlines()
    nres = len([x for x in lines if x.strip().lower().startswith('assi')])
    outputState=xSim.disableOutput()
    xSim.command("""
    hbda
    nres %d
    class back
    @%s
    force %f
    end
    """ % (nres,filename,scale))
    xSim.enableOutput(outputState)
    return

def initHBDB(selection="not pseudo",
             restraintFile=None,
             prnfrq=0):
    """Initialize the XPLOR 
    <l http://nmr.cit.nih.gov/xplor-nih/xplorMan/hbdb.html hbdb>
    potential term database hydrogen-bonding term - h-bonds are determined
    dynamically. 

    The selection argument specifies which atoms in included in the
    calculation. Unfortunately, the current implementation allows only contiguous
    ranges of residue numbers. So the selection specifies only the segment name
    and the minimum and maximum residue number.

    prnfrq can be changed from zero to periodically print out HBDB info during
    dynamics and minimization.

    By default, HBDB runs in ``free'' mode, in which hydrogen bonds are
    determined on the fly. If you instead wish to use a fixed list of
    h-bonds, specify a filename for the restraintFile argument. 

    The HBDB term may have issues with hetereo multimers. Please use with
    caution.
    
    """

    from selectTools import minResid, maxResid, getSegids, convertToAtomSel
    selection = convertToAtomSel(selection)
    segids = getSegids(selection)

    simulation = selection.simulation()
    xSim = getXplorSimulation(simulation)
    outputState=xSim.disableOutput()
    cmd="""
      hbdb
       kdir = 0.20   !force constant for directional term
       klin = 0.08   !force constant for linear term (ca. Nico's hbda)
       nseg = %d      ! number of segments that hbdb term is active on
       """ % len(segids)
    from atomSel import AtomSel, intersection
    for segid in segids:
        minRes = minResid(intersection(selection,
                                         AtomSel('segid "%4s"' % segid,
                                                 simulation)))
        maxRes = maxResid(intersection(selection,
                                         AtomSel('segid "%4s"' % segid,
                                                 simulation)))
        print("segid: %s   minResid: %d  maxResid: %d" % (segid,minRes,maxRes))
        cmd += '''
        nmin = %d      !range of residues (1st)
        nmax = %d        !range of residues (last)
        segm = "%s"
        ''' % (minRes, maxRes, segid)
        pass
    cmd += """
       ohcut   =  2.60  !cut-off for detection of h-bonds
       coh1cut = 100.0  !cut-off for c-o-h angle in 3-10 helix
       coh2cut = 100.0  !cut-off for c-o-h angle for everything else
       ohncut  = 100.0  !cut-off for o-h-n angle
       updfrq = 10     !update frequency usually 1000
       prnfrq = %d     !print frequency usually 1000
       freemode  = 1     !mode= 1 free search
       fixedmode = 0     !if you want a fixed list, set fixedmode=1, and freemode=0
       mfdir = 0 ! flag that drives HB's to the minimum of the directional potential
       mflin = 0 ! flag that drives HB's to the minimum of the linearity potential
       kmfd = 10.0 ! corresp force const
       kmfl = 10.0 ! corresp force const
       renf = 2.30 ! forces all found HB's below 2.3 A
       kenf = 30.0 ! corresponding force const
       @HBDB:hbdb_files.inp
      end
      """ % prnfrq
    xSim.command(cmd)
    if restraintFile!=None:
        xSim.command(r'''hbdb
                      freemode = 0
                      fixedmode = 1
                      @%s
                     end''' % restraintFile)
        pass
    xSim.enableOutput(outputState)
    return


def initMinimize(ivm,
                 potList=0,
                 printInterval=10,
                 numSteps=500,
                 maxCalls=20000,
                 dEPred=0.001):
    """Initialize an <m ivm>.IVM object for Powell minimization.

    In addition to the function arguments, the following IVM parameters are
    initialized:
      constrainLengths
      maxDeltaE
      eTolerance
      gTolerance
    """
    ivm.resetReuse()
#    ivm.setConstrainLengths(0)
    ivm.setMaxDeltaE( 10000 )
    ivm.setStepType("powell")
    ivm.setNumSteps( numSteps )
    ivm.setMaxCalls( maxCalls )
    ivm.setPrintInterval( printInterval )
    ivm.setETolerance(1e-7)
    ivm.setGTolerance(1e-8)
    ivm.setDEpred(dEPred)
    if potList: ivm.setPotList( potList )
    return

def initDynamics(ivm,
                 bathTemp=-1,
                 finalTime=0.2,
                 numSteps=0,
                 stepsize=0.001,
                 potList=0,
                 printInterval=50,
                 initVelocities=0,
                 eTol_factor=0.001,
                 eTol_minimum=0
                 ):
    """Initialize an <m ivm>.IVM object for PC6 (6th order predictor-corrector)
    dynamics.

    The default value of bathTemp is the ivm's current value.

    In addition to the function arguments, the following IVM parameters are
    initialized:
      maxDeltaE       -> 10000
      responseTime    -> 5
      adjustStepsize  -> True
      scaleVel        -> True
      resetCMInterval -> 10

    Additionally, if initVelocities!=0, the initial velocities will be
    randomized.

    eTolerance is set by the following formula:
       eTolerance = eTol_factor * bathTemp + eTol_minimum
    """
    from atomAction import randomizeVelocities
    ivm.resetReuse()
    eTol_temp=0
    if bathTemp<0: bathTemp=ivm.bathTemp()
    
    ivm.setBathTemp(bathTemp)
    if initVelocities: randomizeVelocities( bathTemp,
                                            AtomSel("known",ivm.simulation) )
    eTol_temp=eTol_factor*ivm.bathTemp()

#    ivm.setConstrainLengths(0)
    ivm.setMaxDeltaE( 10000 )
    ivm.setStepType("pc6")
    ivm.setResponseTime(5)
    ivm.setStepsize( stepsize ) #initial stepsize value
    ivm.setETolerance( eTol_temp + eTol_minimum)
    ivm.setPrintInterval( printInterval )
    ivm.setAdjustStepsize(1)
    ivm.setScaleVel(1)
    ivm.setResetCMInterval( 10 )
    ivm.setFinalTime(finalTime)
    ivm.setNumSteps(numSteps)
    if potList: ivm.setPotList( potList )
    return

def massSetup(atomicMass=100,
              friction=10):
    """
    set uniform mass and friction values for all atoms.
    Pseudo atoms may have nonuniform masses.
    """
    from atomAction import SetProperty
    AtomSel("all").apply( SetProperty("mass",atomicMass) )
    #TODO: register these functions when loading the respective modules
    # also need relaxRatio in in this function
    import varTensorTools, prePotTools, planeDistTools, ensWeightsTools
    import relaxRatioPotTools
    varTensorTools.massSetup()
    prePotTools.massSetup()
    planeDistTools.massSetup()
    ensWeightsTools.massSetup()
    relaxRatioPotTools.massSetup()
    AtomSel("all").apply( SetProperty("fric",friction) )
    return

topologyFuncs=[]
def addTopologySetup(topologyFunc,name):
    """ register a topology setup function to be called by the
    torsionTopology and cartesianTopology helper functions. The second
    argument should identify the module in which the topology setup
    function is defined.
    """
    topologyFuncs.append( (topologyFunc,name) )
    return


    

def torsionTopology(ivm, fixedOmega=False, oTensors=[],
                    breakRiboseSel="name C4' or name O4'",
                    breakProlineSel='name CD or name CG',
                    flexRiboseRing="[nucleic]"):
    """Configure the <m ivm>.IVM topology for standard torsion angle setup.

    Given input ivm (an <m ivm>.IVM instance), this function groups
    rigid side chains and "breaks" disulfide bonds and proline and
    ribose rings.

    The bond broken in proline and ribose rings is respectively
    specified by arguments breakRiboseSel and breakProlineSel, each an
    XPLOR-style selection string selecting the atoms in the bond.  The
    lengths of broken bonds are to be restrained by the force field as
    they do not take on fixed values like those of "unbroken" bonds.
    Bond angles and improper dihedrals associated with the atoms of
    broken bonds are also flexible, and need to be restrained by the
    force field.  Additionally, all endocyclic bond angles within
    ribose rings can be made flexible by specifying an XPLOR-style
    selection string in flexRiboseRing (default: empty string).  Such
    selection should involve at least one atom from the targeted
    residues, as in the following examples:

    'not pseudo and tag'  -> one atom per residue, excluding pseudoatoms
    'tag and resid 1:33'  -> one atom per residue in the specified residue
                             range
    'tag and segid "A"'   -> one atom per residue in the specified segid

    (Selection of atoms within residues that do not have a ribose ring will
    cause the program to fail.)

    If fixedOmega is set to True, also fix protein omega backbone angles.

    oTensors is a list of <m varTensor>.VarTensor objects used for alignment
    tensors.

    Warning: This function sets all groupings and hinge types which are not
    already specified, so it must be called last in the setup of an IVM's
    topology.
    
    """
    # WARNING: default values of breakRiboseSel and breakProlineSel should
    # match that of breakSelStr in selectTools.IVM_breakProlines and
    # selectTools.IVM_breakRiboses, respectively.

    for action in ivm.configActions:
        action.runAction()
        pass

    for topologyFunc,name in topologyFuncs:
        topologyFunc(ivm)
        pass
    
    import varTensorTools, prePotTools, planeDistTools, diffPotTools
    import relaxRatioPotTools, ensWeightsTools
    varTensorTools.topologySetup(ivm,oTensors)
    prePotTools.topologySetup(ivm)
    planeDistTools.topologySetup(ivm)
    diffPotTools.topologySetup(ivm)
    relaxRatioPotTools.topologySetup(ivm)
    ensWeightsTools.topologySetup(ivm)
    
    import selectTools
    if fixedOmega: selectTools.IVM_groupRigidBackbone(ivm)
    selectTools.IVM_groupRigidSidechain(ivm)
    selectTools.IVM_breakProlines(ivm, breakSelStr=breakProlineSel)
    selectTools.IVM_breakRiboses(ivm, breakSelStr=breakRiboseSel)
    selectTools.IVM_breakDisulfides(ivm)
    if flexRiboseRing:
        import psfGen
        psfGen.initResidueNames() # for [nucleic] abbreviation
        selectTools.IVM_flexibilizeRiboses(ivm, flexRiboseRing)
        pass
    ivm.autoTorsion()
    return

def cartesianTopology(ivm,
                      sel="all",oTensors=[]):
    """configure the <m ivm>.IVM tolopogy for Cartesian dynamics/minimization -
    for the specified selection.

    This consists of breaking topological bonds, and specifying that all atoms
    are tree ``bases.''

    This function should be called after any custom changes are made to the
    ivm's topology setup, but before torsionTopology().

    oTensors is a list of <m varTensor>.VarTensor objects used for alignment
    tensors.
    """
                    
    for action in ivm.configActions:
        action.runAction()
        pass

    for topologyFunc,name in topologyFuncs:
        topologyFunc(ivm)
        pass
    
    import varTensorTools, prePotTools, planeDistTools, diffPotTools
    import relaxRatioPotTools, ensWeightsTools
    varTensorTools.topologySetup(ivm,oTensors)
    prePotTools.topologySetup(ivm)
    planeDistTools.topologySetup(ivm)
    diffPotTools.topologySetup(ivm)
    relaxRatioPotTools.topologySetup(ivm)
    ensWeightsTools.topologySetup(ivm)

    from atomSel import intersection
    updatePseudoAtoms()
    from selectTools import convertToAtomSel
    sel = convertToAtomSel(sel)
    ivm.breakAllBondsIn(intersection(sel,
                                     AtomSel("not pseudo",sel.simulation())))
    ivm.setBaseAtoms(sel)
    return


def genExtendedStructure(pdbFilename=0,
                         sel=0,
                         verbose=0,
                         maxFixupIters=500
                         ):
    """Generate an arbitrary extended conformation with ideal covalent geometry.

    This assigns X, Y, and Z coordinates to each atom, and then corrects the
    covalent geometry with <m protocol>.fixupCovalentGeom(), using maxFixupIters
    (an integer) as the maxIters argument to that function (maxFixupIters=0
    results in non-ideal covalent geometry).  Note that regardless of the
    choice of maxFixupIters, the nonbonded interactions may not be totally
    satisfied.

    The Y and Z coordinates are random (but small enough (within a range of -0.5
    to 0.5) to allow bonded atoms to form their bonds) and the X coordinate is
    the atom number divided by 10.  This will result in an extended
    configuration along the X axis.

    If a string is provided as pdbFilename, the generated coordinates will be
    written to a pdb file named pdbFilename.  (Note that if a pdb file with that
    name already exists, it will be used by <m protocol>.fixupCovalentGeom(), 
    and the resulting coordinates will update the input file.)
    """

    from atomSel import AtomSel
    updatePseudoAtoms()
    if not sel: sel = AtomSel("not pseudo")

    if type(sel)==type('string'): sel = AtomSel(sel)

    xSim = getXplorSimulation(sel.simulation())

    from pdbTool import PDBTool
    import os
    
    if pdbFilename and os.path.exists(pdbFilename):
        PDBTool(pdbFilename,sel).read()
    else:
        from atomSelAction import SetProperty
        import random
        from vec3 import Vec3
        for atom in sel:
            atom.setPos( Vec3(float(atom.index())/10,
                              random.uniform(-0.5,0.5),
                              random.uniform(-0.5,0.5)) )
            pass
        pass
    

    if verbose:
        print("fixing covalent geometry...")
        pass

    if maxFixupIters>0:
        fixupCovalentGeom(useVDW=1,maxIters=maxFixupIters,sel=sel,dynRatio=10,
                          verbose=verbose)

    if pdbFilename: PDBTool(pdbFilename,sel).write()

    return
        

    
def learnParams(selection="all",
                which="bond angle improper",
                bondForceConstant=1000,
                angleForceConstant=500,
                improperForceConstant=500,
                outFilename="new.par"
                ):
    """
    Read covalent parameters from the currently loaded structure. By default,
    equilibrium values for bond, angles, and improper dihedral angles are
    learned, but a subset of these can be specified using the which argument,
    which is a space-separated string of values. This learn prodecure does not
    determine force constants- they are specified explicitly using the
    appropriate argument. This function requires that coorindates be read,
    and that the associated covalent bond, angles, and impropers have been
    defined. 
    """

    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)

    from xplorSimulation import getXplorSimulation
    xsim=getXplorSimulation(selection.simulation())

    selString = selection.string()

    paramCmd=""
    which = which.lower()
    if 'bond' in which :
        paramCmd += "    BOND  (%s) (%s) %.2f TOKEN\n" % \
                    tuple([selString]*2+[bondForceConstant])
        pass
    if 'angle' in which :
        paramCmd += "    ANGLE (%s) (%s) (%s) %.2f TOKEN\n" % \
                    tuple([selString]*3+[angleForceConstant])
        pass
    if 'improper' in which :
        paramCmd += "    IMPR  (%s) (%s) (%s) (%s) %.2f TOKEN TOKEN\n" % \
                    tuple([selString]*4+[improperForceConstant])
        pass
    

    xsim.command("""

    flags exclude * include  bond angles impropers end
    parameters
      learn initiate sele=(%s) mode=nostatistics end
      learn accumulate end
      learn terminate end
    end

    parameters
    
      %s
    
      reduce MODE=AVER end
    ! ?
    end

    energy end
    """ % ( selString,paramCmd ))

    xsim.command('REMARKS  autogenerated by protocol.learnParams')


    xsim.command('write param output="%s" end' % outFilename)
    return
