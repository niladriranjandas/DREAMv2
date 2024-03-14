

"""
Tools for explicit water refinement with additional NMR restraints along the
lines of the RECOORD work:

Aart J. Nederveen, Jurgen F. Doreleijers, Wim Vranken, Zachary Miller,
Chris A.E.M. Spronk, Sander B. Nabuurs, Peter Guentert, Miron Livny,
John L. Markley, Michael Nilges, Eldon L. Ulrich, Robert Kaptein and
Alexandre M.J.J. Bonvin, ``RECOORD: a REcalculated COORdinates
Database of 500+ proteins from the PDB using restraints from the
BioMagResBank,'' Proteins 59, 662-672 (2005)..
"""

import xplor
from potList import PotList

def refine(heatingParams=[],
           highTempParams=[],
           coolingParams=[],
           potList=[],
           outFilename="waterRefine.pdb",
           rigidRegions=(),
           fixedRegions=(),
           keepWaters=False,
           waterResname='TIP3'
           ):
    """

    
    An implementation of the water refinement protocol from Aria 1.2. M.
    Nilges and J. Linge. This version can accommodate additional NMR
    restraints, such as dipolar coupling, etc.
    

    heatingParams 
    highTempParams
    coolingParams  - lists of ramped parameters to use during the heating,
                     high-temperature and cooling phases of the refinement
                     protocol.

    potList      - potential terms to use in water refinement. To these are
                   added the following terms: ELEC, VDW, BOND, ANGL, IMPR,
                   DIHE. 

    rigidRegions - selections of atoms which do not move relative to
                   each other.
    
    fixedRegions - selections of atoms which do not move at all.

    outFilename  - name of pdb file written after water refinement.

    keepWaters   - if True, do not delete water molecules when writing out
                   pdb files.

    waterResname - residue name of water molecules [default: 'TIP3']
    
    """

    finomega = 100

    from xplorPot import XplorPot

    if not potList: potList =PotList()
    for term in ('ELEC',
                 'VDW' ,
                 'BOND',
                 'ANGL',
                 'IMPR',
                 'DIHE',):
        if not term in list(potList.keys()):
            potList.append(XplorPot(term))
            pass
        pass

    xplor.command("""
    !FIX: nonbonded setup with full lj potential
    parameter 
      nbonds
        nbxmod=5 atom cdiel shift 
        cutnb=9.5 ctofnb=8.5 ctonnb=6.5 eps=1.0 e14fac=0.4 inhibit 0.25
        wmin=0.5
        tolerance  0.5
      end
    end""")


    
    #FIX: load parameters??
    import protocol

# fix: is this for ambiguous disulfide bonds?
#    xplor.command("""
#
#    
#    if ($toppar.ss_ambigunambig eq	"ambigous") then
#    parameter  nbfix S  S  462  13.6  462  13.6 end
#    end if
#    
#
#    """)

    xplor.command("delete sele=(resname %s) end" % waterResname)
    protocol.updatePseudoAtoms()
    
    buildShell(waterResname=waterResname)

    protocol.initDihedrals(reload=True)

    #coordinates for harmonic restraint potential
    xplor.command("""
    vector do (refx = x) (all)
    vector do (refy = y) (all)
    vector do (refz = z) (all)
    """)


    from atomSel import AtomSel
    from ivm import IVM

    ivm = IVM()


    waterSel="resname %s" % waterResname
    xplor.command("""
    ! reduce improper and angle force constant for some atoms
    evaluate ($kangle = 200)
    evaluate ($kimpro = 5)
    evaluate ($komega = 5)
    parameter
!    angle    (%s) (%s) (%s) 200 TOKEN
    angle    (not %s)(not %s)(not %s) $kangle  TOKEN
    improper (all)(all)(all)(all) $kimpro  TOKEN TOKEN
!    dihedral (all)(all)(all)(all) 3        TOKEN TOKEN
    end
    """ % tuple([waterSel]*6))

    # fix the protein for initial minimization
    ivm.fix("not %s" % waterSel)
    protocol.cartesianTopology(ivm)
#    #freeze out internal water motion
#    for oAtom in AtomSel("name O* and resname %s" % waterResname):
#        ivm.group("byresidue (resid %d)" % oAtom.residueNum())
#        pass
    protocol.initMinimize(ivm, numSteps=100, dEPred=100,
                          potList=potList)
    ivm.run()

    ivm.reset()
    #freeze out internal water motion
    for oAtom in AtomSel("name O* and %s" % waterSel):
        ivm.group('segid "%s" and resid %d' % (oAtom.segmentName(),
                                               oAtom.residueNum()  ))
        pass

    #configure for Cartesian dynamics/minimization
    for region in rigidRegions:
        ivm.group( region )
        pass
    for region in fixedRegions:
        ivm.fix( region )
        pass
    
    protocol.cartesianTopology(ivm)

    potList.append(XplorPot("HARM"))
    ivm.setPotList( potList )
    xplor.command("""
    ! release protein and restrain harmonically
    vector do (refx=x) (all)
    vector do (refy=y) (all)
    vector do (refz=z) (all)
    restraints harmonic 
       exponent = 2
    end
    vector do (harm = 0)  (all)
    vector do (harm = 10) (not name h*)
    """)
    
    protocol.initMinimize(ivm,
                          numSteps=100,dEPred=10)
    for i in range(10):
        ivm.run()

        xplor.command("""
        vector do (refx=x) (not PSEUDO)
        vector do (refy=y) (not PSEUDO)
        vector do (refz=z) (not PSEUDO)""")
        pass
    
    
    from atomAction import randomizeVelocities, SetProperty
    protocol.massSetup(atomicMass=50,
                       friction=20)

    kharm=50
    # heat to 500 K
    highTemp=500
    for bath in (100, 200, 300, 400, highTemp):
        xplor.command("vector do (harm = %f) (not name h* and not PSEUDO)"
                      %kharm)

        protocol.initDynamics(ivm,numSteps=100,stepsize=0.002,
                              finalTime=0,
#                              eTol_factor=0.0001,
                              bathTemp=bath,
                              printInterval=50)
        ivm.run()
        kharm = max(0, kharm - 4)

        xplor.command("""
        vector do (refx=x) (not PSEUDO)
        vector do (refy=y) (not PSEUDO)
        vector do (refz=z) (not PSEUDO)
        """)
        pass

    
    #
    #
    # refinement at high T:

    potList['DIHE'].setScale(2)

    potList.remove('HARM')

    protocol.initDynamics(ivm,finalTime=25,
                          numSteps=2000,
                          potList=potList,
#                          eTol_factor=0.0001,
                          bathTemp=highTemp,initVelocities=1)
    ivm.run()
                          
    potList['DIHE'].setScale(3)

    from simulationTools import LinRamp, AnnealIVM
    #
    #
    # cool
    bath=500
    rampedParams=coolingParams[:]

    global kbonds; kbonds=1000
    global kangle; kangle= 200
    global kimpro; kimpro=   5
    global kchira; kchira=   5
    global komega; komega=   5
    def stepInit(annealIVM):

        global kbonds; kbonds    = max(500,kbonds   / 1.1)
        global kangle; kangle    = min(500,kangle   * 1.1)
        global kimpro; kimpro    = min(500,kimpro   * 1.4)
        global kchira; kchira    = min(800,kchira   * 1.4)
        global komega; komega    = min(200,komega   * 1.4)

        from atomAction import randomizeVelocities
        randomizeVelocities( annealIVM.ivm.bathTemp() )

#        rampedParams.append(LinRamp(1000 ,500,"global kbonds; kbonds=VALUE") )
#        rampedParams.append(LinRamp( 200 ,500,"global kangle; kangle=VALUE"))
#        rampedParams.append(LinRamp(   5 ,500,"global kimpro; kimpro=VALUE") )
#        rampedParams.append(LinRamp(   5 ,800,"global kchira; kchira=VALUE"))
#        rampedParams.append(LinRamp(   5 ,200,"global komega; komega=VALUE"))
        xplor.command("""
        eval ($kbonds=%f)
        eval ($kangle=%f)
        eval ($kimpro=%f)
        eval ($kchira=%f)
        eval ($komega=%f)
        """ % (kbonds,kangle,kimpro,kchira,komega))

        xplor.command("""
        parameter
  bond     (not %s and not name H*) (not %s and not name H*) $kbonds  TOKEN
  angle    (not %s and not name H*) (not %s and not name H*) (not %s and not name H*) $kangle  TOKEN

  improper (all)(all)(all)(all) $kimpro  TOKEN TOKEN
  improper (name HB and resn VAL)(name CA and resn VAL) (name CG1 and resn VAL)(name CG2 and resn VAL) $kchira TOKEN TOKEN
  improper (name HB and resn THR)(name CA and resn THR) (name OG1 and resn THR)(name CG2 and resn THR) $kchira TOKEN TOKEN
  improper (name HG and resn LEU)(name CB and resn LEU) (name CD1 and resn LEU)(name CD2 and resn LEU) $kchira TOKEN TOKEN
  improper (name HB and resn ILE)(name CA and resn ILE) (name CG2 and resn ILE)(name CG1 and resn ILE) $kchira TOKEN TOKEN
  improper (name HA)(name N)(name C)(name CB) $kchira TOKEN TOKEN

  improper (name O)  (name C) (name N) (name CA) $komega TOKEN TOKEN
  improper (name HN) (name N) (name C) (name CA) $komega TOKEN TOKEN
  improper (name CA) (name C) (name N) (name CA) $komega TOKEN TOKEN
  improper (name CD) (name N) (name C) (name CA) $komega TOKEN TOKEN
        end""" % tuple([waterSel]*5))
        return
    

    protocol.initDynamics(ivm,
#                          eTol_factor=0.0001,
                          numSteps=500,finalTime=10)
    AnnealIVM(initTemp=highTemp,finalTemp=25,
              tempStep=25,
              ivm=ivm,
              rampedParams=rampedParams,
              extraCommands=stepInit).run()
                   
    
        #  aria methylene NOE assignment swapping..
        #   evaluate ($swap = 1.0) 
        #   inline @RUN:protocols/swap.cns
        
        #
    #final minimization:
    protocol.initMinimize(ivm,numSteps=200)
    ivm.run()
    #

    xplor.command("""
       constraints interaction
         (not (%s) and not PSEUDO)
         (not (%s) and not PSEUDO)
       end
       """ % tuple([waterSel]*2))

    sel="all"
    if not keepWaters:
        sel="not (resname %s)" % waterResname
        pass
        
    from pdbTool import PDBTool
    pdbFile = PDBTool(outFilename,sel)

    from simulationTools import analyze
    pdbFile.addRemarks(analyze(potList,
                               outFilename=outFilename+'.viols'))

    pdbFile.write()

    protocol.updatePseudoAtoms()
    protocol.initDihedrals(reload=True)
    #energy end
    #
    #evaluate ($filename= "NEWIT:water/" + $file - "PREVIT:" - ".pdb" + "w.pdb")
    #
    #@RUN:protocols/print_coorheader.cns
    #
    #do (q=1) (all)
    #write coordinates sele= (not resn WAT) output =$filename end
    #      
    #evaluate ($filename= "NEWIT:water/" + $file - "PREVIT:" - ".pdb" + "w.float")
    #
    #set print $filename end
    #aria float store end end
    #
    #""")


def segidCountMap(cnt):
    """Given an integer return a character or throw an exception if out of
    range. The map looks like:
    0..9 --> '0'..'9'
    10..35 --> 'A'..'Z'
    """

    if cnt<10:
        return str(cnt)
    if cnt>35:
        raise Exception("Segid count too large")
        pass
    return chr(cnt-10+ord('A'))
    
    

def buildShell(thickness=8.0,
               pw_dist = 4.0,
               dyncount = 1,
               waterBoxPDB = "TOPPAR:waterRef/boxtyp20.pdb",
               waterBoxLength = 18.856,       # length of Brooks' water box
               waterDiam = 2.4,               # diameter of water molecule
               waterResname='WAT'
               ):
    """
     water soaking protocol from Aria 1.2 (M. Nilges and J. Linge)

     options:
      thickness -  max initial water-protein(heavy atom) distance 
      pw_dist   -  min initial water-protein(heavy atom) distance 
      dyncount  -  iteration number - used to generate a unique segment name

    """

    import xplor
    from atomSel import AtomSel
    from pdbTool import PDBTool
    
    (xplorEcho, xplorMess) = xplor.command("set echo off mess off end",
                                             ('prev_echo','prev_messages'))
    import protocol
    protocol.initTopology("water")
    protocol.initParams("water")



    from atomAction import SetProperty
    from atomSelAction import Translate
    getBounds=GetBounds()
    AtomSel("(not resn %s) and not PSEUDO"%waterResname).apply(getBounds)
    (xmin,xmax,
     ymin,ymax,
     zmin,zmax) = getBounds.bounds()


    # determine how many boxes are necessary in each dimension
    xbox = int( (xmax-xmin + 2*(thickness + waterDiam)) / waterBoxLength + 0.5)
    ybox = int( (ymax-ymin + 2*(thickness + waterDiam)) / waterBoxLength + 0.5)
    zbox = int( (zmax-zmin + 2*(thickness + waterDiam)) / waterBoxLength + 0.5)

    xmtran =  xmax + thickness - waterBoxLength/2 + waterDiam
    ymtran =  ymax + thickness - waterBoxLength/2 + waterDiam
    zmtran =  zmax + thickness - waterBoxLength/2 + waterDiam


    #--------------------------------------------------
    # read in the same box of water several times, and move it around
    # so as to cover all the space around the site of interest.
    # take into account box offset

    xcount=0
    xtrans = xmin - thickness - waterDiam - waterBoxLength/2 
    while xtrans < xmtran:
      xcount += 1
      xtrans += waterBoxLength

      ycount=0
      ytrans = ymin - thickness - waterDiam - waterBoxLength/2 
      while ytrans < ymtran:
        ycount += 1
        ytrans += waterBoxLength

        zcount=0                 
        ztrans = zmin - thickness - waterDiam - waterBoxLength/2
        while ztrans < zmtran:
          zcount += 1
          ztrans += waterBoxLength

          xplor.command('''
          segment
            name="W000"
            chain
              coordinates @%s
            end
          end''' % waterBoxPDB)
          from os import environ
          
          waterBoxPDBfilename=waterBoxPDB.replace('TOPPAR:',
                                                  environ['TOPPAR']+'/')
          PDBTool(waterBoxPDBfilename).read()
          AtomSel("segid W000").apply(Translate((xtrans, ytrans, ztrans)))

          # all new water oxygens
          waterSel = "segid W000 and name oh2"
          # all new water oxygens close to a protein heavy atom
          waterNeighboringProtein = """(%s) and
              (not (resn %s or PSEUDO or name H*)) around %f""" % \
              (waterSel,waterResname,pw_dist)
          # all new water oxygens close to old water oxygens
          waterNeighboringWater = """(%s) and 
              (segid wat# and not name H*) around %f""" % (waterSel,waterDiam)
          # water oxygens further than thickness away from a protein heavy atom
          farWater = """(%s) and
              not (not (resn %s or PSEUDO or name H*)) around %f""" % \
              (waterSel,waterResname,thickness)
          xplor.command("delete sele= (byres ((%s) or (%s) or (%s))) end" %
                     (waterNeighboringProtein,waterNeighboringWater,farWater))


          # give waters unique segid name
          segid = "W%s%s%s" % (segidCountMap(xcount),
                               segidCountMap(ycount),
                               segidCountMap(zcount))
          AtomSel("segid W000").apply(SetProperty('segmentName',segid))

          pass
      pass

    # give waters a unique resid so that we get the segid to play around with
    AtomSel("segid w* and not segid wat#").apply(RenumberResids())
    AtomSel("segid w* and not segid wat#").apply(SetProperty("segmentName",
                                                             "WAT%d" %
                                                                dyncount))

    # shave off any waters that are left
    xplor.command("""delete sele= (byres (name oh2 and
                                    not (not (resn %s or PSEUDO or name H*))
                                    around %f)) end""" % (waterResname,
                                                          thickness))
    xplor.command("set echo %s mess %s end" %(xplorEcho,xplorMess))

    return


from atomSelAction import PyAtomSelAction
from atom import Atom
class RenumberResids(PyAtomSelAction):
    """
    """
    def __init__(s):
        PyAtomSelAction.__init__(s,s)
        return
    def init(s,sel):
        s.count=0
    def run(s,sim,index):
        Atom(sim,index).setResidueNum(s.count // 3 + 1)
        s.count+=1
        return
    pass


class GetBounds(PyAtomSelAction):
    """compute Cartesian bounds of atom selection.
    The method bounds() returns a tuple of
       (xmin,xmax,ymin,ymax,zmin,zmax)
    """
    def __init__(s):
        PyAtomSelAction.__init__(s,s)
        return
    def init(s,sel):
        s.xmin=1e30
        s.xmax=-1e30
        s.ymin=1e30
        s.ymax=-1e30
        s.zmin=1e30
        s.zmax=-1e30
    def run(s,sim,index):
        (x,y,z) = Atom(sim,index).pos()
        
        if x < s.xmin: s.xmin = x
        if x > s.xmax: s.xmax = x
        if y < s.ymin: s.ymin = y
        if y > s.ymax: s.ymax = y
        if z < s.zmin: s.zmin = z
        if z > s.zmax: s.zmax = z

        return
    def bounds(s):
        return ( s.xmin, s.xmax,
                 s.ymin, s.ymax,
                 s.zmin, s.zmax)
    pass

