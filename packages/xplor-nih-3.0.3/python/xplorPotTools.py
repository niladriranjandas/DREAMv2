
"""
Tools to help manipulate and analyze XPLOR energy terms
"""

def calcPotData(term):
    """return PotData values for the given XPLOR term
    """
    data = PotData(term)
    return (data.rmsd, data.violations, data.numRestraints)

            

from xplorSimulation import getXplorSimulation    

from simulationTools import mktemp
class PotData:
    def __init__(s,pot):
        s.name = pot.instanceName()
        s.rmsd = -1
        s.violations = -1
        s.numRestraints = -1
        s.outputString =''

        xSim = getXplorSimulation(pot.simulation())
        
        outputState=xSim.disableOutput()

        import os
    
        cmd=''

        #standard terms
        if s.name in ('BOND','ANGL','IMPR','CDIH','NOE','DIHE'):
            cmd = r'"print thres=%f %s" % (pot.threshold(),s.name)'
            pass
        #inverted form
        if s.name in ('COUP','HBDA','CARB'):
            cmd = r'"%s print thres=%f end end" % (s.name,pot.threshold())'
            pass

        if cmd:
            s.outputString += "\n   Violations in XPLOR term %8s\n" % s.name
            s.outputString += "             threshold: %f\n\n" % \
                              pot.threshold()
            xplorFilename=mktemp('xplor-violations')
            xSim.fastCommand("set print %s end" % xplorFilename)
            (rms,violations,number) = xSim.fastCommand(eval(cmd),
                                          ("result","violations","number"))
            try:
                s.rmsd = float(rms)
                s.violations = int(violations)
                s.numRestraints=int(number)
            except ValueError:
                pass

            xSim.fastCommand("set print $prev_print_file end" )
            xSim.fastCommand("close %s end" % xplorFilename)
            s.outputString += open(xplorFilename).read()
            os.unlink(xplorFilename)
            pass

        if s.name=='VEAN':
            cmd = r'"%s print thres=%f all end" % (s.name,pot.threshold())'

            s.outputString += "\n   Violations in XPLOR term %8s\n" % s.name
            s.outputString += "             threshold: %f\n\n" % \
                              pot.threshold()
            xplorFilename=mktemp('xplor-violations')
            xSim.fastCommand("set print %s end" % xplorFilename)
            (rms,violations) = xSim.fastCommand(eval(cmd),
                                          ("result","violations"))
            try:
                s.rmsd = float(rms)
                s.violations = int(violations)
            except ValueError:
                pass

            xSim.fastCommand("set print $prev_print_file end" )
            xSim.fastCommand("close %s end" % xplorFilename)
            output = open(xplorFilename).read()
            #the next line is rather fragile
            s.numRestraints = int(output.split()[-2])
            s.outputString += output
            os.unlink(xplorFilename)
            pass

        if s.name=='COLL':
            s.outputString += "\n   XPLOR term %8s analysis\n\n" % s.name

            xplorFilename=mktemp('xplor-violations')
            xSim.fastCommand("set print %s end" % xplorFilename)

            xSim.fastCommand("collapse print end")

            xSim.fastCommand("set print $prev_print_file end" )
            xSim.fastCommand("close %s end" % xplorFilename)
            s.outputString += open(xplorFilename).read()
            os.unlink(xplorFilename)
            pass
            
        if s.name=='VDW':
            #vdw
            from nonBondTools import vdwViolations
            viols = vdwViolations(pot.threshold())
            s.violations = len(viols)
            pass

        if s.name=='NCS':
            s.outputString += "\n   XPLOR term %8s \n\n" % s.name

            xplorFilename=mktemp('xplor-violations')
            xSim.fastCommand("set print %s end" % xplorFilename)

            xSim.fastCommand("ncs restraints ? end")

            xSim.fastCommand("set print $prev_print_file end" )
            xSim.fastCommand("close %s end" % xplorFilename)
            s.outputString += open(xplorFilename).read()
            os.unlink(xplorFilename)
            pass
            
        if s.name=='HBDB':
            pot.calcEnergy()
            s.outputString += "\n   XPLOR term %8s \n\n" % s.name

            xplorFilename=mktemp('xplor-violations')
            xSim.fastCommand("set print %s end" % xplorFilename)

            xSim.fastCommand("HBDB print end")

            xSim.fastCommand("set print $prev_print_file end" )
            xSim.fastCommand("close %s end" % xplorFilename)
            s.outputString += open(xplorFilename).read()
            os.unlink(xplorFilename)
            pass
            
        if s.name=='XREF':
            pot.calcEnergy()
            s.outputString += "\n   XPLOR term %8s \n\n" % s.name

            xplorFilename=mktemp('xplor-violations')
            xSim.fastCommand("set print %s end" % xplorFilename)

            R, TEST_R = xSim.fastCommand("""XREF
              print rfactor
              end
            end""",("R","TEST_R"))

            s.Rwork = float(R)
            s.Rfree = float(TEST_R)
            

            xSim.fastCommand("set print $prev_print_file end" )
            xSim.fastCommand("close %s end" % xplorFilename)
            s.outputString += open(xplorFilename).read()
            s.outputString += "\n"
            s.outputString += "R-work: %.3f   R-free: %.3f" % (s.Rwork,
                                                               s.Rfree)
            s.outputString += "\n"
            os.unlink(xplorFilename)
            pass
            


        xSim.enableOutput(outputState)
        return
    pass

        
def addRestraintsMethod(term):
    """
    add a restraints() method
    """
    from xplorPot import realXplorPot
    realXplorPot.restraints = restraints_meth
    return


def restraints_meth(pot):
    #  parse through restraints string
    # for CDIH
    #   each restraint is prefixed by the literal string
    #  ========================================
    # followed by 6 lines
    #  the first four will make up a name
    #  amount of violation is the final number on the 5th line.
    # support for showAll violations will require that a separate printout
    # be performed.

    #ANGL: each line (name) float float diff
    # amount of violation is abs(diff)

    ret=[]

    xSim = getXplorSimulation(pot.simulation())

    outputState=xSim.disableOutput()

    if pot.instanceName() in ("CDIH","ANGL","VEAN","BOND","IMPR","DIHE"):
        if pot.instanceName()=="VEAN":
            cmd=r'"%s print thres=0 all end"' % pot.instanceName()
        else:
            cmd=r'"print thres=-1 %s"' % pot.instanceName()
            pass
        xplorFilename=mktemp('xplor-violations')
        xSim.fastCommand("set print %s end" % xplorFilename)
        xSim.fastCommand(eval(cmd))
        xSim.fastCommand("set print $prev_print_file end" )
        xSim.fastCommand("close %s end" % xplorFilename)
    
        lines=[l[:-1] for l in open(xplorFilename).readlines()]
        import os
        os.unlink(xplorFilename)
        cnt=0
        if pot.instanceName()=="CDIH":
            while cnt < len(lines):
                if '==============' in lines[cnt]:
                    name = lines[cnt+1] + lines[cnt+2] + lines[cnt+3] + lines[cnt+4]
                    diff = float(lines[cnt+5].split()[-1])
                    cnt += 6
                    ret.append(PY_Restraint(pot,name,diff))
                    pass
                cnt += 1
                pass
            pass
        if pot.instanceName() in ("BOND","ANGL","IMPR","DIHE"):
            for line in lines:
                line = line.lstrip()
                if line.startswith('('):
                    pFields = line.split(")")
                    name = pFields[0][1:]
                    try:
                        diff = float(pFields[1].split()[2])
                        ret.append(PY_Restraint(pot,name,diff))
                    except ValueError:
                        pass
                    pass
                pass
            pass
        elif pot.instanceName()=="VEAN":
            for line in lines:
                if not line.startswith('         '):
                    continue
                line = line.strip()
                name=line[:79]
                nums=line[79:].split()
                diff=float(nums[5])
                ret.append(PY_Restraint(pot,name,diff))
                pass
            pass
        pass
    elif pot.instanceName()=="VDW":
        from nonBondTools import vdwViolations
        violations = vdwViolations(pot.threshold())
        
        for index in range(len(violations)):
            v=violations[index]
            ret.append(PY_Restraint(pot,v.name(),v.diff()))
            pass
        pass
    else:
        ret= "not supported"

    xSim.enableOutput(outputState)
    return ret

class PY_Restraint:
    def __init__(s,pot,name,diff):
        s.pot = pot
        s.name_ = name
        s.diff_ = diff
        return
    def name(s):
        return s.name_
    def diff(s):
        return s.diff_
    def violated(s):
        return abs(s.diff())>=s.pot.threshold()
    pass

def getHBDBbonds():
    """Return info from the XPLOR HBDB energy term on each backbone H-bond in a
    list of named tuples. The tuples have entries:
      resid_CO, resid_HN, dir_type, lin_type, dist_OH, ang_COH, ang_CaCO, ang_OHN
      Edir,  Elin"""
    
    from xplorPot import XplorPot
    data = PotData( XplorPot('HBDB') )

    from collections import namedtuple
    RetType = namedtuple('HBDBbond',
                         """resid_CO resid_HN dir_type lin_type dist_OH ang_COH
                         ang_CaCO ang_OHN Edir Elin""".split())
      
    import re

    #helper function from https://stackoverflow.com/questions/1265665/how-can-i-check-if-a-string-represents-an-int-without-using-try-except
    g_intRegex = re.compile(r"[-+]?\d+(\.0*)?$")
    def isInt(v):
        return g_intRegex.match(str(v).strip()) is not None
    
    ret = []
    num = 0
    for line in data.outputString.split('\n'):
        if line.startswith('#'):
            num = int(line.split()[6])
            continue

        fields = line.split()
        if len(fields)<10:
            continue
        try:
            if isInt(fields[0]):
                resid_CO  = int(line[4:8])
                segid_CO  = line[9:12].strip()
                resid_HN  = int(line[12:16])
                segid_HN  = line[18:21].strip()
                dir_type  = line[21:33]
                lin_type  = line[34:39]
                dist_OH   = float(fields[-6])
                ang_COH   = float(fields[-5])
                ang_CaCO  = float(fields[-4])
                ang_OHN   = float(fields[-3])
                Edir      = float(fields[-2])
                Elin      = float(fields[-1])
                ret.append( RetType(resid_CO, resid_HN, dir_type, lin_type,
                                    dist_OH, ang_COH, ang_CaCO, ang_OHN,
                                    Edir, Elin) )
                pass
            pass
        except ValueError:
            print("Error processing HBDB line:",line)
            pass
            
            
        pass
    if len(ret) != num:
        print("Warning. Found mismatch between number of H-bonds (%d) "%num + \
              "and length of list (%d)" %len(ret))
        pass
        
    return ret


def analyze(potList):
    """perform analysis of XplorPot terms and return nicely formatted summary.

    Note on VDW (nonbonded) violations:
      The analysis here does not heed any XPLOR constraints interaction
      statements: all nonbonded violations are always reported.
    """

    ret = ""

    from simulationTools import getPotTerms
    potList = getPotTerms(potList,'XplorPot')

    if not potList: return ret

    potList.sort(key=lambda x: x.instanceName())

#    ret += "%-9s  %6s  %6s  %6s\n" % \
#           ( "", "RMS", "Devia", "Viols")


    for pot in potList:
        name = pot.instanceName()

        if name=='VDW':
            from nonBondTools import vdwViolations
            violations = vdwViolations(pot.threshold())

            print()
            print("  Nonbonded Violations")
            print("        threshold: %f" % pot.threshold())
            print("        XplorPot.scale: %.2f" % pot.scale())
            print()
            print("%26s %26s %6s  %7s" % ("    ", "     ",
                                         "actual", "allowed"))
            print("%26s %26s %6s  %6s  %6s" % ("atom1     ", "atom2     ",
                                               " dist", " dist ", "diff  "))
            print("_"*78)

            for v in violations:
                print("%26s %26s %6.2f %6.2f %6.2f" % (v.atomi.string(),
                                                       v.atomj.string(),
                                                       v.dist, v.dist0,
                                                       v.dist0-v.dist))
                pass
#            ret += "%-9s  %6.3f  %6.3f  %6d\n" % \
#                   (name , 0., 0., len(violations) )
            pass
        else:   #non VDW terms
            
            data = PotData(pot)


            print(data.outputString)
            print("XplorPot.scale: %.2f" % pot.scale())
            pass
            
#            ret += "%-9s  %6.3f  %6.3f  %6d\n" % \
#                   (name , data.rmsd, 0., data.violations )
        if name=='XREF':
            ret += "  XREFINE R-work: %.3f  R-free: %.3f\n" % (data.Rwork,
                                                               data.Rfree)
            pass
        pass
    
    return ret


from simulationTools import registerTerm
registerTerm(analyze,"Xplor terms","Xplor",
r"""
These terms generate no extra header information.
""")

#default XPLOR threshold values
#  --used when XplorPots are constructed
defaultThresholds = {}
defaultThresholds['BOND'] = 0.05
defaultThresholds['ANGL'] = 2
defaultThresholds['IMPR'] = 2
defaultThresholds['DIHE'] = 10
defaultThresholds['CARB'] = 1.
defaultThresholds['CDIH'] = 5
defaultThresholds['NOE']  = 0.5
defaultThresholds['COUP'] = 1.
defaultThresholds['HBDA'] = 0.
defaultThresholds['VDW'] = 0.2
defaultThresholds['VEAN'] = 5
