"""Generate statistics about most-violated restraints over a series of
structures.
"""
from functools import reduce

class RestraintStats:
    def __init__(s):
        """ reset all counters. 
        """
        s.potTypes = {}
        s.structNums=[]
        return

    def accumulate(s,potList,structNum):
        """
        for each restraint in each potential term, determine whether it is
        violated and save per-structure info.
        
        call this routine once for each structure
        
        this processes energy terms with a restraints() method, and also
        handled nonbonded interactions as a special case.
        """

        from simulationTools import potType
        potLists = [p for p in potList if potType(p)=='PotList']

        from potList import PotList
        tot=PotList('TOTAL')
        for term in potList: tot.append(term)
        potList=list(potList)
        potList.append(tot)

        s.structNums.append(structNum) # should warn if structNum already seen

        from simulationTools import flattenPotList
        subTerms = flattenPotList(potLists)

        for (isSubTerm,terms) in ((0,potList),
                                 (1,subTerms)):
            for term in terms:
                potName = potType(term)
                instName = term.instanceName()
                if potName not in s.potTypes: s.potTypes[potName] = {}
                if instName not in s.potTypes[potName]:
                    s.potTypes[potName][instName] = TermEntry(term,isSubTerm)
                    pass

                getViolations(s.potTypes[potName][instName], term, structNum)
                pass
            pass
        
        return
    def collect(s,comm,
                structNums=None,
                deleteUnusedData=True):
        """
        Collect accumulated restraint stats from remote processes into
        process 0. While summarizeTerms is called from process 0, collect
        should be called directly from other processes. If deleteUnusedData
        is True, data corresponding to structures not in the StructNums list
        will be deleted, thus modifying the RestraintStats object.
        """
        
        for potName in list(s.potTypes.keys()):
            for instName in list(s.potTypes[potName].keys()):
                vData = s.potTypes[potName][instName]
                for xStats in list(vData.extraStats.values()):
                    del xStats.function
                    pass
                del vData.pot
                del vData.restraints
                pass
            pass

        structNums = comm.distribute(structNums)
        #delete data for all structures except those specified by structNums
        for potName in list(s.potTypes.keys()):
            for instName in list(s.potTypes[potName].keys()):
                vData = s.potTypes[potName][instName]
                for structNum in list(vData.energies.keys()):
                    if not structNum in structNums:
                        del vData.energies[structNum]
                        del vData.rmsds[structNum]
                        del vData.violations[structNum]
                        if structNum in vData.violList:
                            del vData.violList[structNum]
                            pass
                        for name in list(vData.extraStats.keys()):
                            del vData.extraStats[name].valDict[structNum]
                            pass
                        pass
                    pass
                pass
            pass


        #recollect structNums to obtain a proc->structNum map
        # this is slightly complicated because procs disappear if their
        # computation crashes in some unforeseen manner. So, we really
        # need a map of the currently existing procs and structures. As
        # structuNumsP and comm.procs() are both sorted, we get this in the
        # loop below in the variables index and proc.
        structNumsP = comm.collect(structNums)
        #convert from array of arrays to single array
        s.structNums = reduce(lambda x,y:x+y, structNumsP, [])

        if comm.procNum!=0:
            comm.writeDataTo(0,s.potTypes)
            return

        
        for index,proc in enumerate(comm.procs()):
            rPotTypes = comm.readDataFrom(proc)
            if rPotTypes==None:
                continue

            if len(structNumsP[index])==0:
                # in this case, there's garbage in rPotTypes
                continue

            #if proc 0 calculated no structures, then potTypes will be empty,
            # so fill it with remote data.
            if len(s.potTypes)==0:
                s.potTypes=rPotTypes
                continue
            
            for potName in list(s.potTypes.keys()):
                for instName in list(s.potTypes[potName].keys()):
                    vData = s.potTypes[potName][instName]
                    if ( potName not in rPotTypes or
                         instName not in rPotTypes[potName] ):
                        # this happens if no structures were successfully
                        # calculated by the remote process
                        continue
                    rData = rPotTypes[potName][instName]
                    vData.count += rData.count
                    for (structNum,energy) in list(rData.energies.items()):
                        vData.energies[structNum] = energy
                        pass
                    for (structNum,rmsd) in list(rData.rmsds.items()):
                        vData.rmsds[structNum] = rmsd
                        pass
                    for (structNum,
                         violations) in list(rData.violations.items()):
                        vData.violations[structNum] = violations
                        pass

                    for (structNum,
                         vRestraints) in list(rData.violList.items()):
                        vData.violList[structNum] = vRestraints
                        pass

                    for name in list(vData.extraStats.keys()):
                        valDict = rData.extraStats[name].valDict
                        for (structNum,
                             val) in list(valDict.items()):
                            vData.extraStats[name].valDict[structNum] = val
                            pass
                        pass
                    pass
                pass
            pass
        return
    
                
    def summarizeTerms(s,comm,
                       structNums=None):
        """
        Return a string with statistics for the accumulated potential term for
        structures specified by structNums (all accumulated structures if
        structNum=None).

        When summarizeTerms is called by process 0, collect must be called by
        all other processes.
        
        """
        if structNums==None: structNums=s.structNums

        s.collect(comm,structNums)

        from math import sqrt
        ret = ' Average values for potential terms:\n\n'
        potKeys = [k for k in list(s.potTypes.keys()) if k!='PotList']
        potKeys.sort()
        potKeys = [k for k in list(s.potTypes.keys()) if k=='PotList'] + potKeys
        ret += "%-13s %-9s  %15s  %13s    %12s     Num\n" % \
               ("type", "name","Energy(dev)","RMSD(dev)","viols(dev)")

        def printTermInfo(term):
            ret=''
            if term.noPotSummary:
                return ret

            (energy, sigmaE) = term.energyAveDev(structNums)
            ret += "%-13s %-9s %8.2f(%6.2f)" % \
                   (potKey,instKey,energy, sigmaE)#,term.minE,term.maxE)

            term.structNums = structNums

            (rmsd,sigmaRMSD) = term.rmsdAveDev(structNums)
            if rmsd>=0.:
                ret += "  %7.3f(%6.3f) " % (rmsd, sigmaRMSD)
            else:
                ret += "  %7s %6s  " % ('','')
                pass

            (viols,sigmaViols) = term.violationsAveDev(structNums)
            if viols>=0:
                ret += "  %6.1f(%5.1f)" % (viols, sigmaViols)
            else:
                ret += '  %13s'%' '
                pass

            numRestraints=term.numRestraints
            if numRestraints>0:
                ret += " %5d"%numRestraints
                pass
            ret += '\n'
            return ret
            
        output0=''
        output1=''
        for potKey in potKeys:
            terms = s.potTypes[potKey]
            instKeys = list(terms.keys())
            instKeys.sort()
            if potKey=="PotList":
                instKeys.remove("TOTAL") # TOTAL always comes first
                instKeys = ["TOTAL"] + instKeys
                pass
            for instKey in instKeys:
                term = terms[instKey]
                if term.isSubTerm:
                    output1 += printTermInfo(term)
                else:
                    output0 += printTermInfo(term)
                    pass
                pass
            pass
        ret += output0
        if output1:
            ret += "\nbreakdown of terms in PotLists:\n\n"
            pass
        ret += output1
        return ret
    

    def summarizeExtraQuantities(s,structNums=None):
        """
        Return a string with statistics for extra quantities for supported
        potential terms for the specified structures (all accumulated
        structures if structNums is None).
        """

        if structNums==None: structNums=s.structNums

        from math import sqrt
        ret = ''

        potTypes = list(s.potTypes.keys())
        potTypes.sort()

        def aveDevPotList(potList,terms,statName):

            sum=0
            sum2=0
            
            weightSum=0
            for pot in potList:
                termEntry=terms[pot.instanceName()]
                count = termEntry.count
                if not count: return (0,0)

                weight = pot.scale()
                sum += weight * termEntry.extraStats[statName].sum
                sum2 += weight * termEntry.extraStats[statName].sum2
                weightSum += weight
                pass

            ave=sum
            if weightSum>0:
                ave /= count * weightSum
                pass

            sigma=0
            if count>1 and weightSum>0:
                var = (sum2 -sum*ave) / (count-1) #FIX: check
                var /= weightSum
                if var>0.:
                    sigma = sqrt(var)
                    pass
                pass
            return (ave,sigma)
        
        from termAnalysis import extraStats
        for potType in potTypes:
        
            if potType not in extraStats:
                continue
            eStats = extraStats[potType]

            terms = s.potTypes[potType]

            ret += "\n%-13s " % potType
            names=[tup[0] for tup in eStats]
            for name in names:
                ret += "%s " % name.center(15)
                pass
            ret += '\n'

            #first treat non-PotList terms
            instKeys = list(terms.keys())
            instKeys.sort()
            for instKey in instKeys:
                term = terms[instKey]
                ret += "%-13s" % instKey
                for name in names:
                    if not term.isList or term.extraStats[name].supportsList:
                        aveDev = term.extraStatsAveDev(name,structNums)
                        if hasattr(aveDev[0],'__getitem__'):
                            ret += "{ "
                            ret += " , ".join(["%7.3f(%6.3f) " % (ave,dev) for
                                               ave,dev in aveDev])
                            ret += " }"
                        else:
                            ret += "%7.3f(%6.3f) " % aveDev
                            pass
                        pass
                    else:
                        ret += "%16s" % " "
                    pass
                ret += '\n'
                pass

            #second treat PotLists of this type of term
            if 'PotList' in s.potTypes:
                instKeys = list(s.potTypes['PotList'].keys())
                instKeys.sort()
                for potListName in instKeys:
                    if not s.potTypes['PotList'][potListName].potType==potType:
                        continue

                    term=s.potTypes['PotList'][potListName]
                    contents=""
                    for name in names:
                        if (name in term.extraStats and 
                            term.extraStats[name].supportsList):
                            (ave,dev) = term.extraStatsAveDev(name)
                            contents += "%7.3f(%6.3f) " % (ave,dev)
                        else:
                            contents += "%16s" % " "
                            pass
                        pass
                    if contents.strip():
                        ret += "%-13s%s\n" % (potListName,contents)
                        pass
                    pass
                pass
            pass
                        
            
        return ret
    

    def summarizeViolations(s,structNums=None,cutoff=0):
        """
        Summarize violation statistics for all potential terms for
        structures specified by structNums (all accumulated structures if
        structNum=None).
        
        A string is returned reporting violated restraints, sorted from most
        violated to least. The cutoff argument can be used to specify the
        fraction of violated structures below which to stop printing violated
        restraints. 
        
        In the returned string is the percentage of structures violated,
        the average restraint rmsd (for violated structures), the index of the
        restraint, and an identifying name.
        """

        if structNums==None: structNums=s.structNums
        
        from math import sqrt
        potKeys = [k for k in list(s.potTypes.keys()) if k!='PotList']
        potKeys.sort()
        potKeys = [k for k in list(s.potTypes.keys()) if k=='PotList'] + potKeys
    
        ret = '\n Summary of most violated restraints\n'
        for potKey in potKeys:
            terms = s.potTypes[potKey]
            instKeys = list(terms.keys())
            instKeys.sort()
            for instKey in instKeys:
                if terms[instKey].noViolationStats:
                    continue
                if terms[instKey].notSupported:
                    ret += "\n unsupported term: %s\n" % instKey
                    continue

                violList = getStructureTerms(terms[instKey].violList,
                                             structNums)
                #convert from array of arrays to single array
                violList = reduce(lambda x,y:x+y, violList, [])
   
                #indexed by name instead of index because indices are
                # not necessarily unique (e.g. vdw indices)
                names = [v.name for v in violList]
                rCount={}
                for n in names: rCount[n]=0
                rDiff2={}
                for n in names: rDiff2[n]=0.
                rIndices={}
                maxVal=0
                minVal=0
                for restraint in violList:
                    rCount[restraint.name] += 1
                    rDiff2[restraint.name] += restraint.diff**2
                    maxVal=max(maxVal,restraint.diff)
                    minVal=min(minVal,restraint.diff)
                    rIndices[restraint.name] = restraint.index
                    pass
                    
                nameCountPairs = list(rCount.items())
                # sort by count first
                nameCountPairs.sort(key=lambda x: x[1],reverse=True)
                #FIX: then should be sorted by index
            
                ret += "\nviolations in %s term %s" %(potKey,instKey)
                ret += "  ( %d structures )\n\n" % len(structNums)
                if maxVal>0 or minVal<0:
                    tmp=""
                    tmp += "  Maximum positive "
                    if minVal<0: tmp += "/ negative "
                    tmp += "violation amount: %.3f" % maxVal
                    if minVal<0: tmp += " / %.3f" % minVal
                    tmp += "\n\n"
                    ret += tmp
                    pass
                ret += "viol%  amount  index    restraint name\n"
                ret += "-"*75 + '\n'
    
                for (name,count) in nameCountPairs:
                    percent= float(count)/len(structNums)
                    if percent<=cutoff: break
                    percent *= 100
                    rmsd=0
                    if count:
                        rmsd= sqrt(rDiff2[name] / count)
                        pass
    
                    ret += "%6.1f  %5.2f  %4d  %s\n" % (percent,
                                                        rmsd,
                                                        rIndices[name],
                                                        name)
                    pass
                pass
            pass
    
        return ret

                
        
class RestraintEntry:
    def __init__(s,index,restraint):
        s.index=index
        s.name =restraint.name()
        s.diff=restraint.diff()
        return
    pass

class ExtraStats:
    def __init__(s,fun,supportsList):
        """fun: function to call on each term
        fun: whether this function supports a list of terms.
        """
        s.function = fun
        s.valDict={}
        s.supportsList=supportsList
        return
    pass

class TermEntry:
    def __init__(s,term,isSubTerm):
        """term is a Pot.
        isSubTerm denotes whether the term is a member of the toplevel PotList,
        or a member of a lower-level PotList.
        """
        s.noPotSummary=False
        s.noViolationStats=False
        s.restraints = {}
        s.count = 0
        s.notSupported=0
        s.energies={}
        s.rmsds={}
        s.violations={}
        s.violList={}

        s.pot=term
        s.isSubTerm=isSubTerm
        s.isList=False
        from simulationTools import potType
        from termAnalysis import extraStats
        s.potType=potType(term)
        s.extraStats={}

        #treat potLists
        if potType(term)=='PotList':
            samePot=True
            s.isList=True
            from simulationTools import flattenPotList
            terms = flattenPotList(term)
            if not len(terms):
                return
            s.potType=potType(terms[0])
            for pot in terms:
                if potType(pot)!=s.potType:
                    samePot=False
                    break
                pass
            if not samePot:
                return
            pass

        if s.potType not in extraStats:
            return

        for (name,fun,supportsList) in extraStats[s.potType]:
            if not s.isList or supportsList:
                s.extraStats[name] = ExtraStats(fun,supportsList)
                pass
            pass
        
        return
    def energyAveDev(s,structNums=None):
        """return tuple of mean and std dev of energy of this term for
        structures in structNums.
        """
        if structNums==None: structNums = s.structNums
        return aveDev( getStructureTerms(s.energies,structNums) )
    def extraStatsAveDev(s,name,structNums=None):
        """return tuple of mean and std dev of the corresponding quatity
        for structures in structNums.
        """
        if structNums==None: structNums = s.structNums
        return aveDev( getStructureTerms(s.extraStats[name].valDict,
                                         structNums) )
    def rmsdAveDev(s,structNums=None):
        """return tuple of mean and std dev rmsd of this term for 
        structures in structNums.
        """
        if structNums==None: structNums = s.structNums
        return aveDev( getStructureTerms(s.rmsds,structNums) )
    def violationsAveDev(s,structNums=None):
        """return tuple of mean and std dev number of violations for
        structures in structNums.
        """
        if structNums==None: structNums = s.structNums
        return aveDev( getStructureTerms(s.violations,structNums) )
        
    pass

def getStructureTerms(dict,keys):
    """
    Get all elements of the input dictionary corresponding to the
    given keys and return as a list.
    """
    ret=[]
    for key in keys:
        ret.append(dict[key])
        pass
    return ret


def aveDev(vector):
    """
    Return tuple of ave. and standard dev. of the input list of values, or
    array of values
    """
    if len(vector)==0: return (0,0)

    if hasattr(vector[0],'__getitem__'): # to deal with arrays
        transpose = [list(i) for i in zip(*vector)]
        #python2: list(map(list, map(None, *vector)))
        return [aveDev(v) for v in transpose]

    from math import sqrt
    from cdsVector import CDSVector_double as Vector
    from cdsVector import sum
    vector = Vector(vector)
    ave = sum(vector) / len(vector)
    sum2 = sum(vector**2)
    sigma=0
    if len(vector)>1:
        var = (sum2 - len(vector)*ave**2) / (len(vector)-1)
        if var>0.:
            sigma = sqrt(var)
            pass
        pass
    return (ave,sigma)
    

def getViolations(vdata,term,structNum):
    """
    obtain violation info. From term.restraint() method, or as a special case.
    Place all info in dictionary associated with structure structNum.
    """

    vdata.count += 1

    if hasattr(term,'noPotSummary'):
        vdata.noPotSummary = term.noPotSummary
        pass

    if hasattr(term,'noViolationStats'):
        vdata.noViolationStats = term.noViolationStats
        pass

    vdata.energies[structNum] = term.calcEnergy()

    vdata.numRestraints=term.numRestraints()

    if hasattr(term,"rms"):
        vdata.rmsds[structNum] = term.rms()

    if hasattr(term,"violations"):
        vdata.violations[structNum]  = term.violations()

    if hasattr(term,"restraints") and term.restraints()!="not supported":
        restraints = term.restraints()
        
        vdata.violList[structNum] = []
        for index in range(len(restraints)):
            restraint = restraints[index]
            if hasattr(restraint,'violated') and restraint.violated():
                vdata.violList[structNum].append(
                    RestraintEntry(index,restraint) )
                pass
            pass
        pass
    else:
        vdata.notSupported=1
        pass

    names=list(vdata.extraStats.keys())
    for name in names:
        val = vdata.extraStats[name].function(term)
        vdata.extraStats[name].valDict[structNum] = val
        pass

    return
