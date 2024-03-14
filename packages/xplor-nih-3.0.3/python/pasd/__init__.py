
"""
Stubs for future PASD NOE assignment support.
"""

#global default values

# peaks are binned by intensity into strong, medium, weak and very weak
#  the associated distance bounds for each bin are given below
#
distanceBins=[ (1.8, 2.7),
               (1.8, 3.3),
               (1.8, 5.0),
               (1.8, 6.0) ]

# contacts between residues which differ by more than this number will
# be considered long-range.
longRangeResidCutoff=6

# nMono is set to the number of symmetric multimers
#
nMono=1

# aveExp is the exponent used for distance averaging
#
aveExp=6


def readSTAR(filename,
             addSaveSet=False):
    """Read in the given NMR-STAR file and store the contents in the global
    starData variable.

    If addSaveSet=True, wrap the contents of the file with the lines
    save_dummy
      .
      .
      .
    save_
    """
    import cif

    contents = open(filename).read()
    if addSaveSet: contents = "save_dummy\n" + contents + "\nsave_\n"

    global starData
    starData=cif.Cif()
    starData.parse(contents)
    return

def starShifts(useAmbiguityCodes=False,
               saveSet=None,
               verbose=False,
               segmentName=None,
               residOffset=0,
               tclOutput=True):
    """
    Read chemical shift data from the global starData <m cif>.Cif
    object, which must have been previously populated by calling readStar().

    If useAmbiguityCodes is True, some effort is made to obey the
    included ambiguity codes.

    Use a non-zero residOffset if the residue numbering in the shift table
    differs from that in your PSF.
    """

    if not saveSet:
        saveSets=[]
        for key in list(starData.keys()):
            if "Atom_chem_shift" in list(starData[key].keys()):
                saveSets.append(key)
                pass
            pass
        if len(saveSets)==0:
            raise Exception("starShifts: "
                            "could not locate any chemical shift saveSets.")
            

        if len(saveSets)!=1:
            raise Exception("starShifts: too many chem shift saveSets. "
                            "Please choose one of: " +
                            str(saveSets))
        
        saveSet=saveSets[0]
        print("starShifts: using saveSet:", saveSet)
    else:
        if not saveSet in list(starData.keys()):
            raise Exception("starShifts: saveSet >%s< not found" % saveSet)
        pass
        

    shiftCat = starData[saveSet].Atom_chem_shift

    #
    # The NMR-STAR ambiguity code is
    #
    # 1:  unique
    # 2:  ambiguity of geminal atoms or geminal methyl proton groups
    # 3:  ambiguity of aromatic atoms on opposite sides of the ring (eg., phe HD1 | HD2)
    # 4:  intraresidue ambiguities (eg., lys HG vs HD protons
    # 5:  interresidue ambiguities (eg., lys 12 vs lys 27)
    # 9:  other ambiguity
    #
    # see http://www.bmrb.wisc.edu/elec_dep/gen_aa.html
    #
    # If there are any entries with ambiguity codes 4 or 5, we need to find 
    # the corresponding ambiguity detail loop
    #

    ambiguityKey= "Ambiguity_code"
    residKey="Seq_ID"
    hasAmbiguityCode= ambiguityKey in shiftCat

    residKey=None
    for key in "Seq_ID Comp_index_ID".split():
        if key in shiftCat:
            residKey = key
            break
        pass
    if not residKey:
        raise Exception("cannot identify column correspnding to residue number")
    
    if hasAmbiguityCode:
        hasHighAmbiguity=False
        for ambigCode in shiftCat[ambiguityKey]:
            if ambigCode == '4' or ambigCode == '5':
                hasHighAmbiguity=True
                pass
            pass
    

        if hasHighAmbiguity:

            #
            # NMR-STAR shift ambiguity entries are a comma-delimited list of atom shift assign IDs.
            # Each is the assignment ID for a chemical shift assignment that has been given 
            # an ambiguity code of 4 or 5.  Each set indicates that the observed chemical shifts
            # are related to the defined atoms, but have not been assigned uniquely to a specific
            # atom in the set.
            # 

            raise Exception("need to implement high-ambiguity chemical shift "
                            "assignments")
        pass
    else:
        useAmbiguityCodes=False
        ambiguityKey= "ID"
        pass    
   
    #	set rawAmbiguityData [readNMRSTARLoop $inUnit [list "_Atom-shift_assign_ID_ambiguity"]]
    #
    #	set ambiguitySets [list]
    #	foreach elem $rawAmbiguityData {
    #	    lappend ambiguitySets [split $elem ","]

    shiftData=[]
    for (assignID,
         resNum,
         resName,
         atomName,
         shift,
         ambiguity) in shiftCat.asTable(("ID",
                                         residKey,
                                         "Comp_ID",
                                         "Atom_ID",
                                         "Val",
                                         ambiguityKey )):
        resNum=int(resNum)+residOffset
        if len(resName)==1:
            import selectTools
            resName = selectTools.renameResidues(resName)[0]
            pass
        
	#
	# ambiguity codes are often wrong, so default behavior is 
	# to make everything unique and let the standard 
	# non-stereo processing handle it
	#
        
        if not useAmbiguityCodes: 
            ambiguity='1'
            pass

        if ambiguity=='1':
            sels=[correctedSelection(resNum,resName,atomName)]
        elif ambiguity=='2':
            sels=geminalAmbigSels(resNum,resName,atomName)
        elif ambiguity=='3':
            sels=aromaticAmbigSels(resNum,resName,atomName)
        else:
            raise KeyError("Can't yet deal with ambiguity type " +
                           ambiguity)
        
        if segmentName!=None:
            sels = ["segid %s and (%s)" %(segmentName,sel) for sel in sels]
            pass
        entry =  [shift,sels]

        shiftData.append(entry)
        pass

    if verbose:
        print("%d shifts were read from NMR-STAR data." % len(shiftData))
        pass
    #convert to TCL string
    if tclOutput:
        ret=""
        for val,sels in shiftData:
            ret += "{" + val
            # PASD selections want outer parentheses
            for sel in sels:
                ret += ' {(' + sel + ')}'
                pass
            ret += "} "
            pass
        ret+=""
    else:
        ret = shiftData
        pass
    
    return ret

def nefShifts(nef,
              name=None,
              verbose=True,
              tclOutput=False):
    """
    Read chemical shift data from the global NEF <m cif>.Cif
    object, which must have been previously populated by calling readNEF().

    """

    from nefTools import shiftsPrefix, fromNefAtomname, getBlock

    table = getBlock(nef,"shifts",name)

    name = table.nef_chemical_shift_list.sf_framecode[0][len(shiftsPrefix):]

    if verbose: print("nefShifts: reading table:", name)

# field no longer supported
#    if table.nef_chemical_shift_list.atom_chemical_shift_units[0] != "ppm":
#        raise Exception("nefShifts: Only ppm units supported")
#
    shiftCat = table.nef_chemical_shift

    shiftData=[]
    for (segid,
         resid,
         resName,
         atomName,
         value,
         error) in shiftCat.asTable(("chain_code",
                                     "sequence_code",    
                                     "residue_name",     
                                     "atom_name",        
                                     "value",            
                                     "value_uncertainty",
                                     )):
        try:
            resid=int(resid)
        except ValueError:
            print("nefShifts: ignoring shift entry with residue id", resid)
            continue
#        if len(resName)==1:
#            import selectTools
#            resName = selectTools.renameResidues(resName)[0]
#            pass
        
        segidStr = segid
        try:
            sel = "resid %d and resname %s and name %s" % (resid,
                                                           resName,
                                                           fromNefAtomname(atomName))
        except ValueError:
            print("nefShifts: ignoring shift entry with atomname", atomName)
            continue
            
        if not segid in ("","."):
            sel = "segid " + segid + " and " + sel
            pass
        entry =  [value, sel]
                                              

        shiftData.append(entry)
        pass

    if verbose:
        print("%d shifts were read from NEF data." % len(shiftData))
        pass
    #convert to TCL string
    if tclOutput:
        ret=""
        for entry in shiftData:
            ret += "{" + entry[0]
            for sel in entry[1:]:
                ret += ' {' + sel + '}'
                pass
            ret += "} "
            pass
        ret+=""
    else:
        ret = shiftData
        pass
    
    return ret


def nefPeaks(nef,
             noePot,
             fromProtonDim=None,
             toProtonDim=None,
             fromHeavyatomDim=None,
             toHeavyatomDim=None,
             name=None,
             verbose=True,
             tclOutput=False):
    """
    Read NEF peak table.

    The from/toProtonDim and from/toHeavyatomDim are used to
    manually specify appropriate columns for chemical shift values.

    Specify the name of the NEF saveSet where the experimental data is located,
    and the identities of labels for the dimensions of the experiment.
    """

    from nefTools import peaksPrefix

    sname=name #saveset name

    snames = [key for key in list(nef.keys()) if key.startswith(peaksPrefix)]
    if not sname:
        sname=snames[0]
    else:
        sname = peaksPrefix + sname
        pass

    if not sname in snames:
        raise Exception("nefPeaks: could not find saveset named ",
                        sname)
    name=sname[len(peaksPrefix):]

    if verbose:
        print("nefPeaks: reading peak information from " + name)
        pass

    peakTable = nef[sname]

    dimMap={}
    if fromProtonDim!=None:
        numDims = 2
        if fromHeavyatomDim!=None:
            numDims+=1
            pass
        if toHeavyatomDim!=None:
            numDims+=1
            pass
    
        if numDims != int(peakTable.nef_nmr_spectrum.num_dimensions[0]):
            raise Exception("nefPeaks: specified dimenensions not consistent " +
                            "with actual number of dimensions")
        pass
    else:  #auto-determine dimensions
        numDims = int(peakTable.nef_nmr_spectrum.num_dimensions[0])
        if numDims<2 or numDims>4:
            raise Exception("nefPeaks: cannot handle this number of "
                            "dimensions: %d" % numDims)
        fromHeavyatomDim=None
        toHeavyatomDim=None
        if numDims==2:
            fromProtonDim=1
            toProtonDim=1
            dimMap[fromProtonDim] = "fromProton"
            dimMap[toProtonDim] = "toProton"
            pass
        else: #dim>2
            #FIX: will not work for non-proton detected experiments
            dimAtoms=dict((i+1,code[-1]) for i,code in 
                          enumerate(peakTable.nef_spectrum_dimension.axis_code))
            protDims=[t for t in dimAtoms if dimAtoms[t]=='H']
            
            for dim1,dim2,type in zip(
                peakTable.nef_spectrum_dimension_transfer.dimension_1,
                peakTable.nef_spectrum_dimension_transfer.dimension_2,
                peakTable.nef_spectrum_dimension_transfer.transfer_type):
                dim1=int(dim1)
                dim2=int(dim2)
                if type=="onebond":
                    if fromHeavyatomDim==None:
                        if dim1 in protDims:
                            fromProtonDim=dim1
                            fromHeavyatomDim=dim2
                        else:
                            fromProtonDim=dim2
                            fromHeavyatomDim=dim1
                            pass
                        toProtonDim=protDims[0] \
                                     if fromProtonDim == protDims[1] \
                                     else protDims[1]
                        pass
                    else: # 4D
                        if dim1 in protDims:
                            toProtonDim=dim1
                            toHeavyatomDim=dim2
                        else:
                            toProtonDim=dim2
                            toHeavyatomDim=dim1
                            pass
                        pass
                    pass
                elif type=="through-space":
                    if not dim1 in protDims and dim2 in protDims:
                        raise Exception("nefPeaks: do not understand a "
                                        "through-space "
                                        "nef_spectrum_dimension_transfer "
                                        "entry")
                    pass
                else:
                    raise Exception("nefPeaks: do not understand the " +
                                    "{} ".format(type) +
                                    "nef_spectrum_dimension_transfer " +
                                    "entry")
                pass
            pass

        if verbose:
            print("  auto-detected a {}-D spectrum ".format(numDims) +
                  "with the following dimension identities:")
            print("    from Proton:    {}".format(fromProtonDim))
            print("    from Heavyatom: {}".format(fromHeavyatomDim))
            print("    to Proton:      {}".format(toProtonDim))
            print("    to Heavyatom:   {}".format(toHeavyatomDim))
            pass

        pass

    dimMap[fromProtonDim] = "fromProton"
    dimMap[toProtonDim] = "toProton"
    if fromHeavyatomDim!=None: dimMap[fromHeavyatomDim] = "fromHeavyatom"
    if toHeavyatomDim!=None: dimMap[toHeavyatomDim] = "toHeavyatom"
    
    peakCat = peakTable.nef_peak

    diffTol=1e-6 # peaks can be listed multiple times to itemize different
                 #possible assignments. However, the peak intensities must
                 #be the same to within diffTol
    from pasdPeak import Peak
    for i,peakID in enumerate(peakCat.peak_id):
        peak = Peak("%s%s" % (name,peakID))
        peak.setIntensity( float(peakCat.volume[i]) )
        peak.setFromProtonShift( float(peakCat["position_%d" %
                                               fromProtonDim][i]) )
        peak.setToProtonShift( float(peakCat["position_%d" %
                                             toProtonDim][i]) )
        if fromHeavyatomDim!=None:
            peak.setFromHeavyatomShift( float(peakCat["position_%d" %
                                                      fromHeavyatomDim][i]))
            pass
        if toHeavyatomDim!=None:
            peak.setToHeavyatomShift( float(peakCat["position_%d" %
                                                    toHeavyatomDim][i]) )
            pass
        peak.appendToNote("from NEF table %s peak %s" % (name,peakID))
        if noePot.hasPeakNamed( peak.name() ):
            if abs(noePot.peakNamed(peak.name()).intensity()-
                   peak.intensity()                           )>diffTol:
                raise Exception("nefPeaks: found two peaks with name %s, but" +
                                "with different intensities")
            pass
        else:
            
            peak.thisown=False # will be cleaned-up by noePot
            noePot.addPeak( peak )
            pass
        pass

    if verbose:
        print("  read %d peaks" % noePot.numPeaks())
        pass

    ret={}
    ret['numDims'] = numDims
    axis_code = peakTable.nef_spectrum_dimension.axis_code
    dimAtoms=dict((i+1,code[-1]) for i,code in enumerate(axis_code))

    ret['fromProton'] = axis_code[ fromProtonDim-1 ] 
    ret['toProton'] = axis_code[ toProtonDim-1 ] 
    ret['fromHeavyatom'] = axis_code[
        fromHeavyatomDim-1] if fromHeavyatomDim else ""
    ret['toHeavyatom'] = axis_code[
        toHeavyatomDim-1] if toHeavyatomDim else ""
        
    dimInfo = peakTable.nef_spectrum_dimension
    for (id,
         units,
         axisName,
         spectralWidth,
         firstPoint,
         foldingType,
         foldingNotNeeded,
         ) in dimInfo.asTable(["dimension_id",
                               "axis_unit",
                               "axis_code",
                               "spectral_width",
                               "value_first_point",
                               "folding",
                               "absolute_peak_positions",
                               ]):
        foldingNeeded=True if foldingNotNeeded=="false" else False
        if units!="ppm":
            raise Exception("only support units of ppm")
        id=int(id)
        if id not in list(dimMap.keys()):
            raise Exception("nefPeaks: dimension %s not specified in arguments"%
                            id)

        if not foldingType in ('circular','none'):
            raise Exception("nefPeaks: folding type %s not supported" % foldingType)
            
        retKey = dimMap[id] + "SpectralRange"
        if foldingType=="none": foldingNeeded=False
            
        start = float(firstPoint) if foldingNeeded else -1e30
        end   = float(firstPoint) + \
                float(spectralWidth) if foldingNeeded else 1e30
        ret[retKey] = (start,end)
        if verbose:
            if foldingNeeded:
                spectralRangeString = "spectral range %f -> %f" % (start,end)
            else:
                spectralRangeString = "no folding"
            print("  %s Dimension has name: %s with %s" % (dimMap[id],
                                                           axisName,
                                                           spectralRangeString))
            pass
        pass
    
    if tclOutput:
        for key in list(ret.keys()):
            try:
                if len(ret[key])==2:
                    ret[key] = "%f %f" % ret[key]
                    pass
                pass
            except:
                ret[key] = "%s" % str(ret[key])
            pass
        pass

    return ret    
   
    
def pipeShifts(pipe,
               verbose=True,
               tclOutput=False):
    """Read chemical shift data from NMRPipe table format.

    pipe is a string specifying either the path of the table or, directly, its
    contents.

    """
    import os
    if os.path.isfile(pipe):  # is pipe filename?
        infile = open(pipe, 'r')
        pipe = infile.read()
        infile.close()
    
    sel_template = '(segid "%s" and resid %s and resname %s and name %s)'

    code = {'P': 'PRO', 'G': 'GLY', 'A': 'ALA', 'R': 'ARG', 'N': 'ASN',
            'D': 'ASP', 'C': 'CYS', 'Q': 'GLN', 'E': 'GLU', 'H': 'HIS',
            'I': 'ILE', 'L': 'LEU', 'K': 'LYS', 'M': 'MET', 'F': 'PHE',
            'S': 'SER', 'T': 'THR', 'W': 'TRP', 'Y': 'TYR', 'V': 'VAL'}

    shiftData = []

    # Read input table.
    lines = pipe.splitlines()
    for line in lines:
        
        # Dismiss empty and comment lines.
        if line and not line.strip().startswith('#'):  
            
            words = line.split()
            
            # Get column definitions.
            if words[0] == 'VARS':
                
                columns = words[1:]
                
                shift_idx = columns.index('SHIFT')
                resid_idx = columns.index('RESID')
                resname_idx = columns.index('RESNAME')
                name_idx = columns.index('ATOMNAME')

                segid_idx = None
                try:
                    segid_idx = columns.index('SEGNAME')
                except ValueError:
                    continue

                try:
                    segid_idx = columns.index('CHAINNAME')
                except ValueError:
                    continue
                
            # Line with data entry: line that doesn't start with special word.
            #                  Is list complete ??!!
            if words[0] not in ['VARS', 'FORMAT', 'REMARK', 'DATA']:

                shift = words[shift_idx]
                resid = words[resid_idx]
                resname = words[resname_idx]
                if len(resname) == 1: resname = code[resname] # 1->3-letter code
                name = words[name_idx]

                if segid_idx is None:
                    segid = ''
                else:
                    segid = words[segid_idx]

                sel = sel_template % (segid, resid, resname, name)
                
                shiftData.append((shift, sel))

    if verbose: print('%i shifts read from NMRPipe table' % len(shiftData))

    # Convert to TCL string (adapted from nefShifts function).
    if tclOutput:
        result = ''
        for entry in shiftData:
            result += '{' + entry[0]
            for sel in entry[1:]:
                result += ' {' + sel + '}'
            result += '} '
    else:
        result = shiftData

    return result          


# format for mappings is resType, nmrStarAtomName, xplorAtomName

mappings=[["*"  , "h"   , "hn"]   ,
          ["gly", "ha3" , "ha2"]  ,
          ["gly", "ha2" , "ha1"]  ,
          ["ile", "hg13", "hg12"] ,
          ["ile", "hg12", "hg11"] ,
          ["phe", "hb3" , "hb2"]  ,
          ["phe", "hb2" , "hb1"]  ,
          ["phe", "cd"  , "cd#"]  ,
          ["phe", "ce"  , "ce#"]  , 
          ["trp", "hb3" , "hb2"]  ,
          ["trp", "hb2" , "hb1"]  ,
          ["cys", "hb3" , "hb2"]  ,
          ["cys", "hb2" , "hb1"]  ,
          ["ser", "hb3" , "hb2"]  ,
          ["ser", "hb2" , "hb1"]  ,
          ["asn", "hb3" , "hb2"]  ,
          ["asn", "hb2" , "hb1"]  ,
          ["tyr", "hb3" , "hb2"]  ,
          ["tyr", "hb2" , "hb1"]  ,
          ["tyr", "cd"  , "cd#"]  ,
          ["tyr", "ce"  , "ce#"]  , 
          ["his", "hb3" , "hb2"]  ,
          ["his", "hb2" , "hb1"]  ,
          ["asp", "hb3" , "hb2"]  ,
          ["asp", "hb2" , "hb1"]  ,
          ["met", "hb3" , "hb2"]  ,
          ["met", "hb2" , "hb1"]  ,
          ["met", "hg3" , "hg2"]  ,
          ["met", "hg2" , "hg1"]  ,
          ["gln", "hb3" , "hb2"]  ,
          ["gln", "hb2" , "hb1"]  ,
          ["gln", "hg3" , "hg2"]  ,
          ["gln", "hg2" , "hg1"]  ,
          ["glu", "hb3" , "hb2"]  ,
          ["glu", "hb2" , "hb1"]  ,
          ["glu", "hg3" , "hg2"]  ,
          ["glu", "hg2" , "hg1"]  ,
          ["pro", "hb3" , "hb2"]  ,
          ["pro", "hb2" , "hb1"]  ,
          ["pro", "hg3" , "hg2"]  ,
          ["pro", "hg2" , "hg1"]  ,
          ["pro", "hd3" , "hd2"]  ,
          ["pro", "hd2" , "hd1"]  ,
          ["arg", "hb3" , "hb2"]  ,
          ["arg", "hb2" , "hb1"]  ,
          ["arg", "hg3" , "hg2"]  ,
          ["arg", "hg2" , "hg1"]  ,
          ["arg", "hd3" , "hd2"]  ,
          ["arg", "hd2" , "hd1"]  ,
          ["lys", "hb3" , "hb2"]  ,
          ["lys", "hb2" , "hb1"]  ,
          ["lys", "hg3" , "hg2"]  ,
          ["lys", "hg2" , "hg1"]  ,
          ["lys", "hd3" , "hd2"]  ,
          ["lys", "hd2" , "hd1"]  ,
          ["lys", "he3" , "he2"]  ,
          ["lys", "he2" , "he1"]  ,
          ["ala", "hb*" , "hb#"] ,
          ["val", "hg1*", "hg1#"],
          ["val", "hg2*", "hg2#"],
          ["val", "cg"  , "cg#"],
          ["ile", "hg2*", "hg2#"],
          ["ile", "hd*" , "hd*"]  ,
          ["leu", "hb3" , "hb2"]  ,
          ["leu", "hb2" , "hb1"]  ,
          ["leu", "hd1*", "hd1#"],
          ["leu", "hd2*", "hd2#"],
          ["leu", "cd"  , "cd#"],
          ["met", "he*" , "he*"]  ,
          ["thr", "hg2*", "hg2#"]]

def correctedSelection(resNum, resName, atomName):
    """

 NMR-STAR formatted shift tables typically have three 
 types of atom name problems:

 1.  Backbone amides are named H, rather than HN.  

 2.  methylenes are typically named H*2 and H*3 instead of 
 H*1 and H*2.  H*3 maps to H*2 and H*2 maps to H*1.  

 3.  Selections involving methyls usually only select one
 of the three protons.  

 Returns an xplor selection with the atom name corrected/expanded 
 as necessary.
 """
    for mResName,oldAtomName,newAtomName in mappings:
        if (nameMatches(mResName, resName) and
            nameMatches(oldAtomName, atomName)):
            atomName=newAtomName
            break
        pass

    return "(resid %d and resn %s and name %s)" % \
           (resNum, resName.lower(), atomName.lower())


def nameMatches(pattern,name):
    pattern=pattern.lower()
    pattern = pattern.replace("*",".*") + "$"
    name = name.lower()
    import re
    if re.match(pattern,name):
        return True
    return False
    

def aromaticAmbigSels(resNum, resName, atomName):
    

    if not resName in ("PHE","TYR"):
        raise Exception("can't parse aromatic ambiguity resid "
                        "%d resn %s atomname %s" % (resNum,resName, atomName))
    
    for match in ("hd*", "he*", "cd*", "ce*"):
        if nameMatches( match, atomName):
            return [correctedSelection(resNum,resName, match)]
        pass

    raise Exception("can't parse aromatic ambiguity resid "
                    "%d resn %d atomname %d" % (resNum, resName, atomName))

geminalAmbigTab= {
    "GLY" : [ "ha*" ],
    "ALA" : [],
    "VAL" : [ "hg*", "cg*" ],
    "ILE" : [ "hg1*", "cg1*" ],
    "LEU" : [ "hb*", "hd*", "cd*" ],
    "PRO" : [ "hb*", "hg*", "hd*" ],
    "THR" : [],
    "GLU" : [ "hb*", "hg*" ],
    "MET" : [ "hb*", "hg*" ],
    "GLN" : [ "hb*", "hg*", "he2*"],
    "ASN" : [ "hb*", "hd2*" ],
    "PHE" : [ "hb*" ],
    "TRP" : [ "hb*" ],
    "CYS" : [ "hb*" ],
    "SER" : [ "hb*" ],
    "TYR" : [ "hb*" ],
    "HIS" : [ "hb*" ],
    "ASP" : [ "hb*" ],
    "LYS" : [ "hb*", "hg*", "hd*", "he*"],
    "ARG" : [ "hb*", "hg*", "hd*"]
    }

def geminalAmbigSels(resNum, resName, atomName):

    for match in geminalAmbigTab[resName]:
        if nameMatches(match,atomName):
            return [correctedSelection(resNum,resName,match)]
        pass
    raise Exception("can't parse geminal ambiguity resid "
                    "%d resn %s atomname %s" % (resNum, resName, atomName))
    pass

    
def starPeaks(fromProtonColumnName,
              toProtonColumnName,
              fromHeavyatomColumnName=None,
              toHeavyatomColumnName=None,
              saveSet=None,
 #    set idName     [flagVal $args -peakIDcolumnName "PkID"]
#    set intName    [flagVal $args -intensityColumnName "Intensity"]
                tclOutput=True):
    """
    Read nmr-star table for 3d experiments.

    Specify NMR-STAR saveSet where the experimental data is located,
    and the identities of labels for the dimensions of the experiment.
    These names correspond to NMR-STAR items Spectral_dim.Dimension_name.
    """

    if not saveSet:
        saveSets=[]
        for key in list(starData.keys()):
            if "Spectral_dim" in list(starData[key].keys()):
                saveSets.append(key)
                pass
            pass
        if len(saveSets)==0:
            raise Exception("starPeaks: "
                            "could not locate any NOE Peak saveSets.")
        if len(saveSets)>1:
            raise Exception("starPeaks: " +
                            "too many NOE Peak saveSets. Please choose one of"+
                            str(saveSets))
        saveSet=saveSets[0]
        print("starPeaks: using saveSet:", saveSet)
    else:
        if saveSet not in list(starData.keys()):
            raise Exception("starPeaks: could not find saveSet named >%s<" %
                            saveSet)
        pass
        
        

    dimCat=starData[saveSet].Spectral_dim

    #FIX: spectral width units? Hz?
    for i in range(len(dimCat.ID)):
        dimName=dimCat.Dimension_name[i]
        if dimName==fromProtonColumnName:
            h1Dim=dimCat.Dataset_dimension[i]
            h1Width=dimCat.Sweep_width[i]

        elif dimName==toProtonColumnName:
            h2Dim=dimCat.Dataset_dimension[i]
            h2Width=dimCat.Sweep_width[i]
        elif dimName==fromHeavyatomColumnName:
            c1Dim=dimCat.Dataset_dimension[i]
            c1Width=dimCat.Sweep_width[i]
        elif dimName==toHeavyatomColumnName:
            c2Dim=dimCat.Dataset_dimension[i]
            c2Width=dimCat.Sweep_width[i]
        else:
            raise Exception("found unrecognized Dimension_name: " + dimName)
        pass


    peaks={}
    peakGeneralCat=starData[saveSet].Peak_general_char
    for (id,inten,meth) in peakGeneralCat.asTable(("Peak_ID",
                                                   'Intensity_val',
                                                   'Measurement_method')):
        peak={}
        if meth!="height":
            continue
        peak["id"]    = id
        peak["inten"] = inten
        peaks[id] = peak
        pass


    peakCharCat = starData[saveSet].Peak_char
    for (id,dim,shift) in peakCharCat.asTable(("Peak_ID",
                                               'Spectral_dim_ID',
                                               'Chem_shift_val')):
        if dim==h1Dim:
            peaks[id]["h1Shift"] = shift
        elif dim==h2Dim:
            peaks[id]["h2Shift"] = shift
        elif dim==c1Dim:
            peaks[id]["c1Shift"] = shift
        elif dim==c2Dim:
            peaks[id]["c2Shift"] = shift
        else:
            raise Exception("unexpected dimension value: " + dim)
        pass

    if tclOutput:
        ret=""
        for entry in list(peaks.values()):
            ret += "{" + entry["id"]
            ret += ' ' + entry["inten"]
            ret += ' ' + entry["h1Shift"]
            ret += ' ' + entry["h2Shift"]
            if fromHeavyatomColumnName:
                ret += ' ' + entry["c1Shift"]
                pass
            if toHeavyatomColumnName:
                ret += ' ' + entry["c2Shift"]
                pass
            ret += "} "
            pass
        pass
    else:
        ret=peaks
        pass

    return ret    
    
def sparkyPeaks(filename,
                noePot,
                fromProtonColumnName,
                toProtonColumnName,
                fromHeavyatomColumnName=None,
                toHeavyatomColumnName=None,
                name=None,
                verbose=True,
                tclOutput=True):
    """Read Sparky format into <m pasdPot>.PASDPot

    arguments:

      filename                - name of sparky file
      pasdPot                 - a PASDPot instance
      fromProtonColumnName    - string label for from- Proton
      toProtonColumnName      - string label for to- Proton
      fromHeavyatomColumnName - string label for from- heavyatom
      toHeavyatomColumnName   - string label for to- heavyatom
      name                    - prefix to use for peak names. In PASD,
                                these must be different for each peaklist.
      verbose                 - whether to produce verbose output
      tclOutput               - If True, return a valid TCL list string
    """
    lines = open(filename).readlines()
    columnNames=lines[0].split()

    namePrefix = name if name else noePot.instanceName()

    fromProtonColumn=None
    toProtonColumn=None
    fromHeavyatomColumn=None
    toHeavyatomColumn=None

    for i,name in enumerate(columnNames):
        if name == fromProtonColumnName    : fromProtonColumn = i
        if name == toProtonColumnName      : toProtonColumn = i
        if name == fromHeavyatomColumnName : fromHeavyatomColumn = i
        if name == toHeavyatomColumnName   : toHeavyatomColumn = i
        pass

    if fromProtonColumnName != None and fromProtonColumn == None:
        raise Exception("column name not found: %s" % fromProtonColumnName)
    if toProtonColumnName != None and toProtonColumn == None:
        raise Exception("column name not found: %s" % toProtonColumnName)
    if fromHeavyatomColumnName != None and fromHeavyatomColumn == None:
        raise Exception("column name not found: %s" % fromHeavyatomColumnName)
    if toHeavyatomColumnName != None and toHeavyatomColumn == None:
        raise Exception("column name not found: %s" % toHeavyatomColumnName)
    
    heightColumn = -1 if "Data Height" in lines[0] else None

    peaks={}
    peakNum=1

    diffTol=1e-6 # peaks can be listed multiple times to itemize different
                 #possible assignments. However, the peak intensities must
                 #be the same to within diffTol

    from pasdPeak import Peak
    for line in lines[2:]:

        peakInfo = {}
        peakInfo["id"] = str(peakNum)
        fields = line.split()
        if not fields:
            continue
        peakInfo["inten"] = fields[heightColumn] if heightColumn else "1"

        peakInfo["h1Shift"] = fields[ fromProtonColumn ]
        peakInfo["h2Shift"] = fields[ toProtonColumn ]
        if fromHeavyatomColumn!=None:
            peakInfo["c1Shift"] = fields[ fromHeavyatomColumn ]
            pass
        if toHeavyatomColumn!=None:
            peakInfo["c2Shift"] = fields[ toHeavyatomColumn ]
            pass
        peaks[id] = peakInfo

        peak = Peak("%s%s" % (namePrefix,peakInfo["id"]))
        peak.setIntensity( float(peakInfo["inten"]) )
        peak.setFromProtonShift( float(peakInfo["h1Shift"]) )
        peak.setToProtonShift( float(peakInfo["h2Shift"]) )
        if fromHeavyatomColumn != None:
            peak.setFromHeavyatomShift( float(peakInfo["c1Shift"]) )
            pass
        if toHeavyatomColumn != None:
            peak.setToHeavyatomShift( float(peakInfo["c2Shift"]) )
            pass
        peak.appendToNote("from set %s, peak %s" % (namePrefix,peakNum))
        if noePot.hasPeakNamed( peak.name() ):
            if abs(noePot.peakNamed(peak.name()).intensity()-
                   peak.intensity()                           )>diffTol:
                raise Exception("sparkyPeaks: " +
                                "found two peaks with name " +
                                peak.name() +
                                ", but" +
                                "with different intensities")
            pass
        else:
            
            peak.thisown=False # will be cleaned-up by noePot
            noePot.addPeak( peak )
            pass
        peakNum += 1
        pass

    if verbose:
        print("  read %d peaks" % noePot.numPeaks())
        pass
    
    if tclOutput:
        ret=""
        for entry in list(peaks.values()):
            ret += "{" + entry["id"]
            ret += ' ' + entry["inten"]
            ret += ' ' + entry["h1Shift"]
            ret += ' ' + entry["h2Shift"]
            if fromHeavyatomColumnName:
                ret += ' ' + entry["c1Shift"]
                pass
            if toHeavyatomColumnName:
                ret += ' ' + entry["c2Shift"]
                pass
            ret += "} "
            pass
        pass
    else:
        ret=peaks
        pass

    return ret

def xeasyShifts(filename,
                verbose=True,
                tclOutput=True,
                simulation=None,
                segmentName=None):
    """Read XEASY format chemical shifts.
    """
    lines = open(filename).readlines()

    from simulation import currentSimulation
    if not simulation: simulation = currentSimulation()

    nullValue = "999.000"
    shiftData = []
    from atomSel import AtomSel
    # assume column ids: ID shift val, err, atom name, residue num
    for line in lines:
        id, shift,err,name,resid = line.split()

        if shift==nullValue:
            continue

        entry = [shift]
        sel = "name %s and resid %d" % (name,int(resid))
        if segmentName!=None:
            sel += ' and segid "%s"' % segmentName
            pass
        if len(AtomSel(sel,simulation))<1:
               raise Exception("selection (%s) should select at least one atom." %
                               sel)
        entry.append(sel)

        shiftData.append(entry)
        pass

    if verbose:
        print("%d shifts were read from XEASY data." % len(shiftData))
        pass
    #convert to TCL string
    if tclOutput:
        ret=""
        for entry in shiftData:
            ret += "{" + entry[0]
            for sel in entry[1:]:
                ret += ' {' + sel + '}'
                pass
            ret += "} "
            pass
        ret+=""
    else:
        ret = shiftData
        pass
    
    return ret

def findUnassignedAtoms(shiftList,
                        flexibleRegion="",
                        tclOutput=False,
                        simulation=None,
                        expectedNonassigned="""
     (name o*) or 
     (name c) or 
     (name s*) or
     (resn thr and name hg1) or 
     (resn ser and name hg) or 
     (resn lys and (name hz* or name nz)) or 
     (resn tyr and (name hh or name cg or name cz)) or 
     (resn arg and (name he or name hh* or name cz or name nh* or name ne)) or
     (resn gln and name cd) or 
     (resn glu and name cd) or 
     (resn asn and name cg) or 
     (resn asp and name cg) or
     (resn phe and name cg) or 
     (resn trp and (name cg or name cd2 or name ce2)) or 
     (resn his and name cg) or
     (resn pro and name n) or
     (name ht*)""",
                        ):
    """

    Given a list of shift entries, find atoms with no assignment

    expectedNonassigned - default atoms which we don't expect to see
                          (exchangable protons, oxygens, sulfurs, or
                          heavyatoms with no attached protons)
    flexibleRegion      - selection specifying atoms which we don't expect
                          to see due to flexibility.
    tclOutput           - if True, return an informational string,
                          else return a named tuple with string member
                          message a list of atoms missingAtoms.
    """

    if not simulation:
        from simulation import currentSimulation
        simulation=currentSimulation()
        pass

    from atomSel import AtomSel
    from selectTools import convertToAtomSel
    expectedNonassigned = convertToAtomSel(expectedNonassigned,simulation)
    flexibleRegion = convertToAtomSel(flexibleRegion,simulation)

    # init
    isAssigned = [0]*simulation.numAtoms()

    # record each selected atom in each shift entry

    for val, selList in shiftList:
        if type(selList)!=type([]): selList = [selList]
        for sel in selList:
            for index in AtomSel(sel,simulation).indices():
                isAssigned[index] += 1
                pass
            pass
        pass

    # complain about any missing expected shifts

    unexpectedIndices=expectedNonassigned.indices() + flexibleRegion.indices()

    unexpectedIndices = sorted(set(unexpectedIndices))
    
    missingAtoms=[]

    for i in range(simulation.numAtoms()):
        if isAssigned[i]==0 and not i in unexpectedIndices:
            missingAtoms.append( simulation.atomByID(i) )
            pass
        pass
    

        msg  = "Of %d total atoms," % simulation.numAtoms()
        msg += " %d " % (simulation.numAtoms()-len(unexpectedIndices))
        msg += "atoms were expected to be assigned.\n"   
        msg += "   This does not include atoms: "
        msg += "%s\n" % expectedNonassigned.string()

        if len(flexibleRegion) != 0:
            msg += "   This also does not include atoms in the flexible "
            msg += "regions: %s\n" % flexibleRegion.string()
            pass

        msg += "Of them, %d atoms were missing assignments\n"%len(missingAtoms)

        if len(missingAtoms)>0:
            msg += "Atoms with missing selections are:\n"
            for atom in missingAtoms:
                msg += atom.string() + '\n'
                pass
            pass
        pass
    if tclOutput:
        ret = msg
    else:
        from collections import namedtuple
        RetType = namedtuple('MsgMissingAtoms',['message','missingAtoms'])
        ret = RetType(msg,missingAtoms)
        pass        

    return ret

def convertShiftsFromTCL(tclString):
    shifts=[]
    startIndex=0
    from parseTools import findNested
    while True:
        startIndex=tclString.find('{',startIndex)+1
        stopIndex = findNested('{','}',startIndex,tclString)
        if startIndex == stopIndex:
            break
        elem = tclString[startIndex:stopIndex].strip()
        val = float(elem.split()[0])
        selString= " ".join(elem.split()[1:])
        selStart=0
        sels=[]
        while True:
            selStart=selString.find('{',selStart)+1
            selStop = findNested('{','}',selStart,selString)
            if selStart==selStop:
                break
            sel = selString[selStart:selStop]
            sels.append(sel)
            selStart = selStop+1
            pass
        shifts.append((val,sels))
        startIndex = stopIndex+1
        pass
#    print "happy:", len(shifts), shifts
    return shifts


def createShiftAssignments(
    shifts
    ):
    """
    Convert all the shift list entries from strings to AtomSels

    The shifts argument is a list of tuples (value, sel).
   
    ONLY CONVERT THE FIRST SELECTION OF EACH ENTRY!  because with
    PIPP, I only get multiple selections from HB1 | HB2 style
    entries.  Since these always come in pairs, I end up doubling
    the number of shift assigns from these non-stereo methylene
    entries
   
    This will eventually break processing of more sophisticated
    shift tables' data!
   

   
    complain both here and in the report about empty atom selections
       """
    #FIX!

    return

def structureMetrics(coordsByFilename,
                       peaks,
                       pasdPots,
                       violCutoff=0.5,
                       completenessWeight=0.,
                       inverseBound=4.,
                       inverseMethylCorrection=0.,
                       verbose=True,
                       sim=None):
    """
    Given a list of <filename, atomPosArray> pairs from grabPDBfiles, 
    check each one's violations of a given set of peaks
    and return a list of the named tuple StructureMetrics for each structure.
    This list is sorted by the structure's score.
    """

    if not sim:
        from simulation import currentSimulation
        sim = currentSimulation()
        pass

    nViols   =[]
    ncs      =[]
    violRepts=[]

    for name,coords in list(coordsByFilename.items()):
        sim.setAtomPosArr( coords )
        numNonViolated = len([peak for peak in peaks if
                              peak.isAssigned() and
                              peak.lowestViolation()<violCutoff])
        
        [pot.updateActivation() for pot in pasdPots]
        totClose = sum([pot.numCloseShiftAssignmentPairs() for
                        pot in pasdPots])
        totCloseAcc = sum([pot.numCloseAccountedShiftAssignmentPairs() for
                           pot in pasdPots])

        nc = float(totCloseAcc) / totClose

        violRepts.append( (name,coords,numNonViolated,nc) )
        nViols.append( numNonViolated )
        ncs.append( nc )
        pass

    minNViols = min( nViols )
    maxNViols = max( nViols )

    nViolRange = max(maxNViols - minNViols,1)

    minNC = min( ncs )
    maxNC = max( ncs )

    ncRange = max(maxNC - minNC,1)

    ret=[]
    
    from collections import namedtuple
    StructureMetrics = namedtuple('StructureMetrics',
  "filename coords numViols nc violScore ncScore score".split())
    for name,coords,curNViols,curNC in violRepts:
        curViolScore = 1- float(curNViols - minNViols) / nViolRange

        curNCScore =  1 - float(curNC - minNC) / ncRange
        curScore = ((1 - completenessWeight) * curViolScore + 
                    completenessWeight * curNCScore)
	
        ret.append( StructureMetrics(name,coords,curNViols,curNC,
                                     curViolScore,curNCScore,curScore) )
        pass

    ret.sort(lambda x,y: cmp(x.score,y.score))

    if verbose:
        for entry in ret:
            print(entry)
            pass
        pass

    return ret

def asciiHistogram(data,
                   minVal=None,
                   maxVal=None,
                   numBins=20,
                   numRows=20,
                   title="Histogram",
                   numTrailingLines=2,
                   useIntegerBins=False):
    """
    Return a string containing an ASCII-art histogram.
    """
    if minVal==None: minVal = min(data)
    if maxVal==None: maxVal = max(data)

    from math import ceil, floor
    if useIntegerBins: numBins = int(ceil(maxVal) - floor(minVal)) + 1


    def whichBin(val):
        #
        # defend against zero binWidth arising from maxVal == minVal
        #
        binWidth = float(maxVal - minVal) / numBins if maxVal != minVal else 1

        from math import floor
        aBin =  min(max(0,floor((val - minVal) / binWidth)),numBins-1)
        return int(aBin)

    numInBin=[0]*numBins
    for datum in data:
        numInBin[ whichBin(datum) ] += 1
        pass

    if len(data)==0:
        return ""


    fracInBin = [float(num)/len(data) for num in numInBin]
    fracPerRow = 1.0 / numRows
    numCols=[round(frac/fracPerRow) for frac in fracInBin]

    #
    # draw the bottom axis
    #

    ret = "0.0 +" + "-"*numBins

    #
    # now draw the histo
    #

    for rowCnt in range(numRows):
        row = "    |" if rowCnt<numRows-1 else "1.0 |"
        
        for binCnt in range(numBins):
            row += "X" if numCols[binCnt]>rowCnt else " "
            pass
        ret = row+'\n'+ret
        pass
    
    #
    # add the vert axis label
    #
    ret = "fraction of\n%d datapoints\n%s" % (len(data),ret)

    #
    # add the title
    #
    nSpaces = (numBins - len(title)) // 2
    title = " "*nSpaces + title
    ret = title + "\n\n" + ret

    #
    # add labels to the bottom axis
    # Integer values need both the check and the 
    # rounding to avoid problems with eg., minVal = 4.0
    #
    minAsString = str(minVal)
    maxAsString = str(maxVal)
    nSpaces = max(1, numBins - len(minAsString) - len(maxAsString))

    botLine = "     %s" % minAsString
    botLine += " "*nSpaces
    botLine += maxAsString

    ret += "\n" + botLine + "\n"*numTrailingLines

    return ret

def peakAssignLikelihoodsFromStructs(structureData,
                                     peaks,
                                     violCutoff=0.5,
                                     sim=None):
    """#
# Given a list of <filename, atomPosArray> pairs of converged structures 
# and a list of peaks, 
# use the coords to determine the likelihood of each peakAssignment
#
# Could be replaced with something that would calculate previous likelihood
# of any pair of ShiftAssignments, whether they're used in a PeakAssignment
# or not.
#
"""
    
    if not sim:
        from simulation import currentSimulation
        sim = currentSimulation()
        pass

    #
    # init array of times each PeakAssignment is violated
    #
    timesViolated={}
    for peak in peaks:
        for pa in peak.peakAssignments():
            name ="%s--%s" % (peak.name(), pa.name())
            timesViolated[name]=0
            pass
        pass
    
    
    #
    # record how often each PeakAssignment is violated 
    # in the converged structures
    #
    fileCount=0
    for structureDatum in structureData:
        sim.setAtomPosArr(structureDatum.coords)

        for peak in peaks:
            upBound  = peak.upBound()
            lowBound = peak.lowBound()
            
            for pa in peak.peakAssignments():
                viol = pa.violation(upBound, lowBound)
                if viol > violCutoff:
                    name ="%s--%s" % (peak.name(), pa.name())
                    timesViolated[name] += 1
                    pass
                pass
            pass
        pass
    
    #
    # convert timesViolated to likelihood for each peakAssignment
    #
    numStructs = len(structureData)
    for peak in peaks:
        for pa in peak.peakAssignments():
            name ="%s--%s" % (peak.name(), pa.name())
            newL = (numStructs - timesViolated[name]) / numStructs
            pa.setPreviousLikelihood( newL )
            pass
        pass
    return

longRangeCut=6

def reportNOEprecision(pot,
                       heavyatoms,
                       highLikelihoodCut,
                       lowLikelihoodCut=None):
    """Return list of strings.
    """
    from selectTools import convertToAtomSel
    heavyatoms = convertToAtomSel(heavyatoms)

    if lowLikelihoodCut==None: lowLikelihoodCutoff=1.0-highLikelihoodCut

    #
    # extract the longrange likelihoods and plot their histogram
    #
    longRange= [peak for peak in pot.peaks() if peak.isLongRange(longRangeCut)]

    temp = [peak.prevLikelihood() for peak in longRange]

    ret=[]
    ret.append(asciiHistogram(temp,
                              minVal=0,maxVal=1,
                              title="Longrange peak likelihoods"))
    #
    # compute NOE precision (# high-likelihood peaks / residue)
    #
    nHighLlongRange=len([peak for peak in longRange if
                         peak.prevLikelihood() >= highLikelihoodCut])
    noePrec=float(nHighLlongRange)/ len(heavyatoms)

    #
    # compute NOE discrimination (fraction of LR peaks w/ high or
    # low likelihood)
    #
    nLowLlongRange = len([peak for peak in longRange if 
                          peak.prevLikelihood() <= lowLikelihoodCut])
    noeDisc = float(nHighLlongRange +
                     nLowLlongRange) / len(longRange) if len(longRange) else 0
   
    line1 = "Defined region of structure is %s" % heavyatoms.string()
    line2 = "Number of high-likelihood long range peaks/residue in "
    line2 += "defined region: %f" % noePrec
    ret.append( line1 + '\n' + line2 )
    ret.append("Long-range NOE discrimination: %f %%" % (100.0 * noeDisc))
    return ret


def reportNOEaccuracy(pot,
                      refStructFile,
                      violCutoff,
                      highLikeliCut,
                      lowLikeliCut=None):
    
    lowLikeliCut = 1. - highLikeliCut if lowLikeliCut==None else lowLikeliCut

    #
    # read in the reference coords
    #
    from . import protocol
    protocol.initCoords(refStructFile)
	
    retVal=[]

    #
    # select longrange peaks
    #
    longRange = [peak for peak in peaks if peak.isLongRange(longRangeCut)]
    nLR = len(longRange)

    #
    # Report fraction of high-likelihood longrange peaks that are good.
    #
    highLR=len([peak for peak in longRange if
                         peak.prevLikelihood() >= highLikelihoodCut])
    nHighLR = len(highLR)

    nHighLRGood = len([peak for peak in highLR
                       if peak.lowestViolation() < violCut])
    nHighLRBad = len([peak for peak in highLR
                      if peak.lowestViolation() > violCut])
    if nHighLR>0:
        fracHighLRGood = 100. * nHighLRGood / nHighLR
        fracHighLRBad  = 100. * nHighLRBad  / nHighLR
        
        line1="%d of %d long range peaks have likelihood >= %f" % (nHighLR,
                                                                   nLR,
                                                        highLikelihoodCutoff)
        line2="%f %% of them agree with reference structure %s to within %f A"\
               %(fracHighLRGood, refStructFile,violCutoff)
        line3="%f %% of them do not agree with " %fracHighLRBad
        line3 += "reference structure %s to within %f A" %(refeStructFile,
                                                           violCutoff)
        
        ret.append( line1+'\n'+line2+'\n'+line3 )
    else:
        ret.append("No long range peaks have likelihood > %f" %
                   highLikelihoodCutoff)
        pass
    

    #
    # Report fraction of low-likelihood longrange peaks that are bad.
    #

    lowLR=[peak for peak in longRange if
           peak.prevLikelihood() <= lowLikeliCutoff]
    nLowLR=len(lowLR)

    nLowLRGood=len([peak for peak in lowLR if peak.violation()<violCut])
    nLowLRBad=len(lowLR) - nLowLRGood

    if nLowLR > 0:
        fracLowLRGood= 100. * nLowLRGood / nLowLR
        fracLowLRBad = 100. * nLowLRBad  / nLowLR
        
        line1="%d of %d long range peaks have likelihood <= %f" %\
                   (nLowLR, nLR,lowLikelihoodCutoff)
        line2="%f %% of them agree with reference structure %s to within %f A"%\
                 (fracLowLRGood, referenceStructFile, violCutoff)
        line3="%f %% of them do not agree with "  % fracLowLRBad
        line3 += "reference structure %s to within %f A"%(fracLowLRBad,
                                                          refStructFile,
                                                          violCutoff)
        
        ret.append( line1+'\n'+line2+'\n'+line3 )
    else:
        ret.append( "No long range peaks have likelihood < %f" % lowLikeliCut )
        pass

    #
    # Report fraction of medium-likelihood longrange peaks that are good.
    #
    medLR=[peak for peak in longRange
           if (peak.prevLikelihood() > lowLikelihoodCutoff and
               peak.prevLikelihood() < highLikelihoodCutoff)]
    nMedLR=len(medLR)

    nMedLRGood=len([peak for peak in medLR if peak.violation()< violCutoff])
    nMedLRBad = len(medLR) - nMedLRGood

    if nMedLR>0:
        fracMedLRGood=100. * nMedLRGood / nMedLR
        fracMedLRBad =100. * nMedLRBad  / nMedLR
        
        line1="%d of %d long range peaks have likelihood between %f and  %f"%\
               (nMedLR, nLR, lowLikelihoodCutoff, highLikelihoodCutoff)
        line2="%f %% of them agree with reference structure %s to within %f A"%\
               (fracMedLRGood, referenceStructFile, violCutoff)
        line3="%f %% of them do not agree with reference " % fracMedLRBad
        line3+="structure %s to within %f A" % (refStructFile, violCutoff)
        
        ret.append( line1+'\n'+line2+'\n'+line3 )
    else:
        ret.append("No long range peaks have likelihood between %f and %f" %
                   (lowLikelihoodCutoff,highLikelihoodCutoff))
        pass

    #
    # Report fraction of good longrange peaks that are high-likelihood
    #
    goodLR=[peak for peak in longRange if peak.violation()<violCutoff]
    nGoodLR=len(goodLR)

    nGoodLRHigh=len([peak for peak in goodLR if
                     peak.prevLikelihood()>highLikelihoodCutoff])
    nGoodLRLow=len([peak for peak in goodLR if
                    peak.precLikelihood()<= lowLikelihoodCutoff])
    nGoodLRMed=nGoodLR - (nGoodLRHigh + nGoodLRLow)

    if nGoodLR>0:
        fracGoodLRHigh= 100. * nGoodLRHigh / nGoodLR
        fracGoodLRLow = 100. * nGoodLRLow  / nGoodLR
        fracGoodLRMed = 100. * nGoodLRMed  / nGoodLR
        
        line1="%d of %d long range peaks agree " %(nGoodLR,nLR)
        line1+= "with reference structure %s to within %f A" % (refeStructFile,
                                                                violCutoff)
        line2="%f %% of them have likelihood >= %f" %(fracGoodLRHigh,
                                                      highLikelihoodCutoff)
        line3="%f %% of them have likelihood <= %f" %(fracGoodLRLow,
                                                      lowLikelihoodCutoff)
        line4="%f %% of them have likelihood " % fracGoodLRMed
        line4 += "between %f and %f" % (lowLikelihoodCutoff,
                                        highLikelihoodCutoff)
        
        ret.append( line1+'\n'+line2+'\n'+line3+'\n'+line4 )
    else:
        ret.append("No long range peaks agree with reference structure " +
                   "%s to within %f A" % (refStructFile, violCutoff))
        pass

    #
    # Report fraction of bad longrange peaks that are low-likelihood
    #
    badLR=[peak for peak in longRange if peak.violation()>violCutoff]
    nBadLR=len(badLR)

    nBadLRHigh=len([peak for peak in badLR if
                    peak.prevLikelihood()>=highLikelihoodCutoff])
    nBadLRLow=len([peak for peak in badLR if
                   peak.prevLikelihood()< lowLikelihoodCutoff])
    nBadLRMed=nBadLR - (nBadLRHigh + nBadLRLow)

    if nBadLR>0:
        fracBadLRHigh = 100. * nBadLRHigh / nBadLR
        fracBadLRLow  = 100. * nBadLRLow  / nBadLR
        fracBadLRMed  = 100. * nBadLRMed  / nBadLR
        
        line1 = "%d of %d long range peaks do not " %(BadLR, nLR)
        line1+="agree with reference structure %s to within %f A" %\
                (refStructFile, violCutoff)
        line2 = "%f %% of them have likelihood <= %f" % (fracBadLRLow,
                                                         lowLikelihoodCutoff)
        line3 = "%f %% of them have likelihood >= %f" % (fracBadLRHigh,
                                                         highLikelihoodCutoff)
        line4 = "%f %% of them have likelihood " % fracBadLRMed
        line4 += "between %f and %f" % (lowLikelihoodCutoff,
                                        highLikelihoodCutoff)
        
        ret.append( line1+'\n'+line2+'\n'+line3+'\n'+line4 )
    else:
        ret.append("No long range peaks do not agree with " +
                   "reference structure %s to within %f A" % (refStructFile,
                                                              violCutoff))
        pass


    #
    # simulate fbad-long w/ current likelihoods
    #
#FIX: add this back
#    ret.append("Frac bad long range forces = %f" %
#               [fracBadLRInfo -peakList [$pot peaks] -violCutoff $violCutoff -numIterations 10 -useLikelihoods]]
    
    #
    # select shortrange peaks
    #

    shortRange=[peak for peak in peaks if peak.isShortRange(lrCutoff)]
    nLR=len(shortRange)

    #
    # Report fraction of high-likelihood shortrange peaks that are good.
    #

    highLR = [peak for peak in shortRange if
              peak.prevLikelihood()>=highLikelihoodCutoff]
    nHighLR = len(highLR)

    nHighLRGood = len([peak for peak in highLR if peak.violation()<violCutoff])
    nHighLRBad  = len(highLR) - nHighLRGood

    if nHighLR>0:
        fracHighLRGood = 100. * nHighLRGood / nHighLR
        fracHighLRBad  = 100. * nHighLRBad  / nHighLR
        
        line1="%d of %d short range peaks have likelihood >= %f" % \
               (nHighLR, nLR, highLikelihoodCutoff)
        line2="%f %% of them agree with reference " % fracHighLRGood
        line2+="structure %s to within %f A" % (refStructFile, violCutoff)
        line3="%f %% of them do not agree with reference " % fracHighLRBad
        line3+="structure %s to within %f A" % (referenceStructFile,violCutoff)
        
        ret.append( line1+'\n'+line2+'\n'+line3 )
    else:
        ret.append("No short range peaks have likelihood > %f" %
                   highLikelihoodCutoff)
        pass
    #
    # Report fraction of low-likelihood shortrange peaks that are bad.
    #
    lowLR=[peak for peak in shortRange if peak<=lowLikelihoodCutoff]
    nLowLR=len(lowLR)
    
    nLowLRGood=len([peak for peak in lowLR if
                    peak.prevLikelihood()<violCutoff])
    nLowLRBad=len(lowLR) - nLowLRGood

    if nLowLR>0:
        fracLowLRGood= 100. * nLowLRGood / nLowLR
        fracLowLRBad = 100. * nLowLRBad  / nLowLR
        
        line1= "%d of %d short range peaks " % (LowLR,nLR)
        line1 +="have likelihood <= %f" %lowLikelihoodCutoff
        line2= "%f %% of them agree with reference " % fracLowLRGood
        line2+="structure %s to within %f A" % (refStructFile, violCutoff)
        line3= "%f %% of them do not agree with reference " % fracLowLRBad
        line3+="structure %s to within %f A" %(refStructFile,violCutoff)
        
        ret.append( line1+'\n'+line2+'\n'+line3 )
    else:
        ret.append("No short range peaks have likelihood < %f" %
                   lowLikelihoodCutoff)
        pass

    #
    # Report fraction of medium-likelihood shortrange peaks that are good.
    #
    medLR=[peak for peak in shortRange if
           (peak.prevLikelihood()>lowLikelihoodCutoff and
           peak.prevLikelihood()<highLikelihoodCutoff)]
    nMedLR=len(medLR)

    nMedLRGood=len([peak for peak in medLR if peak.violation()<violCutoff])
    nMedLRBad=len(medLR) - nMedLRGood

    if nMedLR>0:
        fracMedLRGood= 100. * nMedLRGood / nMedLR
        fracMedLRBad = 100. * nMedLRBad  / nMedLR
        
        line1="%d of %d short range peaks have " % (nMedLR, nLR)
        line1+="likelihood between %f and  %f" % (lowLikelihoodCutoff,
                                                  highLikelihoodCutoff)
        line2="%f %% of them agree with reference " % fracMedLRGood
        line2+="structure %s to within %f A" % (referenceStructFile,
                                                violCutoff)
        line3="%f %% of them do not agree with reference " % fracMedLRBad
        line3 += "structure %s to within %f A" % (referenceStructFile,
                                                  violCutoff)
        
        ret.append( line1+'\n'+line2+'\n'+line3 )
    else:
        ret.append("No short range peaks have likelihood between" +
                   "%f and %f" % (lowLikelihoodCutoff,highLikelihoodCutoff))
        pass

    #
    # Report fraction of good shortrange peaks that are high-likelihood
    #
    goodLR=[peak for peak in shortRange if peak.violation()<violCutoff]
    nGoodLR=len(goodLR)

    nGoodLRHigh=len([peak for peak in goodLR if
                     peak.prevLikelihood()>=highLikelihoodCutoff])
    nGoodLRLow =len([peak for peak in goodLR if
                     peak.prevLikelihood()<=lowLikelihoodCutoff])
    nGoodLRMed =nGoodLR - (nGoodLRHigh + nGoodLRLow)

    if nGoodLR>0:
        fracGoodLRHigh = 100. * nGoodLRHigh / nGoodLR
        fracGoodLRLow  = 100. * nGoodLRLow  / nGoodLR
        fracGoodLRMed  = 100. * nGoodLRMed  / nGoodLR
        
        line1 = "%d of %d short range peaks agree with " % (nGoodLR, nLR)
        line1+= "reference structure %s to within %f A" % (refStructFile,
                                                           violCutoff)
        line2="%f %% of them have likelihood >= %f"%(fracGoodLRHigh,
                                                     highLikelihoodCutoff)
        line3="%f %% of them have likelihood <= %f"%(fracGoodLRLow,
                                                     lowLikelihoodCutoff)
        line4="%f %% of them have likelihood " % fracGoodLRMed
        line4+="between %f and %f" % (lowLikelihoodCutoff,highLikelihoodCutoff)
        
        ret.append( line1+'\n'+line2+'\n'+line3+'\n'+line4 )
    else:
        ret.append("No short range peaks agree with reference structure " +
                   "%s to within %f A" % (referenceStructFile, violCutoff))
        pass

    #
    # Report fraction of bad shortrange peaks that are low-likelihood
    #
    badLR  = [peak for peak in shortRange if peak.violation()>violCutoff]
    nBadLR = len(badLR)

    nBadLRHigh =len([peak for peak in badLR if
                     peak.prevLikelihood()>=highLikelihoodCutoff])
    nBadLRLow  =len([peak for peak in badLR if
                     peak.prevLikelihood()<=lowLikelihoodCutoff])
    nBadLRMed  = nBadLR - (nBadLRHigh + nBadLRLow)

    if nBadLR>0:
        fracBadLRHigh = 100. * nBadLRHigh / nBadLR
        fracBadLRLow  = 100. * nBadLRLow  / nBadLR
        fracBadLRMed  = 100. * nBadLRMed  / nBadLR
        
        line1="%d of %d short range peaks do not " % (nBadLR, nLR)
        line1+="agree with reference structure %s to within %f A" % \
                (refStructFile, violCutoff)
        line2="%f %% of them have likelihood <= %f" %(fracBadLRLow,
                                                      lowLikelihoodCutoff)
        line3="%f %% of them have likelihood >= %f" %(fracBadLRHigh,
                                                      highLikelihoodCutoff)
        line4="%f %% of them have likelihood between %f and %f" %\
               (fracBadLRMed,  lowLikelihoodCutoff, highLikelihoodCutoff)
        
        ret.append( line1+'\n'+line2+'\n'+line3+'\n'+line4 )
    else:
        ret.append("No short range peaks do not agree with reference " +
                   "structure %s to within %f A" %(refStructFile,violCutoff))
        pass


    #
    # simulate fbad-short w/ current likelihoods
    #
    #FIX: add this:
    #lappend retVal [format "Frac bad short range forces = %f" [fracBadLRInfo -peakList [$pot peaks] -violCutoff $violCutoff -numIterations 10 -useLikelihoods]]
    
    return ret

def removeIntraSegidPAs(pasdPot,
                        segidSel):
    from selectTools import convertToAtomSel
    segidSel = convertToAtomSel(segidSel)

    cnt=0
    from atomSel import intersection
    for peak in pasdPot.peaks():
        for pa in peak.peakAssignments():
            if ( len(intersection(pa.fromProtonSelection(),segidSel)) and
                 len(intersection(pa.toProtonSelection(),segidSel))      ):
                peak.removePeakAssignmentNamed(pa.name())
                cnt += 1
                pass
            pass
        pass
    print("removeIntraSegidPAs: Removed %d peak assignments" % cnt)
    pass

def removeInterSegidPAs(pasdPot,
                        segidSel1,
                        segidSel2):
    from selectTools import convertToAtomSel
    segidSel1 = convertToAtomSel(segidSel1)
    segidSel2 = convertToAtomSel(segidSel2)

    cnt=0
    from atomSel import intersection
    for peak in pasdPot.peaks():
        for pa in peak.peakAssignments():
            if ( (len(intersection(pa.fromProtonSelection(),segidSel1)) and
                  len(intersection(pa.toProtonSelection(),segidSel2))      ) or
                 (len(intersection(pa.fromProtonSelection(),segidSel2)) and
                  len(intersection(pa.toProtonSelection(),segidSel1))      ) ):
                peak.removePeakAssignmentNamed(pa.name())
                cnt += 1
                pass
            pass
        pass
    print("removeInterSegidPAs: Removed %d peak assignments" % cnt)
    pass

def removeDiagPAs(pasdPot):
    cnt=0
    from atomSel import intersection
    for peak in pasdPot.peaks():
        for pa in peak.peakAssignments():
            if ( len(intersection(pa.fromProtonSelection(),
                                  pa.toProtonSelection())) ):
                peak.removePeakAssignmentNamed(pa.name())
                cnt += 1
                pass
            pass
        pass
    print("removeDiagPAs: Removed %d peak assignments" % cnt)
    return

def covalentNeighbors(atom,range):
    """ return list of atoms within range bonds of argument atom
    """
    ret=[atom]
    dist=0
    from atomSel import AtomSel
    while True:
        if dist>=range:
            break
        dist += 1
        more=[]
        for atom in ret:
            more += AtomSel('bondedto atom "%s" %d %s' % (atom.segmentName(),
                                                          atom.residueNum(),
                                                          atom.atomName()))
            pass
        ret += more
        pass
    return ret

def hasAssignmentWithinBondedNeighborhood(peak,range):
    """ Return True if the specified peak contains any assignments containing
    a from/to pair separated by range bonds or fewer.
    """
    for pa in peak.peakAssignments():
        for atom in pa.fromAssignment().protonSelection():
            for covalentNeighbor in covalentNeighbors(atom,range):
                if covalentNeighbor in pa.toAssignment().protonSelection():
                    return True
                pass
            pass
        pass
    return False

from trace import notrace_decorate
@notrace_decorate
def removeBondedPeaks(pot,bondRange=2):
    """Remove all peaks from the < pasdPot>.PASDPot argument pot which
    contain any assignments with one or more from/to pair separated by
    bondRange or few bonds.
    """
    numRemoved=0
    for peak in pot.peaks():
        if hasAssignmentWithinBondedNeighborhood(peak,bondRange):
            pot.removePeakNamed( peak.name() )
            numRemoved += 1
            break
        pass
    print("removeBondedPeaks: removed %d peaks with" % numRemoved, end=' ')
    print("assignments between atoms < %d bonds distant" %(bondRange+1))
    return
