"""
Routines to convert to/from Xplor-NIH/IUPAC atom naming
"""

#map IUPAC --> Xplor-NIH
fromRestype={} #FIX: check these
fromRestype['ARG']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     'HG3': 'HG2' , 'HG2': 'HG1' ,
                     'HD3': 'HD2' , 'HD2': 'HD1' ,
                     }
fromRestype['ASN']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     }
fromRestype['ASP']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     }
fromRestype['CYS']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     }
fromRestype['GLN']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     'HG3': 'HG2' , 'HG2': 'HG1' ,
                     }
fromRestype['GLU']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     'HG3': 'HG2' , 'HG2': 'HG1' ,
                     }
fromRestype['GLY']={ 'HA3': 'HA2' , 'HA2': 'HA1' ,
                     }
fromRestype['HIS']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     'HG3': 'HG2' , 'HG2': 'HG1' ,
                     }
fromRestype['ILE'] = {'HG13': 'HG12' , 'HG12': 'HG11' ,
                      }
fromRestype['LEU']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     }
fromRestype['LYS']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     'HG3': 'HG2' , 'HG2': 'HG1' ,
                     'HD3': 'HD2' , 'HD2': 'HD1' ,
                     'HE3': 'HE2' , 'HE2': 'HE1' ,
                     }
fromRestype['MET']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     'HG3': 'HG2' , 'HG2': 'HG1' ,
                     }
fromRestype['PHE']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     }
fromRestype['PRO']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     'HG3': 'HG2' , 'HG2': 'HG1' ,
                     'HD3': 'HD2' , 'HD2': 'HD1' ,
                     }
fromRestype['SER']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     }
fromRestype['TRP']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     }
fromRestype['TYR']={ 'HB3': 'HB2' , 'HB2': 'HB1' ,
                     }

toRestype={}
for restype in list(fromRestype.keys()):
    toRestype[restype]={}
    for fromname in list(fromRestype[restype].keys()):
        toRestype[restype][ fromRestype[restype][fromname] ] = fromname
        pass
    pass

def toIUPAC(selection="all"):
    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)
    for atom in selection:
        name = atom.atomName()
        if name=="HN":
            name = "H"
        else:
            resname = atom.residueName()
            try:
                name = toRestype[resname][name]
            except KeyError:
                pass
            pass
        atom.setAtomName( name )
        pass
    return

def fromIUPAC(selection="all"):
    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)
    for atom in selection:
        name = atom.atomName()
        if name=="H":
            name = "HN"
        else:
            resname = atom.residueName()
            try:
                name = fromRestype[resname][name]
            except KeyError:
                pass
            pass
        atom.setAtomName( name )
        pass
    return

