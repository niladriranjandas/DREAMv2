
"""
Field definitions of ATOM and other records from the PDB format
"""

atomCols={'serial'     : ( 6,12),
          'name'       : (12,16),
          'altLoc'     : (16,17),
          'resName'    : (17,21),
          'chainID'    : (21,22),
          'resSeq'     : (22,26),   #resid
          'iCode'      : (26,27),          
          'x'          : (30,38),
          'y'          : (38,46),
          'z'          : (46,54),
          'occupancy'  : (54,60),
          'tempFactor' : (60,66),
          'segID'      : (72,76),
          'element'    : (76,78),
          'charge'     : (78,80)
          }

def getAtomEntry(name,record,ignoreShortEntries=False):
    try:
        return record[ slice(*atomCols[name]) ]
    except KeyError:
        print("Valid entry names: " + " ".join(list(atomCols.keys())))
        raise
    except IndexError:
        if ignoreShortEntries:
            return ""
        else:
            print("record is too short")
            raise
        pass
    return ""

#def insertAtomEntries(record):
#    """
#    Insert all PDB entries for the given ATOM record into the namespace
#    of the calling function. Note: if any of these values are already present
#    they will be overridden.
#    """
#    import inspect
#
#    def calling_scope_variable(name):
#        frame = inspect.stack()[1][0]
#        
#        while name not in frame.f_locals:
#            frame = frame.f_back
#            if frame is None:
#                return None
#            return frame.f_locals[name]    
#
