"""
 Tools to help with configuration and analysis of SPARTA+ (Shifts
   Prediction from Analogue of Residue type and Torsion Angle).
 """

def create_Sparta(selection="all",
                  annSet="20110211",
                  verbose=False):
    """
    Create a SPARTA object given an atom selection.
    The annSet argument specifies the neural network model used.
    Valid choices are
       '20110211' - published SPARTA model with these modifications:
                      smoothed hbond, electric field contributions
                      order parameter term removed.
       '20110225' - published SPARTA model with these modifications:
                      smoothed electric field contribution, order parameter
                      term removed, hbond term replaced with sum of all
                      (smoothed) h-bond energies. This model does not have
                      discontinuities due to the change of identity of the
                      strongest H-bond (present in the 20110211 model).
       '20110324' - the 20110225 model but with the 20-dimensional BLOSUM
                      residue identification metric with a 5-dimensional
                      version introduced here:
                      http://www.pnas.org/content/102/18/6395.long

    If True the verbose argument causes informational messages to be printed.
    """
    from selectTools import convertToAtomSel
    selection = convertToAtomSel(selection)

    from atomSel import AtomSel, intersection
    selection = intersection(selection,AtomSel("not PSEUDO"))

    from sparta import rc_ptr_SPARTA, SPARTA
    rawPointer=SPARTA(selection)
    rawPointer.this.disown()
    s=rc_ptr_SPARTA(rawPointer)

    s.setVerbose( verbose )
    
    from os.path import join
    import os
    spartaPref=os.environ['SPARTA']
    if annSet=='20110211' or annSet=='20110225' or annSet=='20110324':
        s.setAnnSet(annSet)
        s.setANNLevel1Pattern(join(spartaPref,annSet,"ATOM.level1.PARAM.tab"))
    else:
        raise Exception("bad value for annSet: " + str(annSet))

    s.setTripFileName(os.path.join(spartaPref,"sparta.tab"))
    s.setWeightFileName(os.path.join(os.environ['SPARTA'],"weight.tab"))
    s.setHomoFileName(os.path.join(os.environ['SPARTA'],"homology.tab"))
    s.setFitFileName(os.path.join(os.environ['SPARTA'],"fitting.tab"))
    s.setRCFileName(os.path.join(os.environ['SPARTA'],"randcoil.tab"))
    s.setAdjFileName(os.path.join(os.environ['SPARTA'],"rcadj.tab"))
    s.setPrevFileName(os.path.join(os.environ['SPARTA'],"rcprev.tab"))
    s.setNextFileName(os.path.join(os.environ['SPARTA'],"rcnext.tab"))
    s.setB62FileName(join(os.environ['SPARTA'],annSet,"BLOSUM62.tab"))
    s.setSurfPattern(os.path.join(os.environ['SPARTA'],"errorSurface",
                                  "ATOM","AA..A450.S5.RMS.tab"));
    # ATOM -> atom name
    # PARAM -> 'WI', 'WL1', 'WL2'
#    if annSet=='20110225':
#        s.setAnnInputDim(95)
#        s.setAnnHiddenDim(20)
#        pass
#    if annSet=='20110324':
#        s.setAnnInputDim(50)
#        s.setAnnHiddenDim(20)
#        pass
        

    return s
