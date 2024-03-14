"""code used to downcast rc_Pot objects to the correct specific Pot

this is lowlevel code which should not be normally used in scripts.
"""

def createDerivedPot(pot):
    """given a <m potList>.rc_Pot, cast the result to the correct type.
	This function works for all potential terms created within the Python
        interface."""
    try:
        return pot.instanceData()
    except:
        raise Exception("could not downcast term: (%s, %s)" %
                        (pot.potName(), pot.instanceName()))
    pass
