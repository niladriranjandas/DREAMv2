def plotStdVectors(*list):
    import numpy as np
    import matplotlib.pyplot as plt
    print "Printing from Python in plotStdVectors()"
    
    plt.ylim(-180,180)
    plt.title("Confplot for torsional angle")
    plt.xlabel("residues")
    for ensemble in list:        
        if isinstance(ensemble, basestring):
            print "string here"
            print ensemble
            plt.ylabel(ensemble)
        else:
            x = range(1,len(ensemble)+1)
            plt.plot(x,ensemble,'-o')
        
    plt.show()
    
    return 0
