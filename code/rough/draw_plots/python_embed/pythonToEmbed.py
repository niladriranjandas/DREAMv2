def plotStdVectors(*list):
    import numpy as np
    import matplotlib.pyplot as plt
    print "Printing from Python in plotStdVectors()"
    print list[0]
    print len(list[0])
    print list[1]
    print len(list[1])
    plt.plot(list[0], list[1])
    plt.show()
    
    return 0
