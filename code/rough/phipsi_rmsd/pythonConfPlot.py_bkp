def plotStdVectors(*list):
    import numpy as np
    import matplotlib.pyplot as plt
    print "Printing from Python in plotStdVectors()"
    
    title = list[0]
    plt.figure(figsize=(17, 7), dpi=600)
    for ensemble in list:        
        if isinstance(ensemble, basestring):            
            if (ensemble == "phi"):
                count=0
                ax1 = plt.subplot(2,1,1)                
                plt.title("Confplot for torsional angles (PDB ID:"+title+")")
                plt.ylabel(ensemble)
                plt.ylim(-180,180)                       
                plt.setp(ax1.get_xticklabels(),visible=False)
            if (ensemble == "psi"):
                count=0
                ax2=plt.subplot(2,1,2)                     
                plt.ylabel(ensemble)
                plt.ylim(-180,180)
                plt.xlabel("residues")            
        else:
            count=count+1            
            x = range(1,len(ensemble)+1)
            plt.plot(x,ensemble,'-o',label=count)        
    
# Option 2
# TkAgg backend
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())

    plt.subplots_adjust(wspace=0.02, hspace=0.01, top=0.98, bottom=0.05, left=0.04, right=0.98)
    ax1.legend(loc='upper center')

    #plt.show()
    plt.tight_layout()
    #plt.show()
    ax1.set_xlim(xmin = 1) 
    ax2.set_xlim(xmin = 1) 
    figure_name = title+'.eps'
    plt.savefig(figure_name, format='eps', dpi=1000)
    
    return 0
