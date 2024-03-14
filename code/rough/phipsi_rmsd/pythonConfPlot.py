def plotStdVectors(*list):
    import numpy as np
    import matplotlib.pyplot as plt
    print "Printing from Python in plotStdVectors()"
    
    title = list[0]
    plt.figure(figsize=(18, 10), dpi=1000)
    dummylabels = -1
    for ensemble in list[1:]:        
        if isinstance(ensemble, basestring):            
            if (ensemble == "phi"):
                count=0
                ax1 = plt.subplot(2,1,1)                
                plt.title("Confplot for torsional angles (PDB ID:"+title+")")
                plt.ylabel(ensemble)
                plt.ylim(-180,180)                       
                plt.setp(ax1.get_xticklabels(),visible=False)
            elif (ensemble == "psi"):
                count=0
                ax2=plt.subplot(2,1,2)                     
                plt.ylabel(ensemble)
                plt.ylim(-180,180)
                plt.xlabel("residues")            
            else:
                count=0
                dummylabels = populateXdummy(ensemble, len(x))
        else:
            count=count+1            
            x = range(1,len(ensemble)+1)
            plt.plot(x,ensemble,'-o',label=count)        

# Option 2
# TkAgg backend
    #manager = plt.get_current_fig_manager()
    #manager.resize(*manager.window.maxsize())

    #plt.subplots_adjust(wspace=0.02, hspace=0.01, top=0.98, bottom=0.05, left=0.04, right=0.98)
    #ax1.legend(loc='upper center')

    if (dummylabels != -1):
        ax3 = ax2.twiny()
        ax3.xaxis.set_ticks_position("bottom")
        ax3.xaxis.set_label_position("bottom")

        # Offset the twin axis below the host
        ax3.spines["bottom"].set_position(("axes", -0.2))

        # Turn on the frame for the twin axis, but then hide all
        # but the bottom spine
        ax3.set_frame_on(True)
        ax3.patch.set_visible(False)

        ax3.set_xticks(x)
        ax3.set_xticklabels(dummylabels)
        ax3.set_xlabel("secondary structure")
        
        for tick in ax3.xaxis.get_major_ticks():
            tick.label.set_fontsize(8)

# Option 2
# TkAgg backend
    manager = plt.get_current_fig_manager()
    manager.resize(*manager.window.maxsize())

    plt.subplots_adjust(wspace=0.02, hspace=0.01, top=0.98, bottom=0.05, left=0.04, right=0.98)
    ax1.legend(loc='upper center')

    plt.tight_layout()
    #plt.show()
    ax1.set_xlim(xmin = 1) 
    ax2.set_xlim(xmin = 1) 
    figure_name = title+'.eps'
    #plt.show()
    plt.savefig(figure_name, format='eps', dpi=1000)
    
    return 0

def populateXdummy(filename, n_resi):
    labels = [""] * n_resi
    count=0;
    try:
        for line in open(filename).xreadlines():
            labels[count] = line.strip('\n')
            labels[count] = labels[count].replace('\t','\n')
            count = count+1
    except IOError:
        print "Could not read file ", filename
        return -1

    return labels
