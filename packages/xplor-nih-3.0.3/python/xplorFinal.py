
#
# and EnsembleSimulation has been shut down- don't need to worry about
# multiple (smp) threads

import xplor
from functools import reduce

##the rest is to collect and print CPU info at exit

xplor.exitMessage = ""
if xplor.exitStatus!=0:
    extra=''
    mess = "process %d%s: irregular exit.\n" % (xplor.p_processID,extra)

    if xplor.p_numProcs>1:
        xplor.exitMessage = mess
        pass
    pass
else:

    from ensembleSimulation import EnsembleSimulation_currentSimulation
    esim = EnsembleSimulation_currentSimulation()
    if esim:
        esim.barrier()
        pass
        
    if xplor.p_numProcs>1:

        state = xplor.disableOutput()
        (startTime,startDate,
         time,date          ) = xplor.command("""eval ($time=$time)
                                                 eval ($date=$date)""",
                                              ("start_time","start_date",
                                               "time","date")            )
        xplor.enableOutput(state)

        if esim:
            from ensembleSimulation import Comm
            comm=Comm(xplor.p_comm,esim)
        else:
            comm=xplor.p_comm
            pass
    
        cpu = simWorld.cpuTime()
    
        comm.barrier()
        
        cpus = comm.collect(cpu)
            
        if xplor.p_processID==0:
            from ensembleSimulation import singleThread, multiThread
            if singleThread():
                totCPU=reduce(lambda x,y: x+y, cpus)
                xplor.exitMessage  = "  total CPU time= %f\n" % totCPU
                xplor.exitMessage += "  entry time at   %s %s\n" % (startTime,
                                                                    startDate)
                xplor.exitMessage += "  exit  time at   %s %s\n" % (time, date)
                pass
            multiThread()
            pass
        pass
    pass
