"""Produce XPLOR DCD files with the <m ivm>.IVM objects.
"""


from xplorDCD import openFile, closeFile
import simulation

class TrajFile:
    """a trajectory object to be used with an <m ivm>.IVM object, used to
    write coordinates to a XPLOR-style DCD file during a trajectory run.
    Use the IVM setTrajectory method to write a given ivm
    trajectory. Example:

    from ivm import IVM
    integrator=IVM()

    from trajFile import TrajFile
    traj=TrajFile("traj.dcd",saveInterval=2)
    integrator.setTrajectory( traj)

    integrator.run()
    
    """
    def __init__(s,filename,
                 start=0,
                 saveInterval=1,
                 freeAtoms='all',
                 isFormatted=False,
                 printCoords=False,
                 fixXY=False,
                 read=False,
                 constantTime=False,
                 delta_t=0.02
                 ):
        """Specify the filename for the output trajectory
        file. freeAtoms is an <m atomSel>.AtomSel specifying which
        atoms are moving. start specifies at which frame trajectory
        writing commences, while saveInterval specifies how many
        frames to skip between writing. If isFormatted=True, an ASCII
        trajectory file is written. If printCoords=True, in addition
        to writing the trajectory file, all freeAtoms coordinates are
        printed at each specified saveInterval. If fixXY is set to true, the
        x- and y- coordinates of the entire trajectory are centered on their
        center of mass values.

        If constantTime=True, frames will be 
        interpolated and written at equal time increments such that the
          time = index * delta_t
        
        """
        from selectTools import convertToAtomSel
        s.fileName        = filename
        s.unit            =-1
        s.freeAtoms       = convertToAtomSel(freeAtoms)
        s.isFirst         = True
        s.scale           = 10000.
        s.offset          = 800 # 300 for vel
        # used on first write
        s.saveInterval    = saveInterval
        s.start           = start
        s.delta_t         = delta_t
        s.header          = ""
        s.isFormatted     = isFormatted
        s.printCoords     = printCoords
        s.fixXY           = fixXY
        s.constantTime    = constantTime
        s.nextIndex       = start

        s.create(read)
        return
    def __del__(s):
        s.close()
        return

    def create(s,read):
        if s.fileName == "" or s.isFirst==False: return # first logic?
        
        if not read:
            print("TrajFile.create: creating %s" % s.fileName)
            pass
        
        s.unit = openFile(s.fileName,
                          mode="READ" if read else "WRITE",
                          isFormatted=s.isFormatted)
        return

    def write(s,iter,time):
        """Write the current coordinates specified by the freeAtoms argument
        of the constructor.

        The time argument should be the current time for iter>0.
        """
        if s.unit < 0: return

        if not s.constantTime and iter%s.saveInterval != 0: return

        sim=s.freeAtoms.simulation()
        curPos = sim.atomPosArr()
        from vec3 import Vec3
        xy=Vec3(0,0,0)
        if s.fixXY:
            from atomAction import centerOfMass
            cm = centerOfMass(s.freeAtoms)
            xy=Vec3(cm[0],cm[1],0)
            curPos -= xy
            pass

        #iterate until nextTime>time
        while True:
            nextTime = s.nextIndex * s.delta_t
            if s.constantTime and time<nextTime:
                break

            stepPos = curPos

            if s.constantTime and iter>0:
                a = (time-nextTime) / (time-s.prevTime)
                b = (nextTime-s.prevTime) / (time-s.prevTime)
                stepPos = a * s.savedPos + b * curPos
                pass
            
            if s.printCoords:
                from atomAction import PrintPos
                print("TrajFile.write: step:", iter)
                s.freeAtoms.apply( PrintPos() )
                pass

            from xplorDCD import write
            write(s.freeAtoms,stepPos,s.start,s.isFirst,s.delta_t,
                  s.saveInterval,s.unit,s.isFormatted)
            s.isFirst=False

            s.nextIndex += s.saveInterval
            if not s.constantTime:
                break
            pass

        s.savedPos = curPos
        s.prevTime = time

        return

    def read(s,istep=1,start=1,stop=999999,skip=1):
        """Read the next frame from the trajectory file into the simulation
        specified in the constructor, optionally skipping the specified number
        of frames. On return the time member is set to the current time."""

        from xplorDCD import read
        (s.isDone, s.isEOF, s.isError,
         s.start, s.stop, s.skip, s.istep,
         s.delta_t, s.header) = read(s.freeAtoms,s.isFirst,s.unit,
                                     start,stop,skip,
                                     istep,s.delta_t,s.header,
                                     s.isFormatted)

#        print "read: status: ",
#        print "isDone:", s.isDone
#        print "isEOF:", s.isEOF
#        print "isError:", s.isError
#        print "istep:", s.istep
#        print "header:", s.header
#        print "delta:", s.delta_t
#        print "start:", s.start
#        print "stop:", s.stop
#        print "skip:", s.skip

        if s.isFirst:
            s.curIndex=0
            s.isFirst=False
        else:
            s.curIndex += s.skip
            pass

        s.time = (s.start + s.curIndex) * s.delta_t
        
        return not s.isDone

    def readNext(s,skip=None):
        """Read the next trajectory frame. If the skip argument is specified,
        skip by that number of steps.
        """
        if skip:
            s.skip=skip
            pass
        import sys
        if s.isFirst:
            return s.read(1)
        return s.read(s.istep+s.skip,s.start,s.stop,s.skip)
        

    def close(s):
        if s.unit < 0: return
        closeFile(s.unit)
        s.isFirst=1
        s.unit=-1
        return
    pass
