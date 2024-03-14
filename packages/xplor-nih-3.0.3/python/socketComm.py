"""
  Communicate via TCP sockets.

  This module defines a Comm object allowing communications between
  slave and master processes. Currently, the slaves cannot directly
  communicate.

  This module implements barrier waits, collection of data from slaves and
  sending messages to all slaves.

  Barrier is robust with respect to machine and processes crashes.
  
"""
import socket
import sys

retryInterval=1

# classes needed:
#  constructor: starts up aborts it timeout is reached
#    need to call functions already written
#  barrier
#    master: block waiting for barrier string from each slave, sends string.
#    slave: writes barrier string to master, waits for string.
#   need timeout for master, slave
#   need to catch closed connections
#   return value: master: those slaves which responded.
#  collect
#    slaves: write value to master
#    master: read value fro each slave
#
# low-level:
#    readln, writeln
#    readData, writeData - these first send size w/ writeln, then send
#                          data.

class Comm:
    def __init__(s,numProcs,procNum,host,port=9021,timeout=138240,
                 startupTimeout=50,allowStartupFailures=1,
                 startupDelay=None):
        """
        conns gets updated to reflect processes which have died.
        conns0 does not get updated.

        if allowStartupFailures is set, this is passed to startServer. If
        any failures occur on startup, the processes are renumbered, and
        numProcs is appropriately reset.

        If specified, the startupDelay argument specifies a time (in seconds)
        to wait before connecting to the server. For the server (procNum 0),
        this has no effect.
        """
        startupTimeout += numProcs
        s.numProcs = numProcs
        s.procNum = procNum
        s.barrierCnt=0
        if procNum==0:
            s.conns = startServer(numProcs,port,startupTimeout,timeout,
                                  allowStartupFailures)
            procs = list(s.conns.keys())
            if len(procs) != s.numProcs-1:
                procs.sort()
                newConns={}
                for i in range(len(procs)):
                    newConns[i+1] = s.conns[procs[i]]
                    pass
                s.conns = newConns
                pass
            s.conns0 = s.conns
            s.numProcs = len(procs)+1
            for i in range(1,len(procs)+1):
                s.conns[i].writeln("%d"%i)
                s.conns[i].writeln("%d"%s.numProcs)
                pass
            pass
        else:
            if startupDelay:
                import time
                time.sleep(startupDelay)
                pass
            s.conn = startClient(host,port,procNum,startupTimeout,timeout)
            newProcNum = -1
            data = s.conn.readln()
            try:
                newProcNum = int(data)
            except ValueError:
                print("could not convert %s to integer." % data)
                raise
            if newProcNum != s.procNum:
                s.procNum = newProcNum
                writeDebug("process number reassigned to: %d\n" % s.procNum)
                pass
            s.numProcs = -1
            s.numProcs = int(s.conn.readln())
            pass
        #wait til everyone starts
        #FIX: should use barrier's return value
        s.barrier(timeout=startupTimeout)
        return
    def info(s,procNum):
        if s.procNum>0:
            return
        return s.conns0[procNum]
    def procs(s):
        """
        Return the a sorted list of the current proc numbers. Only to be
        called by proc 0. As it is only consumed by proc 0, 0 is not included
        in the returned list.
        """
        return sorted( s.conns.keys() )
        
    def barrier(s,timeout=-1):
        """wait for all slaves. On master, return a list of connections which
        reached the barrier. If a connection is found to be dead (process or
        machine crashed), it is not included in the returned connections.
        """
        s.barrierCnt += 1
        import inspect
        #get the line number of the calling frame
        # this is a check that the barriers come from the same line of
        # source code
        cookie = inspect.stack()[1][2]
        if s.procNum>0:
            if timeout==None or timeout>=0:
                oldtimeout = s.conn.sock.gettimeout()
                s.conn.sock.settimeout(timeout)
                pass
            import time
            time0 = time.time()
            try:
                s.conn.writeln("barrier-out-%d-%d"%(s.procNum,cookie))
                input = s.conn.readln()
            except socket.timeout:
                print("timeout value:", s.conn.sock.gettimeout())
                print("elapsed time: " , time.time()-time0)
                raise
            if input!=b"barrier-in-%d-%d"%(s.procNum,cookie):
                raise Exception("barrier: unexpected read in barrier: " +
                                str(input))
            if timeout==None or timeout>=0:
                s.conn.sock.settimeout(oldtimeout)
            return []
        procs=[0]
        newconns={}
        keys = list(s.conns.keys())
        keys.sort(key=lambda x: int(x))
        for key in keys:
            conn = s.conns[key]
            if timeout==None or timeout>=0:
                oldtimeout = conn.sock.gettimeout()
                conn.sock.settimeout(timeout)
            try:
                input=conn.readln()
                if input!=b"barrier-out-%d-%d"%(key,cookie):
                    raise Exception("barrier: unexpected read in barrier " +
                                    "for process %d" % key)
                else:
                    conn.writeln("barrier-in-%d-%d"%(key,cookie))
                    procs.append(key)
                    newconns[key] = conn
                    pass
                pass
            except:
                print("barrier: error reading from process %d" % key)
            if timeout==None or timeout>=0:
                conn.sock.settimeout(oldtimeout)
            pass
        s.conns = newconns
        return procs
    def collect(s,msg):
        """ send message to procNum 0- returned in sorted fashion
        """
        s.barrier()
        if s.procNum>0:
            s.conn.writeData(msg)
            return []
        ret = [msg]
        keys= list(s.conns.keys())
        keys.sort(key=lambda x: int(x))
        for key in keys:
            conn = s.conns[key]
            try:
                ret.append( conn.readData() )
            except:
                print("error in socketComm.collect. key: ", key)
                raise
            pass
        return ret
    def distribute(s,msg):
        """send message from procNum 0 to all others. Returns this message
        """
        s.barrier()
        if s.procNum==0:
            for key in list(s.conns.keys()):
                s.conns[key].writeData(msg)
        else:
            msg = s.conn.readData()
        return msg
    def writeDataTo(s,proc,msg):
        """write python data to the specified process. Can write any
        pickle-able object.
        """
        if proc==0:
            if not hasattr(s,"conn"):
                raise Exception("write from proc 0 to proc 0 not supported")
            s.conn.writeData(msg)
        else:
            if not hasattr(s,"conns"):
                raise Exception("write to a sister node not supported")
            s.conns[proc].writeData(msg)
            pass
        return
    def readDataFrom(s,proc):
        """read python data from the specified process. The return value is
        a Python object. If the specified process does not exist, return
        None.
        """
        if proc not in s.conns:
            return None
        if proc==0:
            if not hasattr(s,"conn"):
                raise Exception("read from proc 0 to proc 0 not supported")
            return s.conn.readData()
        else:
            if not hasattr(s,"conns"):
                raise Exception("read from a sister node not supported")
            return s.conns[proc].readData()
            pass
        return
    pass
        
        
import pickle
class Connection:
    def __init__(s,sock):
        s.sock   = sock
        s.input  = sock.makefile("rb")
        s.output = sock.makefile("wb",0)
        # detect crashed hosts - within two hours
        s.sock.setsockopt(socket.SOL_SOCKET,socket.SO_KEEPALIVE,1)
        try:
            # if the following works, connections will be closed 6.5 min
            # after a machine crashes
            #
            #idle time before keepalive packets are sent: 5min
            s.sock.setsockopt(socket.SOL_TCP,socket.TCP_KEEPIDLE,300)
            #packets to send before closing a connection: 3
            s.sock.setsockopt(socket.SOL_TCP,socket.TCP_KEEPCNT,3)
            #how often to send packets: 30 sec
            s.sock.setsockopt(socket.SOL_TCP,socket.TCP_KEEPINTVL,30)
        except AttributeError:
            # not support with current Mac OS.
            #writeDebug('Notice: Setting tcp keepalive parameters not supported.')
            pass
        peername=sock.getpeername()[0]
        s.remoteHost=peername
        try:
            s.remoteHost = socket.gethostbyaddr(peername)[0]
        except socket.herror:
            print("socketCom.Connection: warning: ", end=' ')
            print("could not obtain a name for", peername)
            pass
        pass
    def readln(s):
        """ does not return the trailing newline character.
        """
        r=s.input.readline()[:-1]
        return r
    def writeln(s,str):
        str += '\n'
        str = bytes(str,'utf8')
        lenToSend = len(str)
        lenSent=0
        #note: the file interface (s.output) is not safe against interrupted
        #syscalls - hence the use of the send function
        while lenSent < lenToSend:
            sent = s.sock.send( str[lenSent:] )
            if sent == 0:
                raise RuntimeError("socket connection broken")
            lenSent += sent
            pass
        return
    def writeData(s,msg):
        """ msg can be any picklable Python object
        """
        pickled = pickle.dumps(msg)
        s.writeln("%d"%len(pickled))
        toSend=len(pickled)
        sent=0
        while sent<toSend:
            sent += s.output.write(pickled[sent:])
            pass
        s.output.flush()
    def readData(s):
        len = int(s.readln())
        pickled = s.input.read(len)
        return pickle.loads(pickled)
    def closed(s):
        return s.input.closed
    pass


def writeDebug(msg):
    #FIX: should olnly print if logLevel!=none
    print(msg)
    return

def startServer(numConnect,port,startupTimeout,timeout,
                allowStartupFailures=0):
    """
    start server and wait for numConnect connections.
    Return dictionary of connections whose keys correspond to process numbers.

    if allowIncomplete is set, reaching startupTimeout does not throw an
    exception: just that process is skipped.
    """
    if numConnect<2:
        return {}
    try:
        server = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        server.setsockopt(socket.SOL_SOCKET,socket.SO_REUSEADDR,1)
        server.bind(('',port))
        server.listen(5)
    except socket.error as why:
        import xplor
        xplor.writeConsole("can't setup server on port %d: %s." % (port,why) 
                           + " Try a different port.\n")
        sys.exit("can't setup server on port %d: %s." % (port,why) 
                 + " Try a different port.")
        pass
    writeDebug("[Server accepting clients]")
    
    server.settimeout(startupTimeout)
      
    proc={}
    import os
    
    nodeFilename=os.path.join(os.environ['NODE_FILENAME'])
    nodeFile = open(nodeFilename,"w")
    nodeFile.write("%s %s\n" % (os.environ['XPLOR_PHOST'],
                                os.getpid()))
    for i in range(numConnect-1):
        try:
            (sock, rem_addr) = server.accept()
        except socket.timeout:
            writeDebug("startServer: timeout waiting for connection")
            server.close()
            for p in range(1,numConnect):
                if not p in list(proc.keys()):
                    print("\tprocess not heard from:", p)
                    pass
                pass

            if allowStartupFailures:
                break
            else:
                raise("startServer: timeout waiting for connection")
            


#zero sndbuf not supported under Darwin, Irix
#        sock.setsockopt(socket.SOL_SOCKET,socket.SO_SNDBUF,0)

        sock.settimeout(timeout)

        conn = Connection(sock)
        
        procNum = int( conn.readln() )
        processID = int( conn.readln() )
        conn.writeln("%d" % procNum)
        writeDebug("[Connect from %s (process %d)]\n" %
                   (conn.remoteHost,procNum) )
        nodeFile.write("%s %s\n" % (conn.remoteHost,processID))
        
        proc[procNum] = conn
        
        pass

    # the close command stops this server from listening on the port
    #
    server.close()
    
    return proc

def startClient(host,port,procNum,startupTimeout,timeout):
    """connect to port on host
    """
    import time
    for i in range( int(startupTimeout) ):
        writeDebug("connect attempt #%d to %s:%d..." %
                   (i,host,port))
        try:
            sock = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
            #attempt to set close on exec property - so
            # children don't inherent the file descriptor
            try:
                import fcntl
                fcntl.fcntl(sock,fcntl.F_SETFD,fcntl.FD_CLOEXEC)
            except:
                writeDebug('unable to set close-on-exec')
                pass
            sock.settimeout(startupTimeout)
            sock.connect( (host,port) )
            
            sock.settimeout(timeout)
            conn = Connection(sock)
            conn.writeln("%d" % procNum)
            import os
            conn.writeln("%d" % os.getpid())
            if procNum != int(conn.readln()):
                raise Exception("connection fault: communication garbled")
            break

        except socket.error as why:
            writeDebug("connection fault: " + str(why))
            time.sleep(retryInterval)
            sock=0
            pass
        pass
    if not sock: raise Exception("%s:%d NOT CURRENTLY REACHABLE" %
                                 (host,port))

    return conn

        
