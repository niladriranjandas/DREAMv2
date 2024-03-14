#!/usr/bin/env python
"""convert PIPP PCK files into XPLOR assignment statements
"""



diag = {}
cross = {}
cross1 = {}
cross2 = {}

C_RELAX = 1.25
DEPHASING_TIME = 0.0261

def convert(filename,
            c_relax        = C_RELAX,
            dephasing_time = DEPHASING_TIME):
    """ first argument: a list of lines from a pck file
    second argument: c_relax            [default value: 1.25      ]
    third argument: dephasing_time (2T) [default value: 0.0261 sec]
    output: assignment statements appropriate for xplor
    
     Input to this should be a .PCK file
     with Assignments in Column 7 
     and the .PCK and Intensity in Column 6 and ID in column 1
     Assignment MUST be of form I14.HA or I14.HN to that this
     script can differentiate cross and diagonal peaks
     There Can be anything between the . and the Atom name as long
     as the atom name is the last part of the assignment string
     The period MUST be immediately following the residue number
     script calculates Corresponding 3-bond coupling constant
     writes to stdout a Table containing:
		Residue   Diag Intens (Pk-ID)   X Peak Intens (Pk-ID)     Jhnha
    """
    class Peak:
        def __init__(s,id,index,hy,hn,n,intensity,label):
            s.id        = id
            s.index     = index
            s.hy        = hy
            s.hn        = hn
            s.n         = n
            s.intensity = intensity
            s.label     = label
            return
        pass

    def addPeak(hash,index,peak):
        #
        # perform simple heuristics 
        #
        coupledTo = peak.label.split(",")[-1]
        import sys
        err = sys.stdout.write
        if coupledTo.find("HA")>0 and (peak.hy<3 or peak.hy>6.5):
            err("addPeak: peak #%4d has an unusual chemical shift for an HA.\n" %
                peak.id)
            pass
        if index not in hash:
            hash[index] = peak
            return
        
        if hash[index].intensity < peak.intensity:
            hash[index] = peak
            pass
        return


    from time import time,ctime
    from math import atan,sqrt,pi
    import re, os
    ret = ""
    diag = {}
    cross = {}
    cross1 = {}
    cross2 = {}
    for line in open(filename).readlines():
        match = re.search(r"""^\s*
        ([0-9]+)\s+         # Peak ID
        ([0-9]+)\s+         # index(Z)
        ([0-9.+-eE]+)\s+     # Hy
        ([0-9.+-eE]+)\s+     # HN
        ([0-9.+-eE]+)\s+     # N
        ([0-9.+-eE]+)\s+     # Intensity
        ([^ \t\n]+)\s*         # Label
        """,line, re.VERBOSE)
        if match:
            mGroups = match.groups()
            id        = int(mGroups[0])
            index     = int(mGroups[1])
            hy        = float(mGroups[2])
            hn        = float(mGroups[3])
            n         = float(mGroups[4])
            intensity = abs( float(mGroups[5]) )
            label     = mGroups[6]
            if label[0] == '#' or label[0] == '*':
                continue
            resid     = int( label[1:].split(".")[0] )
            coupledTo = label.split(",")[-1]
            peak = Peak(id,index,hy,hn,n,intensity,label)
            if coupledTo == "HN":
                addPeak(diag,resid,peak)
            elif coupledTo[-2:] == "A1":
                addPeak(cross1,resid,peak)
                pass
            elif coupledTo[-2:] == "A2":
                addPeak(cross2,resid,peak)
                pass
            elif coupledTo.find("HA")>=0:
                addPeak(cross,resid,peak)
                pass
            pass
        pass

    #peaks.sort( (lambda p1,p2: cmp(p1.resid,p2.resid)) )



    ret += "! hnha_pck2xplr.nawk Table; Input File: %s\n" % filename
    ret += "!   Date: %s\n!   Path: %s\n" % (ctime(time()), os.getcwd() )
    ret += "!   c_relax        = %f\n" % c_relax
    ret += "!   dephasing_time = %f\n" % dephasing_time
    ret += "!\n! No GLY Table\n!\n"


    pi_2_t = pi * DEPHASING_TIME

    resids = list(diag.keys())
    resids.sort()

    print("found %5d diagonal peaks" % len(resids))

    nonGpeaks=0
    for id in resids:
        if diag[id].label[0] == "G": continue
        ret += "! %s%1d\n" % (diag[id].label[0], id)
        if id in list(cross.keys()):
            #not GLY
            nonGpeaks += 1
            jhnha = c_relax* \
                    atan( sqrt(cross[id].intensity/diag[id].intensity)) \
                    / pi_2_t
            ret += "assign (resid %3d and name c ) (resid %3d and name n )\n"\
                   % (id-1, id)
            ret += "       (resid %3d and name ca) (resid %3d and name c )"\
                   % (id,id) + \
                   "       %4.1f 0.5\n" % jhnha
            pass
        pass

    print("      %5d non-glycines contributing assignments peaks" % nonGpeaks)


    #            pass # end if Correct rs and RES



    ret +=  "\n\n"
    ret +=  "! hnha_pck2xplr.nawk Table; Input File: %s\n" % filename
    ret +=  "!   Date: %s\n!   Path: %s\n" % (ctime(time()), os.getcwd() )
    ret +=  "!   c_relax        = %f\n" % c_relax
    ret +=  "!   dephasing_time = %f\n" % dephasing_time
    ret +=  "!\n"
    ret +=  "! GLY only Table\n"
    ret +=  "!\n"

    gPeaks=0
    for id in resids:
        if diag[id].label[0] != "G": continue
        #not GLY
        ret +=  "! %s%1d\n" % (diag[id].label[0], id)
        found=0
        if id in list(cross1.keys()):
            found=1
            num="1"
            coupledTo = cross1[id].label.split(",")[-1]
            if coupledTo.find("HA2")>=0:
                #ambiguous
                num="#"
                pass
            jhnha1 = c_relax*\
                     atan( sqrt(cross1[id].intensity/diag[id].intensity)) /\
                     pi_2_t
            ret += "assign (resid %3d and name ha%s) " % (id,num) + \
                   "(resid %3d and name ca)\n" % id
            ret += "       (resid %3d and name n ) (resid %3d and name hn)"\
                   %   (id, id) + \
                   " %5.1f 0.5\n" % jhnha1
            pass
        if id in list(cross2.keys()):
            found=1
            num="2"
            coupledTo = cross2[id].label.split(",")[-1]
            if coupledTo.find("HA1")>=0:
                #ambiguous
                num="#"
                pass
            jhnha2 = c_relax*\
                     atan( sqrt(cross2[id].intensity/diag[id].intensity)) /\
                     pi_2_t

            ret += "assign (resid %3d and name ha%s) " % (id,num) + \
                   "(resid %3d and name ca)\n" % id
            ret += "       (resid %3d and name n ) " % id+ \
                   "(resid %3d and name hn)" % id + \
                   " %5.1f 0.5\n" % jhnha2
            pass
        if found:
            gPeaks += 1
            pass
        pass
    print("      %5d glycines contributing assignments" % gPeaks)
    return ret


if __name__ == "__main__":
    import sys
    filename = sys.argv[1]
    print(convert(open(filename).readlines()))
    pass

