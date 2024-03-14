#!/usr/bin/env python3

# scan swig-generated c wrapper files
# check that appropriate new/delete pairs are present

# this is somewhat fragile. Run with -v to make sure variable names
# are being detected.
#
#

import re
import sys

from re import findall

mode="python"
verbose=False
files=[]
argv=sys.argv
i=1
while i<len(argv):
    if argv[i] == '-tcl':
        mode="tcl"
    elif argv[i] == '-py':
        mode="python"
    elif argv[i] == '-v':
        verbose=True
    else:
        files.append( argv[i] )
        pass
    i += 1
    pass

if mode == "tcl":
    convertFunc = "fromTcl"
else:
    convertFunc = "fromPy"
    pass   

                        

exitStatus=0

def findNested(start,stop,startIndx,buf,initDepth=1):
    indx=startIndx
    depth=initDepth
    started=False
    if depth>0: started=True
    while indx<len(buf):
        if buf[indx] == stop: depth-= 1
        if buf[indx] == start:
            started=True
            depth+= 1
            pass
        #print indx, depth, buf[indx]
        if depth<1 and started: break
        indx += 1
        pass
    if indx>=len(buf): indx=startIndx
    return indx


for fileName in files:
    
    lines = open(fileName).readlines()

    #strip // comments
    cnt=0
    while cnt<len(lines):
        lines[cnt] = re.sub(r"//.*$","",lines[cnt])
        cnt+=1
        pass
    
    buf = ''.join(lines)

    #strip /* *. comments
    cont=1
    while cont:
        m1 = re.search(r"/\*",buf)
        m2 = re.search(r"\*/",buf)
        if m1 and m2:
            nls = buf.count("\n", m1.start(), m2.span()[-1])
            buf = buf[:m1.start()] +  "\n"*nls + buf[m2.span()[-1]:]
        else:
            cont=None
            pass
        pass
    
    #condense strings to ""
    buf = re.sub(r'""','',buf,0)
    buf = re.sub(r'"(\\"|[^"])+"','',buf,0)
     
    #print buf
    #sys.exit()

    #
    # now process function blocks
    # 
    indx=0
    cont=1
    subBuf = buf
    write = sys.stdout.write
    while cont:
        m = re.search(r"([0-9a-zA-Z_]+)\([^\)]*(\(self\))?[^\)]*\)[ \t\n]*\{",
                      subBuf)

        if m==None:
            cont=0
            continue
        
        name = m.group(1)
        #indx=m.span()[-1]
        indx=m.start(1)
        endx = findNested('{','}',indx,subBuf,0)

        if indx==endx:
            subBuf=subBuf[m.end(1):]
            continue

        block=subBuf[indx:endx]

        #print "name:", name, indx,endx,block

        #only check wrapping functions
        if name.startswith("_wrap_"):

            #print "found block with name: " + name, name.startswith("_wrap_")
            #print block
            
            #remove catch clauses: this script only checks for normal return
            #conditions
            while 1:
                m = re.search(r"\scatch\s[^{]+\{",block)
                if not m: break
                c_indx = m.span()[0]
                c_endx = findNested('{','}',c_indx,block,initDepth=0)
                block = block[:c_indx] + block[c_endx:]
                pass

            
            nList = []
            nmatch = findall(r"([0-9a-zA-Z_]+)\s*=\s*(\([^\)]+\)\s*)?new\s",
                             block)
            
            for n in nmatch: nList.append( n[0] )
                    
            #new in fromPyWithNew function call
            nList += findall(convertFunc + r"WithNew\(\s*([^, ]+)\s*,",
                             block)

            #call constructor added in .i file
            nmatch = findall(r"([0-9a-zA-Z_]+)\s*=\s*(\([^\)]+\)\s*)?" +
                             r"new_[0-9a-zA-z_]+__SWIG",
                             block)
            for n in nmatch: nList.append( n[0] )

            #instances of variable deletions
            dList = findall(r"\sdelete\s+(?:\[\]\s+)?([0-9a-zA-Z_]+)",block)

            # memory management assigned to python
            if mode=="python":
                dList += findall(r"SWIGPY(?:_Python)?"+
                                 r"_NewPointerObj\(\(void \*\)\s*\(?([^, )]+)\)"+
                                 r"\s*,[^,]+,\s*1\)", block)
                pass

            if mode=="tcl":
                dList += findall("SWIGTCL"+
                                 r"_NewPointerObj\(\s*?([^, )]+)", block)

            #print name,dList
            #dList += findall(r"SWIGTCL(_Tcl)?"+_NewPointerObj\(\(void\*\)\s*([^, ]+)",
            #                 block)
            dList += findall(r"SWIGTCL"+
                             r"_NewInstanceObj\((?:interp,)?"+
                             r"\(void\*\)\s*([^, ]+)"+
                             r"\s*,[^,]+,\s*1\)", block)

            tmp = dList
            dList = []
            for el in tmp:
                if not el in dList: dList.append(el)
                pass

            if verbose:
                print(fileName+": in " + name + ": ")
                print("  allocated:", nList)
                print("  deleted  :", dList)
                
            for n in nList:
                if not n in dList and \
                   not (name.startswith("_wrap_new") and n=='result'):
                    write(fileName + ": in " + name + ": ")
                    print(n , " is not deleted")
                    exitStatus=1
                    pass
                pass
            for d in dList:
                if d=="temp": continue # temp objs used in accessors
                if not d in nList and \
                   not name.startswith("_wrap_delete"):
                    write(fileName + ": in " + name + ": ")
                    print(d , " is not allocated")
                    exitStatus=1
                    pass
                #see if it's assigned a default value (not through new)
                assignStatements = findall(d +
                                           r"\s*=\s*\(?:\([^\)]+\)\s*\)?.*",
                                           block)
#                print assignStatements
                for statement in assignStatements:
                    #if d=="arg1": print statement
                    if statement.find('new ')<0 and statement.find('new_')<0 \
                       and not statement.rstrip(' ;').endswith('0'):
                        write(fileName + ": in " + name + ": ")
                        write(d + " is deleted, " +
                              "but assigned from non-heap value\n")
                        print(statement)
                        exitStatus=1
                        pass
                    pass
                pass
            pass
             
        subBuf=subBuf[endx:]
        pass
    pass
        
sys.exit(exitStatus)



#
#
