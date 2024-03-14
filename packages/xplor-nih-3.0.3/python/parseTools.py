
"""
Tools to help with string parsing.
"""

def readFloat(string):
    """
    Read a floating point number from the leading portion of string.
    Return a tuple of the number and the string with the number removed.
    """
    import re
    string=string.lstrip()
    m=re.match("[\-0-9.+efg]*",string,re.IGNORECASE)
    if not m:
        raise ValueError("Could not read leading floating point number: " +
                        string)
    num=float(m.group())
    return (num,string[len(m.group()):])

    
        
def readInt(string):
    """
    Read an integer from the leading portion of string.
    Return a tuple of the number and the string with the number removed.
    """
    import re
    string=string.lstrip()
    m=re.match("[\-0-9]*",string,re.IGNORECASE)
    if not m:
        raise ValueError("Could not read leading floating point number: " +
                        string)
    num=int(m.group())
    return (num,string[len(m.group()):])

    
        
def findNested(start,stop,startIndx,buf,initDepth=1):
    """
    Find nested single character delimiter pairs (start,stop) such as (), [],
    etc, in buf argument, starting at offset startIndx. initDepth specifies
    the number of leading start characters already seen. Returns the index of the
    closing delimiter.
    
    """
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

