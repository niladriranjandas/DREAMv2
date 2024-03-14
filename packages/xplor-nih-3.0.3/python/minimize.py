"""
optimization routines.
"""

from cdsVector import CDSVector_double as vector
from cdsVector import norm, dot
from functools import reduce

class F1Dim:
    def __init__(s,f,p,d):
        s.f = f
        s.p = vector(p)
        s.d = vector(d)
#        print 'F1Dim:', list(d),list(s.d)
        return
    def __call__(s,x):
        v = s.p + float(x) * s.d
#        print 'in __call__', x,list(v),list(s.p),list(s.d)
       
        return s.f(v)
    pass




def linemin(func,
            p,
            xi  ,
            stepsize,
            tolerance=1e-8):
    """
    determine the minimum of func starting at point p, in the direction xi.

    return (pf,pfmin,xif)
      minimum position, fmin, and the actual displacement vector used
    
    """
    ax = 0
    xx = stepsize
#    print 'linemin', list(p),list(xi)
    (ax,xx,bx,fa,fx,fb) = mnbrak(ax,xx,
                                 F1Dim(func,p,xi));
    # if (ax < 0.0) cout << "linemin: ax<0: " << ax << '\n';
    # if (xx < 0.0) cout << "linemin: xx<0: " << xx << '\n';
    # if (bx < 0.0) cout << "linemin: bx<0: " << bx << '\n';
    # if ( !(ax<xx&&xx<bx) )
    #   cout << "linemin: abscissas out of order: " 
    #	  << ax << ' ' << xx << ' ' << bx << '\n';
    # if ( !( fx<fa && fx<fb) )
    #   cout << "linemin: minimum not bracketed: "
    #	  << fa << ' ' << fx << ' ' << fb << '\n';
   
 
#    print 'linemin1'
    (stepsize,fret) = brent(ax,xx,bx,
                            F1Dim(func,p,xi))
#    print 'linemin2', stepsize
    #// cout << "linemin: fret: " << fret << '\n';
    xi *= stepsize;
    p += xi;
#    print 'linemin3', list(p),list(xi)
    return (p,fret,xi,stepsize)

GOLD   = 1.618034;
GLIMIT = 100.0;
TINY   = 1.0e-20;


def SIGN(a,b):
    if b>0:
        return a
    else:
        return -a
    pass


def mnbrak(ax,
           bx,
           func):
    """given function func, and given distinct initial points ax and bx, this
    routine searches in the downhill direction *defined by the function as
    evaluated at the initial points) and returns new points ax, bx, cx which
    bracket a minimum of the function. Also returned are the function values 
    at the three points fa, fb and fc."""

    
    fa=func(ax);
    fb=func(bx);
#    print 'mnbrak', ax,bx,fa,fb
    if fb > fa:
        (ax,bx) = (bx,ax)
        (fa,fb) = (fb,fa)
        pass
    cx = bx + GOLD*(bx-ax)
    # cout << "nbrack: initial step: " << cx << ' '
    fc = func(cx)

    while fb > fc:
        r = (bx-ax) * (fb-fc)
        q = (bx-cx) * (fb-fa)
        u = bx - ((bx-cx)*q - (bx-ax)*r) / \
            (2.0*SIGN(max(abs(q-r),TINY),q-r))
        ulim = bx + GLIMIT*(cx-bx);
        fu=0
        if (bx-u)*(u-cx) > 0.0:
            #//     cout << "nbrack: parabolic step: " << u << ' ';
            fu = func(u);
            if fu < fc:
                (ax,bx) = (bx,u)
                (fa,fb) = (fb,fu)
                return (ax,bx,cx,fa,fb,fc)
            elif fu > fb:
                cx=u
                fc=fu
                return (ax,bx,cx,fa,fb,fc)
            u = cx + GOLD*(cx-bx);
            #//     cout << "nbrack: min eval: " << u << ' ';
            fu = func(u);
        elif (cx-u)*(u-ulim) > 0.0:
            #//     cout << "nbrack: hunt1 eval: " << u << ' ';
            fu=func(u);
            if fu < fc:
                (bx,cx,u) = (cx,u,cx+GOLD*(cx-bx))
                #//       cout << "nbrack: hunt2 eval: " << u << ' ';
                (fb,fc,fu) = (fc,fu,func(u))
                pass
            pass
        elif (u-ulim)*(ulim-cx) >= 0.0:
            u = ulim
            #//     cout << "nbrack: hunt3 eval: " << u << ' ';
            fu = func(u)
        else:
            u  = cx + GOLD*(cx-bx)
            #//     cout << "nbrack: hunt4 eval: " << u << ' ';
            fu = func(u)
            pass
        (ax,bx,cx) = (bx,cx,u)
        (fa,fb,fc) = (fb,fc,fu)
        pass
    return (ax,bx,cx,fa,fb,fc)


CGOLD =  0.3819660;

def brent(ax,
          bx,
          cx,
          f,
          tol=1e-8,
          ITMAX=100,
          ZEPS=1e-10):
    """
    Find minimum bracketed by ax,bx,cx
    return (xmin,f(xmin))
    """
    

    if ax<cx: #//a and b bracket minimum w/
        a=float(ax)
    else:
        a=float(cx)
        pass
    if ax>cx: #//  a<b
        b=float(ax)
    else:
        b=float(cx)
        pass

    #print 'brent0',ax,bx,cx,f(ax),f(bx),f(cx)
    (x,w,v) = (float(bx),float(bx),float(bx))
    fw = f(x);
    fv = float(fw);
    fx = float(fw);

    e=0 #distance moved on the step before last
    d=0;
    for iter in range(1,ITMAX+1):
#        print 'brent1',iter,'a:',a,'b:',b,'v:',v,'w:',w,'x:',x
        xm=0.5*(a+b);                            #midpoint of a and b
        tol1 = tol*abs(x)+ZEPS;
        tol2= 2 * tol1;
        if abs(x-xm) <= (tol2-0.5*(b-a)):     #stop if points close together
            xmin=float(x)
            return (xmin,fx)
        if abs(e) > tol1:
            r=(x-w)*(fx-fv);
            q=(x-v)*(fx-fw);
#            print x,v,q,w,r
            p=(x-v)*q-(x-w)*r;
            #print 'brent2:',p,q,r,e
            q=2*(q-r);
            if q > 0.0: p = -p;
            q=abs(q);
            etemp=float(e);
            e=float(d);
            if (abs(p) >= abs(0.5*q*etemp) or
                p <= q*(a-x)               or 
                p >= q*(b-x)                 ):
#                print "brent: taking section step: ", p, q, fv,fw,fx
                if x>=xm:
                    e=a-x
                else:
                    e=b-x
                    pass
                d=CGOLD*e
            else:
#                print "brent: taking parabolic step: ", p, q, fv,fw,fx
                d=p/q;
                u=x+d;
                if u-a < tol2 or b-u < tol2:
                    d=SIGN(tol1,xm-x)
                    pass
                pass
            pass
        else:
#            print "brent: taking section step2: ", fv,fw,fx
            if x>=xm:
                e=a-x
            else:
                e=b-x
                pass
            d=CGOLD*e
            pass
        if abs(d) >= tol1:
            u = x+d
        else:
            u= x+SIGN(tol1,d)
            pass
        fu=f(u);
#        print 'brent',u,fu
        if fu <= fx:
            if u >= x:
                a=float(x)
            else:
                b=float(x)
                pass
            (v,w,x) = (w,x,u)
            (fv,fw,fx) = (fw,fx,fu)
        else:

            if u < x:
                a=float(u)
            else:
                b=float(u)
                pass
            if fu <= fw or w == x:
                (v,w,fv,fw) = (w,u,fw,fu)
            elif fu <= fv or v == x or v == w:
                (v,fv) = (u,fu)
                pass
            pass
        pass
    raise Exception("Too many iterations in BRENT")
    xmin=float(x);
    return (xmin,fx);


costTol = 1e-8
gradTol = 1e-4

def conmin(x,f,df,costTol=1e-10,gradTol=1e-4,maxIters=200,stepsize=1,
           verbose=0):
    """
    conjugate gradient minimization of function f with gradient df, starting at
    position x.

    stepsize specifies initial stepsize to be taken

    Return value is the tuple (x,f(x),numIters)

    FIX: this routine makes unnecessary function calls.
    """
    

    x = vector(x)

    fx = f(x)
    
    dfx = vector(df(x))
    
    dfx *= -1;
    g = vector(dfx)
    h = vector(dfx)
    
    costPrev = float(fx);
    
    #minimize/fit
    done=0
    iter=0
    while not done:
        cost = float(fx);
    
        (x,fx,xi,stepsize) = linemin(f,x,dfx,stepsize,costTol)
#        fx = f(x)
        iter += 1
     
        if iter>maxIters:
            print("conmin: WARNING: too many iterations taken")
            break


        f(x)  #frequently necessary
        dfx = vector(df(x))
    
        if verbose:
            print("conmin: fx: ", iter,fx, "gg ", norm(dfx)**2, end=' ')
            print("stepsize", stepsize)
            pass

        if norm(dfx) <= gradTol:
            if verbose:
                print("conmin: exit: gradient tolerance achieved.")
            break
        if abs(fx-costPrev) <= costTol:
            if verbose:
                print("conmin: exit: improvement in cost <= costTol.")
            break
        
#        stepsize=abs(stepsize)
#        if stepsize<0:
#            if verbose:
#                print "conmin: resetting..."
#                pass
#            g = vector(dfx)
#            h = vector(dfx)
#            stepsize=1
#            continue


       
        costPrev = float(fx)
    
        gg  = norm(g)**2
        dgg = dot( dfx+g , dfx );
        if gg == 0.0: break;  # unlikely...
        
        gam = dgg / gg;
        g = -dfx;
        h *= gam; h += g;
        dfx = vector(h);
        pass

    return (x, f(x),iter)

def outerProd(x,y):
    """
    return the outer product matrix result of multiplying the vectors x and y
    """
    from cdsMatrix import CDSMatrix_double as matrix
    ret = matrix(len(x),len(y))
    for i in range(len(x)):
        for j in range(len(y)):
            ret[i,j] = x[i] * y[j]
            pass
        pass
    return ret
                   
                 
import sys
def bfgs(x,f,df,costTol=1e-10,gradTol=1e-4,maxIters=200,verbose=0,
         stepsize0=1,resetInterval=None,
         writer=sys.stdout.write):
    """
    Perform Broyden-Fletcher-Goldfarb-Shanno variant of
    Davidson-Fletcher-Powell minimization.

    The writer argument can be used to specify an alternate destination
    for verbose output. By default it goes to sys.stdout.write

    stepsize0 is the initial stepsize as a fraction of the initial
    gradient length.
    
    Return value is the tuple (x,J(x),numIters)
    """

    from cdsMatrix import CDSMatrix_double as matrix

    n=len(x)
    if resetInterval==None:
        resetInterval=2*n
    x=vector(x)

    fx = f(x)
    fp=float(fx)
    g  = vector(df(x))
    normg = max(norm(g),1e-30)
    
    stepsize0 /= normg #distance to search along the gradient
    stepsize = stepsize0
    hess = matrix(n,n,0.)
    hess.setDiag(1.)

    #initial line direction
    xi = -g

    xp = vector(x)

    iters=0
    while 1:
        if verbose:
            writer( "bfgs:  " +str( iters)+"  fx:  "+str(fx) + " gg  " +\
                    str(norm(g)**2) )
            writer(" dx " + str( norm(x-xp)) + "\n" )
            pass
        if verbose>1:
            writer("   x = " + str(x) + '\n')
        iters += 1
        if iters>maxIters:
            fret = fx
            writer("bfgs: WARNING: too many iterations taken\n")
            break

        xp=vector(x)
        while 1:
            try:
                (x,fx,xi,stepsize) = linemin(f,x,xi,stepsize,costTol)
                break
            except OverflowError:
                stepsize/=10
                pass
            pass

        if 2*abs(fx-fp) <= costTol*(abs(fx)+abs(fp)+1e-10):
            fret=fx
            if verbose:
                writer( "bfgs: normal exit: costTol reached\n")
                writer( "bfgs: fx:  " + str( iters) + " " + str(fx) + " gg " +\
                        str(norm(g)**2) + '\n')
                pass
            break

        fp=float(fx) #save function value for exit determination
        og= vector(g)
        #calculate the new gradient
        fret=f(x) #this is often necessary to get correct grad value
        g = vector(df(x))

        if norm(g) <= gradTol:
            fret=fx
            if verbose:
                writer("bfgs: normal exit: gradient tolerance achieved.\n")
            break

#        stepsize=abs(stepsize)
        if stepsize<0 or iters%resetInterval==0:
            #something wrong- reset the Hessian
            if verbose:
                writer("bfgs: resetting... (stepsize=%e)\n"%stepsize)
            from cdsMatrix import svd
            svdResults = svd(hess,'S')
            hess.set(0.)
            hess.setDiag(1.)
            xi = -g
            stepsize=stepsize0
            continue


        dg = g - og      #gradient difference

        hdg = hess * dg



        fac = dot(dg,xi)   #dot products for denominators
        fae = dot(dg,hdg)
        

        try:
            fac = 1. / fac
        except ZeroDivisionError:
            fac = 0.
            pass
        
        try:
            fad = 1. / fae
        except ZeroDivisionError:
            fad = 0.
            pass

        dg = fac*xi - fad*hdg

        

        # BFGS update formula
        hess += fac*outerProd(xi,xi) - fad*outerProd(hdg,hdg) + \
                fae*outerProd(dg,dg)
#        hess += fac*outerProd(xi,xi) - fad*outerProd(hdg,hdg)

#        print 'hess:',hess


        

        xi = -hess * g   #new search direction
        pass

    if verbose:
        writer( "bfgs:  " +str( iters)+"  fx:  "+str(fx) + " gg  " +\
                str(norm(g)**2) )
        writer(" dx " + str( norm(x-xp)) + "\n" )
        pass
        
    return (x,fret,iters)

# the following routines are from scipy optimize.fmin
#
#Copyright (c) 2001, 2002 Enthought, Inc.
#All rights reserved.
#
#Copyright (c) 2003-2009 SciPy Developers.
#All rights reserved.
#
#Redistribution and use in source and binary forms, with or without
#modification, are permitted provided that the following conditions are met:
#
#  a. Redistributions of source code must retain the above copyright notice,
#     this list of conditions and the following disclaimer.
#  b. Redistributions in binary form must reproduce the above copyright
#     notice, this list of conditions and the following disclaimer in the
#     documentation and/or other materials provided with the distribution.
#  c. Neither the name of the Enthought nor the names of its contributors
#     may be used to endorse or promote products derived from this software
#     without specific prior written permission.
#
#
#THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
#AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
#IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
#ARE DISCLAIMED. IN NO EVENT SHALL THE REGENTS OR CONTRIBUTORS BE LIABLE FOR
#
#
def wrap_function(function, args):
    ncalls = [0]
    def function_wrapper(x):
        ncalls[0] += 1
        #print "wrap_function:", x, args
        return function(x, *args)
    return ncalls, function_wrapper

def simplex(func, x0, args=(), xtol=1e-4, ftol=1e-4, maxIters=None, maxfun=None,
         full_output=0, verbose=1, retall=0, callback=None):
    """Minimize a function using the downhill simplex algorithm.

    :Parameters:

      func : callable func(x,*args)
          The objective function to be minimized.
      x0 : ndarray
          Initial guess.
      args : tuple
          Extra arguments passed to func, i.e. ``f(x,*args)``.
      callback : callable
          Called after each iteration, as callback(xk), where xk is the
          current parameter vector.

    :Returns: (xopt, {fopt, iter, funcalls, warnflag})

      xopt : ndarray
          Parameter that minimizes function.
      fopt : float
          Value of function at minimum: ``fopt = func(xopt)``.
      iter : int
          Number of iterations performed.
      funcalls : int
          Number of function calls made.
      warnflag : int
          1 : Maximum number of function evaluations made.
          2 : Maximum number of iterations reached.
      allvecs : list
          Solution at each iteration.

    *Other Parameters*:

      xtol : float
          Relative error in xopt acceptable for convergence.
      ftol : number
          Relative error in func(xopt) acceptable for convergence.
      maxIters : int
          Maximum number of iterations to perform.
      maxfun : number
          Maximum number of function evaluations to make.
      full_output : bool
          Set to True if fval and warnflag outputs are desired.
      verbose : bool
          Set to True to print convergence messages.
      retall : bool
          Set to True to return list of solutions at each iteration.

    :Notes:

        Uses a Nelder-Mead simplex algorithm to find the minimum of
        function of one or more variables.

    """
    fcalls, func = wrap_function(func, args)
    x0 = vector(x0)
    N = len(x0)
    if maxIters is None:
        maxIters = N * 200
    if maxfun is None:
        maxfun = N * 200

    rho = 1; chi = 2; psi = 0.5; sigma = 0.5;
    one2np1 = list(range(1,N+1))

    sim=[(0,0)]*(N+1)

    sim[0] = (x0, func(x0))
    if retall:
        allvecs = [sim]
    nonzdelt = 0.05
    zdelt = 0.00025
    for k in range(0,N):
        y = vector(x0)
        if y[k] != 0:
            y[k] = (1+nonzdelt)*y[k]
        else:
            y[k] = zdelt

        f = func(y)
        sim[k+1] = (y,f)
        pass

    ## sort so sim[0,:] has the lowest function value
    sim.sort( key=lambda x: x[1] )
        
    iterations = 1

    while (fcalls[0] < maxfun and iterations < maxIters):
        if (max([norm(v[0] - sim[0][0]) for v in sim[1:]]) <= xtol \
            and max([abs(x[1]-sim[0][1]) for x in sim[1:]]) <= ftol):
            break

        xbar = vector(reduce(lambda x,y: x+y,
                             [v[0] for v in sim[:-1]])) / N
        xr = (1+rho)*xbar - rho*sim[-1][0]
        fxr = func(xr)
        doshrink = 0

        if fxr < sim[0][1]:
            xe = (1+rho*chi)*xbar - rho*chi*sim[-1][0]
            fxe = func(xe)

            if fxe < fxr:
                sim[-1] = (xe, fxe)
            else:
                sim[-1] = (xr, fxr)
        else: # fsim[0] <= fxr
            if fxr < sim[-2][1]:
                sim[-1] = (xr, fxr)
            else: # fxr >= fsim[-2]
                # Perform contraction
                if fxr < sim[-1][1]:
                    xc = (1+psi*rho)*xbar - psi*rho*sim[-1][0]
                    fxc = func(xc)

                    if fxc <= fxr:
                        sim[-1] = (xc, fxc)
                    else:
                        doshrink=1
                else:
                    # Perform an inside contraction
                    xcc = (1-psi)*xbar + psi*sim[-1][0]
                    fxcc = func(xcc)

                    if fxcc < sim[-1][1]:
                        sim[-1] = (xcc, fxcc)
                    else:
                        doshrink = 1

                if doshrink:
                    for j in one2np1:
                        simj = sim[0][0] + sigma*(sim[j][0] - sim[0][0])
                        fj   = func(simj)
                        sim[j]=(simj,fj)
                        pass
                    pass

        sim.sort( key=lambda x: x[1] )

        if callback is not None:
            callback(sim[0])
        iterations += 1
        if retall:
            allvecs.append(sim[0])

    (x,fval) = min(sim,key=lambda x:x[1])
    warnflag = 0

    if fcalls[0] >= maxfun:
        warnflag = 1
        if verbose:
            print("Warning: Maximum number of function evaluations has "\
                  "been exceeded.")
    elif iterations >= maxIters:
        warnflag = 2
        if verbose:
            print("Warning: Maximum number of iterations has been exceeded")
    else:
        if verbose:
            print("Optimization terminated successfully.")
            print("         Current function value: %f" % fval)
            print("         Iterations: %d" % iterations)
            print("         Function evaluations: %d" % fcalls[0])


    if full_output:
        retlist = x, fval, iterations, fcalls[0], warnflag
        if retall:
            retlist += (allvecs,)
    else:
        retlist = x
        if retall:
            retlist = (x, allvecs)

    return retlist
