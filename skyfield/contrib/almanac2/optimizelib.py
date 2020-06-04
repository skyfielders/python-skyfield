"""
The variable ``f`` is the plain objective function that returns numbers 
directly from skyfield.

The variable ``g`` is ``f`` transformed such that either all the target 
values appear to the secant method as roots, or such that all of the 
extremes appear to Brent's method as minima.
"""

import numpy

def secant(f, jd0, jd1, targets=0, f0=None, f1=None, tol=1e-10):
    """ Performs secant search to find target values in f.
    
    This is a vectorized version of the function newton from pyephem:
    https://github.com/brandon-rhodes/pyephem/blob/f96daf12d4f815be92e0caa52611b444517b9e0d/ephem/__init__.py#L91
    
    ``f`` is the plain objective function, and ``g`` is the objective function 
    translated vertically such that all of the target values appear to the 
    algorithm as roots, and renormalized to -180 to 180 degrees so that the 
    partitions won't contain any discontinuities.
    
    Arguments
    ---------
    f : function
        Objective function. Must accept and return ndarrays.
    jd0 : ndarray
        Left sides of partitions known to contain the ``targets``
    jd1 : ndarray
        Right sides of partitions known to contain the ``target``
    targets : float or ndarray
        Target values corresponding to the partitions defined by jd0 and jd1
    f0 : ndarray
        ``f(jd0)``. Providing this saves an extra function call
    f1 : ndarray
        ``f(jd1)``. Providing this saves an extra function call
    tol : float
        Tolerance used to determine when convergence is complete.
        
    Returns
    -------
    jd : ndarray
        the jd values at which ``f(jd) == targets``
    """
    if numpy.isscalar(targets):
        targets = numpy.ones_like(jd0) * targets
    
    def g(t, targets):
        return (f(t) - targets + 180)%360 - 180
    
    g0 = (f0 - targets + 180)%360 - 180 if f0 is not None else g(jd0, targets)
    g1 = (f1 - targets + 180)%360 - 180 if f1 is not None else g(jd1, targets)
    
    iters = 0
    max_iters = 50
    while iters < max_iters:   
        converged = (g1 == 0) + (abs(jd1 - jd0) < tol) + (g1 == g0)
        if converged.all():
            break
        inds = ~converged

        jd0[inds], jd1[inds] = jd1[inds], jd1[inds] + (jd1[inds] - jd0[inds]) / (g0[inds]/g1[inds] - 1)
        g0[inds], g1[inds] = g1[inds], g(jd1[inds], targets[inds])
        iters += 1
    return jd1


def _bracket(func, xa, xb, multiplier=None, f0=None, f1=None):
    """ This function is only meant to be used by brent_min
    """
    if multiplier is None:
        multiplier = numpy.ones_like(xa)
    
    _gold = 1.618034  # golden ratio: (1.0+sqrt(5.0))/2.0
    _verysmall_num = 1e-21
    grow_limit=110.0
    maxiter=1000
    failing = numpy.ones_like(xa, dtype=bool)
    
    def g(x, multiplier):
        return multiplier*func(x) + (multiplier==1)*360
    
    fa = g(xa, multiplier) if f0 is not None else multiplier*f0 + (multiplier==1)*360
    fb = g(xb, multiplier) if f1 is not None else multiplier*f1 + (multiplier==1)*360
    
    # Switch so fa > fb
    ind1 = fa<fb
    xa[ind1], xb[ind1] = xb[ind1], xa[ind1]
    fa[ind1], fb[ind1] = fb[ind1], fa[ind1]
    
    xc = xb + _gold * (xb - xa)
    fc = g(xc, multiplier)
    iter_ = 0
    
    failing[fc>=fb] = False
    
    while failing.any():
        tmp1 = (xb - xa) * (fb - fc)
        tmp2 = (xb - xc) * (fb - fa)
        val = tmp2 - tmp1
        
        denom = numpy.zeros_like(val)
        ind2 = numpy.abs(val) < _verysmall_num
        denom[ind2] = 2.0 * _verysmall_num
        denom[~ind2] = 2.0 * val
        
        w = xb - ((xb - xc) * tmp2 - (xb - xa) * tmp1) / denom
        wlim = xb + grow_limit * (xc - xb)
        if iter_ > maxiter:
            raise RuntimeError("Too many iterations.")
        iter_ += 1
        
        case1 = (w - xc) * (xb - w) > 0.0
        case2 = ~case1 * ((w - wlim)*(wlim - xc) >= 0.0)
        case3 = (~case1 + ~case2) * ((w - wlim)*(xc - w) > 0.0)
        case4 = ~case1 * ~case2 * ~case3
        
        fw = numpy.zeros_like(xa)
        
        fw[case1] = g(w[case1], multiplier[case1])
        
        ind1 = case1 * (fw < fc) * failing
        ind2 = case1 * ~ind1 * (fw > fb) * failing

        xa[ind1] = xb[ind1]
        xb[ind1] = w[ind1]
        fa[ind1] = fb[ind1]
        fb[ind1] = fw[ind1]
        failing[ind1] = False

        xc[ind2] = w[ind2]
        fc[ind2] = fw[ind2]
        failing[ind2] = False
            
        w[case1] = xc[case1] + _gold * (xc[case1] - xb[case1])
        fw[case1] = g(w[case1], multiplier[case1])

        w[case2] = wlim[case2]
        if case2.any():
            fw[case2] = g(w[case2], multiplier[case2])

        if case3.any():
            fw[case3] = g(w[case3], multiplier[case3])
        
        ind3 = case3 * (fw < fc) * failing
        xb[ind3] = xc[ind3]
        xc[ind3] = w[ind3]
        w[ind3] = xc[ind3] + _gold * (xc[ind3] - xb[ind3])
        fb[ind3] = fc[ind3]
        fc[ind3] = fw[ind3]
        
        if ind3.any():
            fw[ind3] = g(w[ind3], multiplier[ind3])

        w[case4] = xc[case4] + _gold * (xc[case4] - xb[case4])
        
        if case4.any():
            fw[case4] = g(w[case4], multiplier[case4])
        
        xa[failing] = xb[failing]
        xb[failing] = xc[failing]
        xc[failing] = w[failing]
        fa[failing] = fb[failing]
        fb[failing] = fc[failing]
        fc[failing] = fw[failing]
    return xa, xb, xc, fa, fb, fc


def brent_min(f, jd0, jd1, minimum=True, f0=None, f1=None, tol=1.48e-8):
    """ Vectorized version of Brent's method for finding maxima and minima.
    
    This is a vectorized version of scipy.optimize.optimize.Brent.optimize:
    https://github.com/scipy/scipy/blob/de9e0d6022be0319c1ba73c07a0946be46e02a24/scipy/optimize/optimize.py#L1896
    
    ``f`` is the plain objective function, and ``g`` is the objective function 
    optionally mirrored over the horizontal axis such that all of the extrema 
    appear to the algorithm as minima.
    
    Arguments
    ---------
    f : function
        Objective function. Must accept and return ndarrays.
    jd0 : ndarray
        Left sides of partitions known to contain maxima or minima
    jd1 : ndarray
        Right sides of partitions known to contain the target value
    minimum : bool or ndarray
        Whether or not the corresponding partitions contain minima. True means 
        the corresponding partition contains a minimum, False means it 
        contains a maximum.
    f0 : ndarray
        f(jd0). Providing this saves an extra function call
    f1 : ndarray
        f(jd1). Providing this saves an extra function call
    tol : float
        Tolerance used to determine when convergence is complete.
        
    Returns
    -------
    jd : ndarray
        the jd values at which f(jd) is a maximum or minimum.
    """
    if isinstance(minimum, bool):
        if minimum:
            multiplier = numpy.ones_like(jd0)
        else:
            multiplier = -numpy.ones_like(jd0)
    else:
        multiplier = numpy.ones_like(minimum, dtype=int)
        multiplier[~minimum] = -1
    
    maxiter = 500
    mintol = 1e-11
    
    def g(x, multiplier):
        return multiplier*f(x) + (multiplier==1)*360
    
    xa, xb, xc, fa, fb, fc = _bracket(f, jd0, jd1, multiplier, f0, f1)
    _cg = 0.3819660
    x = numpy.copy(xb)
    w = numpy.copy(xb)
    v = numpy.copy(xb)
    fw = numpy.copy(fb)
    fv = numpy.copy(fb)
    fx = numpy.copy(fb)
    
    a = numpy.zeros_like(xa)
    b = numpy.zeros_like(xb)
    case1 = xa < xc
    case2 = ~ case1
    a[case1] = xa[case1]
    b[case1] = xc[case1]
    a[case2] = xc[case2]
    b[case2] = xa[case2]

    deltax =  numpy.zeros_like(a)
    u =  numpy.zeros_like(a)
    fu = numpy.zeros_like(a)
    rat = numpy.zeros_like(a)
    tmp1 = numpy.zeros_like(a)
    tmp2 = numpy.zeros_like(a)
    p = numpy.zeros_like(a)
    dx_temp = numpy.zeros_like(a)

    iter_ = 0
    while (iter_ < maxiter):
        tol1 = tol * numpy.abs(x) + mintol
        tol2 = 2.0 * tol1
        xmid = 0.5 * (a + b)
        
        converged = numpy.abs(x - xmid) < (tol2 - 0.5 * (b - a))
        if converged.all():
            break
        not_converged = numpy.nonzero(~converged)[0]

        for i in not_converged:
            if (numpy.abs(deltax[i]) <= tol1[i]):
                if (x[i] >= xmid[i]):
                    deltax[i] = a[i] - x[i]       # do a golden section step
                else:
                    deltax[i] = b[i] - x[i]
                rat[i] = _cg * deltax[i]
            else:                              # do a parabolic step
                tmp1[i] = (x[i] - w[i]) * (fx[i] - fv[i])
                tmp2[i] = (x[i] - v[i]) * (fx[i] - fw[i])
                p[i] = (x[i] - v[i]) * tmp2[i] - (x[i] - w[i]) * tmp1[i]
                tmp2[i] = 2.0 * (tmp2[i] - tmp1[i])
                if (tmp2[i] > 0.0):
                    p[i] = -p[i]
                tmp2[i] = numpy.abs(tmp2[i])
                dx_temp[i] = deltax[i]
                deltax[i] = rat[i]
                # check parabolic fit
                if ((p[i] > tmp2[i] * (a[i] - x[i])) and (p[i] < tmp2[i] * (b[i] - x[i])) and
                        (numpy.abs(p[i]) < numpy.abs(0.5 * tmp2[i] * dx_temp[i]))):
                    rat[i] = p[i] * 1.0 / tmp2[i]        # if parabolic step is useful.
                    u[i] = x[i] + rat[i]
                    if ((u[i] - a[i]) < tol2[i] or (b[i] - u[i]) < tol2[i]):
                        if xmid[i] - x[i] >= 0:
                            rat[i] = tol1[i]
                        else:
                            rat[i] = -tol1[i]
                else:
                    if (x[i] >= xmid[i]):
                        deltax[i] = a[i] - x[i]  # if it's not do a golden section step
                    else:
                        deltax[i] = b[i] - x[i]
                    rat[i] = _cg * deltax[i]
    
            if (numpy.abs(rat[i]) < tol1[i]):            # update by at least tol1
                if rat[i] >= 0:
                    u[i] = x[i] + tol1[i]
                else:
                    u[i] = x[i] - tol1[i]
            else:
                u[i] = x[i] + rat[i]
                
        fu[not_converged] = g(u[not_converged], multiplier[not_converged])      # calculate new output value

        for i in not_converged:
            if (fu[i] > fx[i]):                 # if it's bigger than current
                if (u[i] < x[i]):
                    a[i] = u[i]
                else:
                    b[i] = u[i]
                if (fu[i] <= fw[i]) or (w[i] == x[i]):
                    v[i] = w[i]
                    w[i] = u[i]
                    fv[i] = fw[i]
                    fw[i] = fu[i]
                elif (fu[i] <= fv[i]) or (v[i] == x[i]) or (v[i] == w[i]):
                    v[i] = u[i]
                    fv[i] = fu[i]
            else:
                if (u[i] >= x[i]):
                    a[i] = x[i]
                else:
                    b[i] = x[i]
                v[i] = w[i]
                w[i] = x[i]
                x[i] = u[i]
                fv[i] = fw[i]
                fw[i] = fx[i]
                fx[i] = fu[i]

        iter_ += 1
    return x
