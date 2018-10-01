import numpy

def newton(f, x0, x1, tol=1e-10):
    max_iters = 50
    f0, f1 = f(x0), f(x1)
    
    iters = 0
    while iters < max_iters:   
        converged = (f1 == 0) + (abs(x1 - x0) < tol) + (f1 == f0)
        if converged.all():
            break
        inds = ~converged

        x0[inds], x1[inds] = x1[inds], x1[inds] + (x1[inds] - x0[inds]) / (f0[inds]/f1[inds] - 1)
        f0[inds], f1[inds] = f1[inds], f(x1[inds])
        iters += 1
    return x1


def bracket(func, xa, xb, grow_limit=110.0, maxiter=1000):
    _gold = 1.618034  # golden ratio: (1.0+sqrt(5.0))/2.0
    _verysmall_num = 1e-21
    failing = numpy.ones_like(xa, dtype=bool)
    
    fa = func(xa)
    fb = func(xb)
    
    # Switch so fa > fb
    ind = fa<fb
    xa[ind], xb[ind] = xb[ind], xa[ind]
    fa[ind], fb[ind] = fb[ind], fa[ind]
    
    xc = xb + _gold * (xb - xa)
    fc = func(xc)
    funcalls = 3
    iter_ = 0
    
    failing[fc>=fb] = False
    
    while failing.any():
        tmp1 = (xb - xa) * (fb - fc)
        tmp2 = (xb - xc) * (fb - fa)
        val = tmp2 - tmp1
        
        denom = numpy.zeros_like(val)
        ind = numpy.abs(val) < _verysmall_num
        denom[ind] = 2.0 * _verysmall_num
        denom[~ind] = 2.0 * val
        
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
        
        fw[case1] = func(w[case1])
        funcalls += 1
        
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
        fw[case1] = func(w[case1])
        funcalls += 1

        w[case2] = wlim[case2]
        if case2.any():
            fw[case2] = func(w[case2])
            funcalls += 1

        if case3.any():
            fw[case3] = func(w[case3])
            funcalls += 1
        
        ind = case3 * (fw < fc) * failing
        xb[ind] = xc[ind]
        xc[ind] = w[ind]
        w[ind] = xc[ind] + _gold * (xc[ind] - xb[ind])
        fb[ind] = fc[ind]
        fc[ind] = fw[ind]
        
        if ind.any():
            fw[ind] = func(w[ind])
            funcalls += 1

        w[case4] = xc[case4] + _gold * (xc[case4] - xb[case4])
        
        if case4.any():
            fw[case4] = func(w[case4])
            funcalls += 1
        
        xa[failing] = xb[failing]
        xb[failing] = xc[failing]
        xc[failing] = w[failing]
        fa[failing] = fb[failing]
        fb[failing] = fc[failing]
        fc[failing] = fw[failing]
    return xa, xb, xc, fa, fb, fc, funcalls


def brent_min(f, left, right, tol=1.48e-8):
    # set up for optimization
    maxiter = 500
    mintol = 1e-11
    
    func = f
    xa, xb, xc, fa, fb, fc, funcalls = bracket(f, left, right)
    _cg = 0.3819660
    x = numpy.copy(xb)
    w = numpy.copy(xb)
    v = numpy.copy(xb)
    fw = func(x)
    fv = numpy.copy(fw)
    fx = numpy.copy(fw)
    
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

    funcalls = 1
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
                
        fu[not_converged] = func(u[not_converged])      # calculate new output value
        funcalls += 1

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