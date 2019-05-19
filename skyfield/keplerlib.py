from skyfield.functions import length_of, dots
import sys
from numpy import (sum, array, arange, power, cross, sqrt, log, exp, cos, sin, 
                   cosh, sinh, zeros_like, abs, tile, repeat, newaxis, amax, 
                   amin, ones_like)
import math

dpmax = sys.float_info.max


def bracket(num, end1, end2):
    num[num<end1] = end1
    num[num>end2] = end2
    return num
    
    
def find_trunc():
    denom = 2
    factr = 2
    trunc = 1
    x = 1 / denom
    while 1+x > 1:
        denom = denom * (2+factr) * (1+factr)
        factr = factr + 2
        trunc = trunc + 1
        x = 1 / denom
    return trunc


trunc = find_trunc()
odd_factorials = array([math.factorial(i) for i in range(3, trunc*2, 2)])
even_factorials = array([math.factorial(i) for i in range(2, trunc*2, 2)])


def stumpff(x):
    if (x < (-(log(2) + log(dpmax))**2)).any():
        raise ValueError('Argument below lower bound')
        
    z = sqrt(abs(x))
    
    c0 = zeros_like(x)
    c1 = zeros_like(x)
    c2 = zeros_like(x)
    c3 = zeros_like(x)
    
    low = x < -1
    c0[low] = cosh(z[low])
    c1[low] = sinh(z[low])/z[low]
    c2[low] = (1 - c0[low])/x[low]
    c3[low] = (1 - c1[low])/x[low]
    
    high = x > 1
    c0[high] = cos(z[high])
    c1[high] = sin(z[high])/z[high]
    c2[high] = (1 - c0[high])/x[high]
    c3[high] = (1 - c1[high])/x[high]
    
    mid = ~low * ~high
    n = sum(mid)
    exponents = tile(arange(0, trunc-1), [n, 1])
    odd_denominators = tile(odd_factorials, [n, 1])
    even_denominators = tile(even_factorials, [n, 1])
    numerators = repeat(x[mid][newaxis].T, trunc-1, axis=1)
    c3[mid] = sum(power(numerators, exponents)/odd_denominators, axis=1)
    c2[mid] = sum(power(numerators, exponents)/even_denominators, axis=1)
    c1[mid] = 1 - x[mid]*c3[mid]
    c0[mid] = 1 - x[mid]*c2[mid]
    
    return c0, c1, c2, c3
    

def propagate(position, velocity, t0, t1, gm):
    if gm <= 0: raise ValueError("'gm' should be positive")
    if length_of(velocity) == 0: raise ValueError('Velocity vector has zero magnitude')
    if length_of(position) == 0: raise ValueError('Position vector has zero magnitude')
        
    r0 = length_of(position)
    rv = dots(position, velocity)
    
    hvec = cross(position, velocity)
    h2 = dots(hvec, hvec)
    
    if h2 == 0: raise ValueError('Motion is not conical')
 
    eqvec = cross(velocity, hvec)/gm + -position/r0
    e = length_of(eqvec)
    q = h2 / (gm * (1+e))
    
    f = 1 - e
    b = sqrt(q/gm)
    
    br0 = b * r0
    b2rv = b**2 * rv
    bq = b * q
    qovr0 = q / r0
    
    
    maxc = max(abs(br0), 
               abs(b2rv),
               abs(bq),
               abs(qovr0/(bq)))
    
    if f < 0:
        fixed = log(dpmax/2) - log(maxc)
        rootf = sqrt(-f)
        logf = log(-f)
        bound = min(fixed/rootf, (fixed + 1.5*logf)/rootf)
    else:
        logbound = (log(1.5) + log(dpmax) - log(maxc)) / 3
        bound = exp(logbound)
        
    def kepler(x):
        c0, c1, c2, c3 = stumpff(f*x*x)
        return x*(br0*c1 + x*(b2rv*c2 + x*(bq*c3)))
    
    dt = t1 - t0
    
    x = bracket(dt/bq, -bound, bound)
    kfun = kepler(x)
    
    past = dt < 0
    future = dt > 0
    
    upper = zeros_like(dt, dtype='float64')
    lower = zeros_like(dt, dtype='float64')
    oldx = zeros_like(dt, dtype='float64')
    
    lower[past] = x[past]
    upper[future] = x[future]
    
    while (kfun[past] > dt[past]).any():
        upper[past] = lower[past]
        lower[past] *= 2
        oldx[past] = x[past]
        x[past] = bracket(lower[past], -bound, bound)
        if (x[past] == oldx[past]).any():
            raise ValueError('The input delta time (dt) has a value of {0}.'
                             'This is beyond the range of DT for which we '
                             'can reliably propagate states. The limits for '
                             'this GM and initial state are from {1}'
                             'to {2}.'.format(dt, kepler(-bound), kepler(bound)))
        kfun[past] = kepler(x[past])
        
    while (kfun[future] < dt[future]).any():
        lower[future] = upper[future]
        upper[future] *= 2
        oldx[future] = x[future]
        x[future] = bracket(upper[future], -bound, bound)
        if (x[future] == oldx[future]).any():
            raise ValueError('The input delta time (dt) has a value of {0}.'
                             'This is beyond the range of DT for which we '
                             'can reliably propagate states. The limits for '
                             'this GM and initial state are from {1} '
                             'to {2}.'.format(dt, kepler(-bound), kepler(bound)))    
        kfun[future] = kepler(x[future])
    
    x = amin(array([upper, amax(array([lower, (lower+upper)/2]), axis=0)]), axis=0)
    
    lcount = zeros_like(dt)
    mostc = ones_like(dt)*1000
    not_done = (lower < x) * (x < upper)
    
    while not_done.any():
        kfun[not_done] = kepler(x[not_done])
        
        high = (kfun > dt) * not_done
        low = (kfun < dt) * not_done
        same = (~high * ~low) * not_done
        
        upper[high] = x[high]
        lower[low] = x[low]
        upper[same] = lower[same] = x[same]        
        
        condition = not_done * (mostc > 64) * (upper != 0) * (lower != 0)
        mostc[condition] = 64
        lcount[condition] = 0
        
        # vectorized version of min(upper, max(lower, (upper + lower)/2))        
        x[not_done] = amin(array([upper[not_done], amax(array([lower[not_done], (lower[not_done]+upper[not_done])/2]), axis=0)]), axis=0)
        lcount += 1
        not_done = (lower < x) * (x < upper) * (lcount < mostc)
    
    c0, c1, c2, c3 = stumpff(f*x*x)
    br = br0*c0 + x*(b2rv*c1 + x*(bq*c2))
    
    pc = 1 - qovr0 * x**2 * c2
    vc = dt - bq * x**3 * c3
    pcdot = -qovr0 / br * x * c1
    vcdot = 1 - bq / br * x**2 * c2
    
    position_prop = pc*tile(position[newaxis].T, dt.size) + vc*tile(velocity[newaxis].T, dt.size)
    velocity_prop = pcdot*tile(position[newaxis].T, dt.size) + vcdot*tile(velocity[newaxis].T, dt.size)
    
    return position_prop, velocity_prop
