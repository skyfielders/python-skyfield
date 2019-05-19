from skyfield.functions import length_of, dots
import sys
from numpy import sum, array, arange, power, cross, sqrt, log, exp, cos, sin, cosh, sinh
import math

dpmax = sys.float_info.max

def bracket(num, end1, end2):
    if end1 <= num <= end2:
        return num
    elif num < end1:
        return end1
    else:
        return end2
    
    
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
nterms = trunc
odd_factorials = array([math.factorial(i) for i in range(3, nterms*2, 2)])
even_factorials = array([math.factorial(i) for i in range(2, nterms*2, 2)])


def stumpff(x):    
    if x < (-(log(2) + log(dpmax))**2):
        raise ValueError('Argument below lower bound')

    if x < -1:
        z = sqrt(-x)
        c0 = cosh(z)
        c1 = sinh(z)/z
        c2 = (1 - c0)/x
        c3 = (1 - c1)/x
    elif x > 1:
        z = sqrt(x)
        c0 = cos(z)
        c1 = sin(z)/z
        c2 = (1 - c0)/x
        c3 = (1 - c1)/x
    else:
        c3 = sum(power(x, arange(0, trunc-1)) / odd_factorials)
        c2 = sum(power(x, arange(0, trunc-1)) / even_factorials)
        c1 = 1 - x*c3
        c0 = 1 - x*c2
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
    
    x = dt/bq
    x = bracket(x, -bound, bound)
    kfun = kepler(x)
    
    if dt < 0:
        upper = 0
        lower = x
        
        while kfun > dt:
            upper = lower
            lower = lower * 2
            oldx = x
            x = bracket(lower, -bound, bound)
            if x == oldx:
                raise ValueError(f'The input delta time (dt) has a value of {dt}.'
                                 'This is beyond the range of DT for which we '
                                 'can reliably propagate states. The limits for '
                                 'this GM and initial state are from {kepler(-bound)} to '
                                 '{kepler(bound)}.')
            kfun = kepler(x)
        
    elif dt > 0:
        lower = 0
        upper = x
        
        while kfun < dt:
            lower = upper
            upper = upper * 2
            oldx = x
            x = bracket(upper, -bound, bound)
            if x == oldx:
                raise ValueError(f'The input delta time (dt) has a value of {dt}.'
                                 'This is beyond the range of DT for which we '
                                 'can reliably propagate states. The limits for '
                                 'this GM and initial state are from {kepler(-bound)} to '
                                 '{kepler(bound)}.')
            kfun = kepler(x)
    else:
        return position, velocity
    
    x = min(upper, max(lower, (lower + upper)/2))
    
    lcount = 0
    mostc = 1000
    
    while (lower < x < upper) and (lcount < mostc):
        kfun = kepler(x)
        if kfun > dt:
            upper = x
        elif kfun < dt:
            lower = x
        else:
            # TODO: when vectorizing make sure to use np.copy here:
            upper = lower = x
        
        if mostc > 64: # this is the hard coded value of the variable maxbit in the original fortran function
            if upper != 0 and lower != 0:
                mostc = 64
                lcount = 0
        
        x = min(upper, max(lower, (upper + lower)/2))
        lcount += 1
    
    c0, c1, c2, c3 = stumpff(f*x*x)
    br = br0*c0 + x*(b2rv*c1 + x*(bq*c2))
    
    pc = 1 - qovr0 * x**2 * c2
    vc = dt - bq * x**3 * c3
    pcdot = -qovr0 / br * x * c1
    vcdot = 1 - bq / br * x**2 * c2
    
    position_prop = pc*position + vc*velocity
    velocity_prop = pcdot*position + vcdot*velocity
    
    return position_prop, velocity_prop
