#!/usr/bin/env python3

"""
Choosing a source for ∆T
========================

1. `ftp://ftp.iers.org/products/eop/rapid/standard/`_
2. `ftp://ftp.iers.org/products/eop/rapid/standard/csv/`_

* The CSV files are slightly larger, both compressed and uncompressed.
* For the position of an object at the celestial equator,
  1 second of time = 15 arcseconds.
  1 millisecond of time = 15 mas.

TODO:

[X] Build reader for long-term file.
[X] Build polynomial interpolation function.
[ ] Build translator between the two.
[ ] Build translator between finals2000A and interpolator.
[ ] Build combiner.
[ ] Save combined table in Skyfield?  Or too big with all those 0's?

"""
import sys
from time import time

import matplotlib.pyplot as plt
import numpy as np
from numpy import concatenate, maximum as clip_upper, nan
from skyfield import functions, timelib
from skyfield.api import load, Loader
from skyfield.data import iers
from skyfield.data.earth_orientation import parse_S15_table

def cat(*args):
    return concatenate(args)

inf = float('inf')

def main(argv):
    try_out_new_class()
    #compare_splines_to_finals2000_error_bars()
    #try_adjusting_spline()
    #try_solving_spline()
    #try_out_different_interpolation_techniques()

def try_out_new_class():
    f = load.open('finals2000A.all')
    mjd_utc, dut1 = iers.parse_dut1_from_finals_all(f)
    delta_t_recent, leap_dates, leap_offsets = (
        iers._build_timescale_arrays(mjd_utc, dut1)
    )
    daily_tt, daily_delta_t = delta_t_recent

    ts = load.timescale()
    #daily_t = ts.tt_jd(daily_tt)

    url = 'http://astro.ukho.gov.uk/nao/lvm/Table-S15.2020.txt'
    with load.open(url) as f:
        names, columns = parse_S15_table(f)
    #s15 = dict(zip(names, columns))
    i, start_year, end_year, a0, a1, a2, a3 = columns

    parabola = stephenson_morrison_hohenkerk_2016_parabola
    print(move_spline_endpoints(100, 200, parabola.left, parabola.right,
                                *parabola.coefficients))
    print(move_spline_endpoints(0, 1, parabola.left, parabola.right,
                                *parabola.coefficients))

    b3, b2, b1, b0 = move_spline_endpoints(
        -721, -720, parabola.left[0], parabola.right[0],
        *[c[0] for c in parabola.coefficients],
    )
    print(b3, b2, b1, b0)

    c3, c2, c1, c0 = move_spline_endpoints(
        2019, 2020, parabola.left[0], parabola.right[0],
        *[c[0] for c in parabola.coefficients],
    )
    print(c3, c2, c1, c0)

    long_term_curve = Splines(
        cat([-721], start_year, [2019]),
        cat([-720], end_year, [2020]),
        cat([b3], a3, [c3]),
        cat([b2], a2, [c2]),
        cat([b1], a1, [c1]),
        cat([b0], a0, [c0]),
    )
    #     s15['K_i'], s15['K_{i+1}'],
    #     s15['a_3'], s15['a_2'], s15['a_1'], s15['a_0'],
    # )

    daily_tt, daily_delta_t = extend(ts, daily_tt, daily_delta_t, long_term_curve, 180)

    delta_t_function = DeltaT(daily_tt, daily_delta_t, long_term_curve)
    #year = np.arange(1980, 2010, 0.1)
    #year = np.arange(-720, 2000)
    #year = np.arange(-720, 2020)
    #year = np.arange(-800, 2020)

    fig, axes = plt.subplots(4, 1)
    axes = iter(axes)

    ax = next(axes)
    #t = ts.J(np.linspace(-2000, 2500))
    t = ts.J(np.linspace(-2000, 2150))
    do_plot(ax, t, delta_t_function)
    ax.plot(t.J, stephenson_morrison_hohenkerk_2016_parabola(t.J), linestyle='--')

    #do_plot(next(axes), ts.J(np.linspace(-800, 2020)), delta_t_function)
    do_plot(next(axes), ts.J(np.linspace(-730, -710)), delta_t_function)

    days = 360

    tt = np.linspace(daily_tt[0] - days, daily_tt[0] + days)
    do_plot(next(axes), ts.tt_jd(tt), delta_t_function)

    tt = np.linspace(daily_tt[-1] - days, daily_tt[-1] + days)
    do_plot(next(axes), ts.tt_jd(tt), delta_t_function)

    # ax.plot(t.J, delta_t)
    # ax.grid()
    # ax2.plot(t.J[1:], np.diff(delta_t))
    # ax2.grid()
    fig.tight_layout()
    fig.savefig('tmp.png')

def do_plot(ax, t, delta_t_function):
    ax.plot(t.J, delta_t_function(t))
    ax.plot(t.J, delta_t_function.long_term_curve(t.J))
    ax.grid()

def compare_splines_to_finals2000_error_bars():
    #url = 'http://astro.ukho.gov.uk/nao/lvm/Table-S15-v18.txt'
    url = 'http://astro.ukho.gov.uk/nao/lvm/Table-S15.2020.txt'
    with load.open(url) as f:
        names, columns = parse_S15_table(f)

    i, start_year, end_year, a0, a1, a2, a3 = columns

    f = load.open('finals2000A.all')
    mjd_utc, dut1 = iers.parse_dut1_from_finals_all(f)
    delta_t_recent, leap_dates, leap_offsets = (
        iers._build_timescale_arrays(mjd_utc, dut1)
    )

    print('Size of IERS table:', delta_t_recent.shape)
    print('Number of splines:', i.shape)

    #year = [-720, 400, 700]
    #year = np.arange(-720, 2010)
    #year = np.arange(1800, 2010)
    #year = np.arange(1980, 2010)
    year = np.arange(1980, 2010, 0.1)
    interpolate = Splines(start_year, end_year, a3, a2, a1, a0)
    s15_curve = interpolate(year)

    finals_tt, finals_delta_t = delta_t_recent
    ts = load.timescale()
    t = ts.utc(year)
    tt = t.tt

    interpolate = Splines(
        finals_tt[:-1],
        finals_tt[1:],
        finals_delta_t[1:] - finals_delta_t[:-1],
        finals_delta_t[:-1],
    )

    T0 = time()
    finals_curve = interpolate(tt)
    print(time() - T0, 's for interpolate()')

    T0 = time()
    finals_curve2 = np.interp(tt, finals_tt, finals_delta_t)
    print(time() - T0, 's for interp()')

    assert (finals_curve == finals_curve2).all()

    diff = max(abs(s15_curve - finals_curve))
    print('Max difference (seconds, arcseconds):', diff, diff * 15)

    compare_interpolations(t, finals_curve, s15_curve)

def compare_interpolations(t, f1, f2):
    dt1 = f1#(t)
    dt2 = f2#(t)

    diff = dt2 - dt1
    print('Biggest difference:', np.max(np.abs(diff)), 'seconds')

    fig, (ax, ax2) = plt.subplots(2, 1)
    ax.plot(t.J, dt1, linestyle='--')
    ax.plot(t.J, dt2)
    ax.grid()

    ax2.plot(t.J, diff)
    ax2.grid()

    fig.savefig('tmp.png')

def extend(ts, daily_tt, daily_delta_t, long_term_curve, days):
    t_left = ts.tt_jd(daily_tt[0] - days)
    t_right = ts.tt_jd(daily_tt[-1] + days)
    y_left = long_term_curve(t_left.J)
    y_right = long_term_curve(t_right.J)
    return (cat([t_left.tt], daily_tt, [t_right.tt]),
            cat([y_left], daily_delta_t, [y_right]))

class DeltaT(object):
    def __init__(self, daily_tt, daily_delta_t, long_term_curve):
        self.daily_tt = daily_tt
        self.daily_delta_t = daily_delta_t
        self.long_term_curve = long_term_curve

    def __call__(self, t):
        delta_t = np.interp(t.tt, self.daily_tt, self.daily_delta_t, nan, nan)
        [nan_indexes] = np.nonzero(np.isnan(delta_t))  # Or np.argwhere()?
        #print('nan_indexes:', nan_indexes)
        if len(nan_indexes):
            delta_t[nan_indexes] = self.long_term_curve(t.J[nan_indexes])
        return delta_t

def _a(a):
    return a if hasattr(a, 'shape') else np.array((a,))

class Splines(object):
    def __init__(self, left, right, *coefficients):
        self.left = left = _a(left)
        self.right = right = _a(right)
        self._width = right - left
        self.n = np.arange(len(left))
        self.coefficients = [_a(c) for c in coefficients]

    def __call__(self, x):
        i = np.interp(x, self.left, self.n)
        i = i.astype(int)
        t = (x - self.left[i]) / self._width[i]
        coefficients = iter(self.coefficients)
        value = next(coefficients)[i]
        for c in coefficients:
            value *= t
            value += c[i]
        return value

stephenson_morrison_hohenkerk_2016_parabola = Splines(
    1825.0, 1925.0, 0.0, 32.5, 0.0, -320.0)

def move_spline_endpoints(new_left, new_right, old_left, old_right,
                          a3, a2, a1, a0):
    k0, k1 = old_left, old_right
    j0, j1 = new_left, new_right

    u0 = a1/(-k0 + k1)
    u1 = j0*u0
    u2 = k0**2
    u3 = k1**2
    u4 = 2*k0
    u5 = a2/(-k1*u4 + u2 + u3)
    u6 = k0**3
    u7 = 3*k0
    u8 = 3*u2
    u9 = a3/(k1**3 + k1*u8 - u3*u7 - u6)
    u10 = j0**3*u9
    u11 = u4*u5
    u12 = j0*u11
    u13 = u8*u9
    u14 = j0*u13
    u15 = j0**2
    u16 = u15*u5
    u17 = u15*u9
    u18 = u16 - u17*u7
    u19 = j0*j1
    u20 = 2*u19*u5
    u21 = 3*u10
    u22 = 6*k0
    u23 = u19*u22*u9
    u24 = j1*u17
    u25 = 3*u24
    u26 = j1**2
    u27 = u26*u9
    u28 = 3*j0*u27

    b0 = a0 - k0*u0 + u1 + u10 - u12 + u14 + u18 + u2*u5 - u6*u9
    b1 = j1*u0 - j1*u11 + j1*u13 - u1 + u12 - u14 - 2*u16 + u17*u22 + u20 - u21 - u23 + u25
    b2 = u18 - u20 + u21 + u23 - 6*u24 + u26*u5 - u27*u7 + u28
    b3 = j1**3*u9 - u10 + u25 - u28
    return b3, b2, b1, b0

def try_solving_spline():
    import sympy as sy
    sy.init_printing()

    a0, a1, a2, a3, k0, k1, j0, j1, new_t = sy.symbols(
        'a0, a1, a2, a3, k0, k1, j0, j1, new_t')

    # Q: How much simpler is it if we only move one end?
    #j0 = k0  # A: Wow, much simpler!  30+ ops instead of 80+
    #j1 = k1  # A: GADS, not much simpler at all, still 80+ operations.

    years = new_t * (j1 - j0) + j0
    old_t = (years - k0) / (k1 - k0)
    d = (((a3 * old_t + a2) * old_t) + a1) * old_t + a0

    #d = sy.factor(d)
    #d = sy.expand(d)
    #d = sy.simplify(d)
    d = sy.expand(d)

    d = sy.collect(d, new_t)
    b0 = d.coeff(new_t, 0)
    b1 = d.coeff(new_t, 1)
    b2 = d.coeff(new_t, 2)
    b3 = d.coeff(new_t, 3)

    commons, outputs = sy.cse(
        [b0, b1, b2, b3],
        sy.numbered_symbols('u'),
        #optimizations='basic',
    )
    n = 0
    for symbol, expr in commons:
        n += sy.count_ops(expr)
        print(symbol, '=', expr)
    print()
    for i, expr in enumerate(outputs):
        n += sy.count_ops(expr)
        print('b{} = {}'.format(i, expr))
    print('Total operations: {}'.format(n))

def try_out_different_interpolation_techniques():
    url = 'http://astro.ukho.gov.uk/nao/lvm/Table-S15-v18.txt'
    with load.open(url) as f:
        names, columns = parse_S15_table(f)
    i, start_year, end_year, a0, a1, a2, a3 = columns
    report = []
    print('Table start and end years:', start_year[0], end_year[-1])

    # Range of years to plot.

    #y = np.arange(start_year[0], end_year[-1] + 0.1, 0.01)
    y = np.arange(start_year[0], end_year[-1] + 0.1, 0.1)
    #y = np.arange(end_year[-1] - 30.0, end_year[-1] + 0.1, 0.01)

    # Skyfield original tables.

    ts = Loader('ci').timescale(builtin=False)
    print('Old shape:', ts.delta_t_table.shape)
    t = ts.J(y)

    T0 = time()
    t.delta_t
    report.append((time() - T0, 's for old interpolation tables'))

    # For perspective.

    T0 = time()
    t.M
    report.append((time() - T0, 's to compute N P B'))

    # Skyfield IERS table-driven interpolation.

    ts = load.timescale()
    print('IERS shape:', ts.delta_t_table.shape)
    t = ts.J(y)

    t.delta_t
    del t.delta_t
    T0 = time()
    old_delta_t = t.delta_t
    report.append((time() - T0, 's for IERS table interpolated delta_t'))

    # New Parabola.

    c1825 = (y - 1825.0) / 100.0
    t0 = time()
    delta_t_parabola = -320.0 + 32.5 * c1825 * c1825
    report.append((time() - t0, 's to compute ∆T with 2018 parabola'))

    # New 2018 splines.

    indexes = np.arange(len(start_year))
    #print(start_year[:10])

    T0 = time()

    #i = np.searchsorted(start_year, y, 'right') - 1
    i = np.interp(y, start_year, indexes)
    t = i
    i = i.astype(int)
    t %= 1.0
    # y0 = start_year[i]
    # y1 = end_year[i]
    # t = (y - y0) / (y1 - y0)
    delta_t = ((a3[i] * t + a2[i]) * t + a1[i]) * t + a0[i]

    report.append((time() - T0, 's to compute ∆T with 2018 splines'))

    # New 2018 splines combined with IERS data.

    arrays = functions.load_bundled_npy('iers.npz')
    iers_tt, iers_delta_t = arrays['delta_t_recent']
    iers_y = (iers_tt - 1721045.0) / 365.25

    cutoff_index = np.searchsorted(end_year, iers_y[0])
    cutoff_year = end_year[cutoff_index - 1]

    print(cutoff_index, len(end_year))
    print(cutoff_year)
    # print(end_year[:cutoff_index])

    def to_tt(year):
        return year * 365.25 + 1721045.0

    # : i, start_year, end_year, a0, a1, a2, a3 :
    z = np.zeros_like(iers_tt)
    #iers_end_year =

    c = cutoff_index
    start_tt = cat([to_tt(start_year[:c]), to_tt(end_year[-1:]), iers_tt[:-1]])
    end_tt = cat(to_tt(end_year[:c]), iers_tt)
    a0 = cat(a0[:c], [0], iers_delta_t[:-1])
    a1 = cat(a1[:c], iers_delta_t - a0[-len(iers_delta_t):])
    a2 = cat(a2[:c], z)
    a3 = cat(a3[:c], z)

    # Try using combined splines.

    tt = to_tt(y)
    f_index = np.arange(len(start_tt))

    A = np.searchsorted(start_tt, tt, 'right') - 1
    T0 = time()
    A = np.searchsorted(start_tt, tt, 'right') - 1
    report.append((time() - T0, 'trial A'))

    B = np.searchsorted(start_tt, tt, 'left')
    T0 = time()
    B = np.searchsorted(start_tt, tt, 'left')
    report.append((time() - T0, 'trial B'))

    print('sizes:', tt.shape, start_tt.shape, f_index.shape)
    interp = np.interp
    C = interp(tt, start_tt, f_index, 1.0, 1.0)
    T0 = time()
    C = interp(tt, start_tt, f_index, 1.0, 1.0)
    report.append((time() - T0, 'trial C'))
    C1 = C.astype(int)

    assert A.shape == B.shape == C.shape

    print(A[:4])
    print(B[:4])
    print(C[:4])
    print(C1[:4])

    a0123 = np.array([a0, a1, a2, a3]).T
    print(a0123.shape)

    I1 = np.interp(tt, start_tt, f_index)

    T0 = time()
    #i = np.searchsorted(start_tt, tt, 'right') - 1
    I1 = np.interp(tt, start_tt, f_index)
    i = I1.astype(int)
    # I2 = I1 % 1.0
    # tt0 = start_tt[i]
    # tt1 = end_tt[i]
    # t = (tt - tt0) / (tt1 - tt0)
    # print('t ', t[:4])
    # print('I2', I2[:4])
    # t = I1 % 1.0
    I1 %= 1.0
    t = I1

    #print(a0123[i].shape)
    # a0, a1, a2, a3 = a0123[i].T
    # experimental_delta_t = ((a3 * t + a2) * t + a1) * t + a0

    experimental_delta_t = ((a3[i] * t + a2[i]) * t + a1[i]) * t + a0[i]
    #experimental_delta_t = a0[i]

    report.append((time() - T0, 's to compute ∆T with combined spline table'))

    experimental_delta_t

    # Hybrid approach: splines for years outside IERS range, but simple
    # linear interpolation inside IERS table.

    print(y)
    y_tt = to_tt(y)
    iers_mask = (iers_tt[0] <= y_tt) & (y_tt <= iers_tt[-1])
    print(iers_mask.shape, sum(iers_mask))

    # generate index: subtract floor, build mask, turn matches into ints
    # and 0.0-0.99 remainders, then do direct indexing into table.
    # for ~mask, do splines, maybe having fake spline where table is?
    # paste them together

    # T0 = time()
    # report.append((time() - T0, 's to compute ∆T with combined spline table'))

    # The plot.

    fig, (ax1, ax2, ax3) = plt.subplots(3, 1)

    ax = ax1

    ax.plot(y, delta_t_parabola, label='2018 parabola')
    ax.plot(y, delta_t, label='2018 splines ∆T')
    ax.plot(y, old_delta_t, label='Old Skyfield ∆T')

    #ax.set_xlim(-720, -718)
    #ax.set_ylim(-400, 200)
    ax.set_xlim(1550, 2020)
    ax.set_ylim(-400, 200)

    ax.legend()
    ax.grid()

    ax = ax2

    i = (y >= 1973)
    ax.plot(y[i], delta_t[i] - old_delta_t[i])
    ax.grid()
    ax.set(ylabel='2018 splines - IERS')

    ax = ax3

    end_year = 1973.1
    i = (y >= 1972.9) & (y <= end_year)
    ax.plot(y[i], delta_t[i])

    i = (iers_y <= end_year)
    ax.plot(iers_y[i], iers_delta_t[i], '.')
    ax.grid()
    #ax.set(ylabel='2018 splines - IERS')

    fig.savefig('tmp.png')

    for args in report:
        print(*args)

if __name__ == '__main__':
    main(sys.argv[1:])
