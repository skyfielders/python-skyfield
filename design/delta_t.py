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

Oh, look, they came out with a new paper!

https://royalsocietypublishing.org/doi/10.1098/rspa.2020.0776

Its data:

https://rs.figshare.com/collections/Supplementary_material_from_Addendum_2020_to_Measurement_of_the_Earth_s_rotation_720_BC_to_AD_2015_/5300925

http://astro.ukho.gov.uk/nao/lvm/Table-S15.2020.txt

Their long-term formula:

lod = +1.72 t − 3.5 sin(2π(t+0.75)/14)

TODO:

[ ] Build reader for long-term file.
[ ] Build polynomial interpolation function.
[ ] Build translator between the two.
[ ] Build translator between finals2000A and interpolator.
[ ] Build combiner.
[ ] Save combined table in Skyfield?  Or too big with all those 0's?
[ ] Can we make the IERS table more sparse without loosing too much precision?

"""
import sys
from time import time

import matplotlib.pyplot as plt
import numpy as np
from numpy import concatenate as cat
from skyfield import functions, timelib
from skyfield.api import load, Loader
from skyfield.data import iers

inf = float('inf')

def main(argv):
    try_cubic_splines_of_various_sparseness()
    #compare_splines_to_finals2000_error_bars()
    #try_adjusting_spline()
    #try_solving_spline()
    #try_out_different_interpolation_techniques()

def try_cubic_splines_of_various_sparseness():
    f = load.open('finals2000A.all')
    mjd_utc, dut1 = iers.parse_dut1_from_finals_all(f)
    delta_t, leap_dates, leap_offsets = (
        iers._build_timescale_arrays(mjd_utc, dut1)
    )

    table_tt, table_delta_t = delta_t

    from scipy import interpolate

    # x_points = [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9]
    # y_points = [12,14,22,39,58, 77,89,95,98,99]

    x_points, y_points = delta_t

    # if 1:
    #     extra = 0
    #     extended_x = x_points
    #     extended_y = y_points
    # else:
    #     extra = 1
    #     extended_x = cat([[-1], x_points, [10]])
    #     extended_y = cat([y_points[:1], y_points, y_points[-1:]])

    # t, c, k = interpolate.splrep(extended_x, extended_y)
    # pp = interpolate.PPoly.from_spline((t, c, k))

    # # Let SciPy do the interpolation.

    # #y = interpolate.splev(x, (t, c, k))
    # #y2 = pp(x)

    # # Do it ourselves.

    # print(pp)
    # print(pp.c)
    # print(pp.x)
    # # print(pp.c.shape, 'vs', len(x_points))
    # trim_n = 3 #+ extra
    # trim = slice(trim_n, -trim_n)
    # breaks = pp.x[trim]
    # coeffs = pp.c[:,trim]
    # print(breaks)
    # print(coeffs)

    step = 2
    dtp = build_splines(x_points[::step], y_points[::step])

    x = np.linspace(x_points[0], x_points[-2-step], 1000)
    y2 = dtp(x)

    np.savez_compressed(
        'tmp.npz',
        b=dtp._b,
        c=dtp._c,
    )

    diff = max(abs(y_points[:-1-step] - dtp(x_points[:-1-step])))

    print('Max difference between y and y2:\n',
          diff, 'seconds =', diff * 15, 'arcseconds')

    return

    fig, ax = plt.subplots()
    ax.plot(x_points, y_points, linestyle='--')
    # ax.plot(x, y)
    ax.plot(x, y2)
    ax.grid()
    fig.savefig('tmp.png')

def build_splines(x, y):
    from scipy import interpolate
    t, c, k = interpolate.splrep(x, y)
    pp = interpolate.PPoly.from_spline((t, c, k))

    # Let SciPy do the interpolation.

    #y = interpolate.splev(x, (t, c, k))
    #y2 = pp(x)

    # Do it ourselves.

    print(pp)
    print(pp.c)
    print(pp.x)
    # print(pp.c.shape, 'vs', len(x_points))
    trim_n = 3 #+ extra
    trim = slice(trim_n, -trim_n)
    breaks = pp.x[trim]
    coeffs = pp.c[:,trim]
    print(breaks)
    print(coeffs)

    return DeltaTPolynomial(breaks, coeffs)

class DeltaTPolynomial(object):
    def __init__(self, breaks, coefficients):
        self._b = breaks
        self._c = coefficients
        self._i = np.arange(len(breaks))

    def __call__(self, x):
        i = np.interp(x, self._b, self._i)
        i, t = divmod(i, 1.0)  # t: polynomial parameter, not raw time
        i = i.astype(int)
        across = x - self._b[i]  # TODO: long-term polynomials scale to [0,1)
        c = iter(self._c)
        y = next(c)[i]
        for ci in c:
            y *= across
            y += ci[i]
        return y

def compare_splines_to_finals2000_error_bars():
    url = 'http://astro.ukho.gov.uk/nao/lvm/Table-S15-v18.txt'
    with load.open(url) as f:
        columns = load_table_S15(f)
    i, start_year, end_year, a0, a1, a2, a3 = columns

    f = load.open('finals2000A.all')
    mjd_utc, dut1 = iers.parse_dut1_from_finals_all(f)
    delta_t_recent, leap_dates, leap_offsets = (
        iers._build_timescale_arrays(mjd_utc, dut1)
    )

    print(delta_t_recent.shape)
    print(i.shape)

    #year = [-720, 400, 700]
    #year = np.arange(-720, 2010)
    #year = np.arange(1800, 2010)
    #year = np.arange(1980, 2010)
    year = np.arange(1980, 2010, 0.1)
    s15_curve = interpolate(year, start_year, a3, a2, a1, a0)

    finals_tt, finals_delta_t = delta_t_recent
    ts = load.timescale()
    t = ts.utc(year)
    tt = t.tt
    print(tt)
    print(finals_tt)
    finals_curve = interpolate(
        tt, finals_tt,
        finals_delta_t[1:] - finals_delta_t[:-1],
        finals_delta_t[:-1],
    )

    print(finals_curve)

    diff = max(abs(s15_curve - finals_curve))
    print('Max difference between long-term splines and finals2000A.all:\n',
          diff, 'seconds =', diff * 15, 'arcseconds')

    if 1:
        fig, (ax, ax2) = plt.subplots(2, 1)
        ax.plot(year, s15_curve, label='label', linestyle='--')
        ax.plot(year, finals_curve, label='label')
        ax.grid()

        ax2.plot(year, finals_curve - s15_curve, label='label')
        ax2.grid()

        #plt.legend()
        fig.savefig('tmp.png')

def interpolate(t, t_barriers, c, *coefficients):
    i_array = np.arange(len(t_barriers))
    i = np.interp(t, t_barriers, i_array)
    print(i)
    # x = i % 1.0
    # print(x)
    i = i.astype(int)
    print(i)
    left = t_barriers[i]
    right = t_barriers[i+1]
    print(left, right)
    little_t = (t - left) / (right - left)
    print(little_t)
    value = c[i]
    for c in coefficients:
        value *= little_t
        value += c[i]
    return value

def try_adjusting_spline():
    k0, k1 = -720.0, -100.0
    a0, a1, a2, a3 = 20371.848, -9999.586, 776.247, 409.160
    j0, j1 = -800, -700

    # Problem: generate new b0..b3 for j0,j1.

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

    print(a0, a1, a2, a3)
    print(b0, b1, b2, b3)

    # Try out the new parameters:

    years = np.array([-720, -710, -700])

    t = (years - k0) / (k1 - k0)
    d = (((a3 * t + a2) * t) + a1) * t + a0
    print(d)

    t = (years - j0) / (j1 - j0)
    d = (((b3 * t + b2) * t) + b1) * t + b0
    print(d)

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
        columns = load_table_S15(f)
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

    cat = np.concatenate
    c = cutoff_index
    start_tt = cat([to_tt(start_year[:c]), to_tt(end_year[-1:]), iers_tt[:-1]])
    end_tt = cat([to_tt(end_year[:c]), iers_tt])
    a0 = cat([a0[:c], [0], iers_delta_t[:-1]])
    a1 = cat([a1[:c], iers_delta_t - a0[-len(iers_delta_t):]])
    a2 = cat([a2[:c], z])
    a3 = cat([a3[:c], z])

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

    return

    # parser = argparse.ArgumentParser(description=put description here)
    # parser.add_argument('integers', metavar='N', type=int, nargs='+',
    #                     help='an integer for the accumulator')
    # parser.add_argument('--sum', dest='accumulate', action='store_const',
    #                     const=sum, default=max,
    #                     help='sum the integers (default: find the max)')
    # args = parser.parse_args(argv)
    # print(args.accumulate(args.integers))

    #draw_plot_comparing_USNO_and_IERS_data()

    f = load.open('finals2000A.all')
    mjd_utc, dut1 = iers.parse_dut1_from_finals_all(f)
    delta_t_recent, leap_dates, leap_offsets = (
        iers._build_timescale_arrays(mjd_utc, dut1)
    )

    ts = load.timescale()

    print('TT 2050:', ts.J(2050).tt)

    fig, ax = plt.subplots()

    # jd = mjd_utc + 2400000.5
    # t = ts.tt_jd(jd)
    # ax.plot(t.J, t.delta_t)

    t = ts.J(range(1600, 2150))
    ax.plot(t.J, t.delta_t, label='Skyfield')

    t = ts.J(range(2150, 2190))
    ms2004 = timelib.delta_t_formula_morrison_and_stephenson_2004(t.tt)
    ax.plot(t.J, ms2004, label='M&S 2004 long-term parabola')

    y = np.arange(2005, 2050)
    t = y - 2000
    ax.plot(y, 62.92 + 0.32217 * t + 0.005589 * t * t,
            label='M&S 2004 polynomial 2005-2050')

    y = np.arange(2050, 2150)
    ax.plot(y, -20 + 32 * ((y-1820)/100)**2 - 0.5628 * (2150 - y),
            label='M&S 2004 polynomial 2050-2150')
    # y = np.array([2050, 2149])
    # ax.plot(y, -20 + 32 * ((y-1820)/100)**2 - 0.5628 * (2150 - y))

    ax.set(xlabel='Year', ylabel='∆T')
    ax.grid()
    ax.legend()
    fig.savefig('tmp.png')

    is_2050_to_2150_polynomial_worth_it()

def load_table_S15(f):
    # http://astro.ukho.gov.uk/nao/lvm/Table-S15-v18.txt
    content = f.read()
    banner = b'- ' * 36 + b'-\n'
    sections = content.split(banner)
    table = np.loadtxt(sections[2].splitlines())
    return table.T

def is_2050_to_2150_polynomial_worth_it():
    # How different is the polynomial at:
    # https://eclipse.gsfc.nasa.gov/SEhelp/deltatpoly2004.html
    # from simply doing a linear interpolation from the endpoints?
    y = np.arange(2050, 2150)
    delta_t = -20 + 32 * ((y-1820)/100)**2 - 0.5628 * (2150 - y)
    linear_approximation = np.linspace(delta_t[0], delta_t[-1], len(delta_t))
    print('Maximum difference:', max(abs(delta_t - linear_approximation)), 's')

def draw_plot_comparing_USNO_and_IERS_data():
    f = load.open('finals.all')
    year, month, day, dut1 = iers.parse_dut1_from_finals_all(f)

    # auto-detect leap seconds
    leap_seconds = np.diff(dut1) > 0.9
    #x = np.append(leap_seconds > 0.9, False)
    #x = np.prepend(leap_seconds > 0.9, False)
    #x = np.concatenate([leap_seconds > 0.9, [False]]).astype(bool)
    x = np.concatenate([[False], leap_seconds, ])
    print(x)
    print(x, sum(x))
    print(year[x])
    print(month[x])
    print(day[x])

    delta_t = dut1 - np.cumsum(x)

    ts = load.timescale()
    t = ts.utc(year, month, day)
    y_new = dut1
    y_old = t.dut1



    # and: figure out error cost of interpolating weeks or months
    all_points = np.arange(len(delta_t))
    samples = all_points[::600] + 0.001
    print(f'{len(samples)} samples')

    for skip in range(1, 40):
        index = np.arange(0, len(delta_t), skip)
        excerpt = delta_t[index]
        interpolated = np.interp(all_points, index, excerpt)
        difference_seconds = abs(interpolated - delta_t).max()
        difference_arcseconds = difference_seconds / 24.0 * 360.0
        storage = 4 * len(excerpt)

        t0 = time()
        np.interp(samples, index, excerpt)
        duration = time() - t0

        fmt = '{:2} days  {:10.6f} s   {:10.6f} arcseconds  {} bytes  {:.6f} s'
        print(fmt.format(
            skip, difference_seconds, difference_arcseconds, storage, duration,
        ))
    #print(dut1)

    fig, ax = plt.subplots()
    ax.set(xlabel='Year', title='UT1 minus UTC')
    # ax.plot(t.J, y_old, '-', label='deltat.dat')
    # ax.plot(t.J, y_new, ',', label='finals.all')
    ax.plot(delta_t)
    ax.plot(np.diff(delta_t))
    ax.grid()
    plt.legend()
    fig.savefig('tmp.png')

    print(len(y_new))
    np.savez_compressed(
        'test.npz',
        time=t.J,
        ut1_minus_utc=y_new,
    )

"""
Example line from finals.all, for crafting our RE:

73 1 2 41684.00 I  0.120733 0.009786  0.136966 0.015902  I 0.8084178 0.0002710  0.0000 0.1916  P    44.969     .500     2.839     .300   .143000   .137000   .8075000      .000      .000  \n\\
"""

if __name__ == '__main__':
    main(sys.argv[1:])

