#!/usr/bin/env python3
"""Build Skyfield's internal table of constellation boundaries.

See:

https://iopscience.iop.org/article/10.1086/132034/pdf
http://cdsarc.u-strasbg.fr/viz-bin/Cat?VI/42

"""
import argparse
import sys
from numpy import array, searchsorted
from skyfield import api

URL = 'http://cdsarc.u-strasbg.fr/ftp/VI/42/data.dat'

def main():
    with api.load.open(URL) as f:
        lines = list(f)

    unique_ra = set()
    unique_dec = set()
    fracs = set()
    boundaries = []

    for line in lines:
        fields = line.split()
        ra_low = extend(fields[0])
        ra_up = extend(fields[1])
        de_low = extend(fields[2])
        const = fields[3].decode('ascii')

        print(ra_low, const)

        #print(ra_int(ra_low))

        #fracs.add(fields[0].split(b'.')[1])
        unique_ra.add(ra_low)
        unique_ra.add(ra_up)
        unique_dec.add(de_low)
        fracs.add(const)
        boundaries.append([ra_low, ra_up, de_low, const])

    print(sorted(fracs))
    print('constellations:', len(fracs))
    print('unique_ra:', len(unique_ra))
    print('unique_dec:', len(unique_dec))

    sorted_consts = array(sorted(fracs))
    sorted_ra = array(sorted(unique_ra))
    sorted_dec = array(sorted(unique_dec))

    assert sorted_ra[0] == 0
    assert sorted_ra[-1] == 24

    assert sorted_dec[0] == -90
    assert sorted_dec[-1] == 88

    sorted_ra = sorted_ra[1:-1]

    print('bytes', sorted_ra.nbytes)
    print('bytes', sorted_dec.nbytes)

    #grid = [[5] * len(unique_dec)] * len(unique_ra)
    #grid = array(grid, 'i1')

    row = [-128] * len(sorted_ra + 1)
    grid = []
    i = 0
    de = -90.0
    for ra_low, ra_up, de_low, const in boundaries[::-1]:
        if de_low > de:
            grid.append(row)
            row = list(row)
            de = de_low
        i0 = searchsorted(sorted_ra, ra_low)
        i1 = searchsorted(sorted_ra, ra_up)
        c = searchsorted(sorted_consts, const)
        for j in range(i0, i1):
            row[j] = c

    grid.append(row)
    grid.append(row)
    #grid = grid[::-1]
    grid = array(grid, 'i1').T

    print(sorted_consts[57])
    print(grid)

    print('shape', grid.shape)
    print('bytes', grid.nbytes)

    for ra, dec in [(0, 0), (0.1, 0.1),
                    (16, 80), (16, 90), (16, -90), (24, 360),
                    ([0, 16], [0, 80])]:
        c = compute_constellation(ra, dec, sorted_ra, sorted_dec,
                                  sorted_consts, grid)
        print('=', ra, dec, c)

def compute_constellation(ra, dec, sorted_ra, sorted_dec, sorted_consts, grid):
    i = searchsorted(sorted_ra, ra) - 1
    j = searchsorted(sorted_dec, dec)
    print("i,j", i, j)
    # print(ra, sorted_ra)
    return sorted_consts[grid[i, j]]

def extend(s):
    """Return a float for `s` extended to machine precision.

    Takes a string like '13.6667', extends it to >30 digits by
    duplicating the second-to-last character ('6' and not '7' because
    the '7' in many cases will have been rounded), and passes it to
    `float()`.

    """
    s = s[:-2] + 30 * s[-2:-1] + s[-1:]
    return float(s)

# Some discarded code that I might want to revive someday: how to grow
# and shrink a list of segments as new ones supersede old ones on the
# way down the sky.

def segment_experiment():
    assert insert_segment([0, 4, 7, 10], 0, 3) == [0, 3, 4, 7, 10]
    assert insert_segment([0, 4, 7, 10], 4, 7) == [0, 4, 7, 10]
    assert insert_segment([0, 4, 7, 10], 6, 9) == [0, 4, 6, 9, 10]
    assert insert_segment([0, 4, 7, 10], 7, 10) == [0, 4, 7, 10]
    assert insert_segment([0, 4, 7, 10], 0, 10) == [0, 10]
    assert insert_segment([0, 10], 4, 7) == [0, 4, 7, 10]
    assert insert_segment([], 4, 7) == [4, 7]

    segments = []
    n = 0
    for ra_low, ra_up, de_low in boundaries[::-1]:
        segments = insert_segment(segments, ra_low, ra_up)
        print(len(segments), end=' ')
        n += len(segments)
    print(n)

def insert_segment(ra_list, ra_low, ra_up):
    new = []
    i = 0
    while i < len(ra_list) and ra_list[i] < ra_low:
        new.append(ra_list[i])
        i += 1
    new.append(ra_low)
    new.append(ra_up)
    while i < len(ra_list) and ra_list[i] <= ra_up:
        i += 1
    while i < len(ra_list):
        new.append(ra_list[i])
        i += 1
    return new

if __name__ == '__main__':
    main()
