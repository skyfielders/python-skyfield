"""Routines to help draw star charts."""

#import numpy as np
from .starlib import Star

def plot_stars(catalog, observer, project, mag1, mag2, ax):


    # from astropy import wcs
    # w = wcs.WCS(naxis=2)
    # w.wcs.crpix = [-234.75, 8.3393]
    # w.wcs.cdelt = np.array([-0.066667, 0.066667])
    # w.wcs.crval = [0, -90]
    # w.wcs.ctype = ["RA---AIR", "DEC--AIR"]
    # w.wcs.set_pv([(2, 1, 45.0)])

    # import matplotlib.pyplot as plt

    # plt.subplot(projection=wcs)
    # #plt.imshow(hdu.data, vmin=-2.e-5, vmax=2.e-4, origin='lower')
    # plt.grid(color='white', ls='solid')
    # plt.xlabel('Galactic Longitude')
    # plt.ylabel('Galactic Latitude')

    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_xlim()
    lim = max(abs(xmin), abs(xmax), abs(ymin), abs(ymax)) * 1.25
    lims = (-lim, lim)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect('equal')

    o = observer[0]

    c = catalog
    c = c[c['magnitude'] <= mag1]
    s = Star(ra_hours=c.ra_hours, dec_degrees=c.dec_degrees)
    spos = o.observe(s)
    x, y = project(spos)
    mmin = c['magnitude'].min()
    mmax = c['magnitude'].max()
    m = (mmax - c['magnitude']) / (mmax - mmin)
    size = 1.0 + 10.0 * m
    size **= 1.8
    ax.scatter(x, y, s=size, c='k')

    c = catalog
    c = c[c['magnitude'] >= mag1]
    c = c[c['magnitude'] <= mag2]
    s = Star(ra_hours=c.ra_hours, dec_degrees=c.dec_degrees)
    spos = o.observe(s)
    x, y = project(spos)
    m = (mag2 - c['magnitude']) / (mag2 - mag1)
    size = 0.1 + m * 0.9
    # Note that "gray_r" is white for 0.0 and black for 1.0
    ax.scatter(x, y, s=1, c=m, cmap='gray_r', vmin=0.0, vmax=1.0)
