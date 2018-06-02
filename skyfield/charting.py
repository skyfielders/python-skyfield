"""Routines to help draw star charts."""

from .starlib import Star

def plot_stars(catalog, observer, project, ax):
    xmin, xmax = ax.get_xlim()
    ymin, ymax = ax.get_xlim()
    lim = max(abs(xmin), abs(xmax), abs(ymin), abs(ymax)) * 2.0
    lims = (-lim, lim)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect('equal')

    c = catalog
    c = c[c['magnitude'] < 6.5]
    s = Star(ra_hours=c.ra_hours, dec_degrees=c.dec_degrees)

    o = observer[0]
    spos = o.observe(s)
    x, y = project(spos)
    size = 1.0
    size = 0.1 + c['magnitude'].max() - c['magnitude']
    size **= 2.5
    ax.scatter(x, y, s=size, c='k')
