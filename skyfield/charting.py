"""Routines to help draw star charts."""

#import numpy as np
from IPython.display import HTML
from matplotlib.animation import FuncAnimation
from .starlib import Star

def _plot_stars(catalog, observer, project, ax, mag1, mag2, margin=1.25):
    """Experiment in progress, hence the underscore; expect changes."""

    art = []

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
    lim = max(abs(xmin), abs(xmax), abs(ymin), abs(ymax)) * margin
    lims = (-lim, lim)
    ax.set_xlim(lims)
    ax.set_ylim(lims)
    ax.set_aspect('equal')

    o = observer[0]

    c = catalog
    c = c[c['magnitude'] <= mag1]
    print('First star group:', len(c))
    s = Star(ra_hours=c.ra_hours, dec_degrees=c.dec_degrees)
    spos = o.observe(s)
    x, y = project(spos)
    scale = 2.0
    size = ((mag1 - c['magnitude']) * scale) ** 2.0
    art.append(ax.scatter(x, y, s=size, c='k'))

    c = catalog
    c = c[c['magnitude'] > mag1]
    c = c[c['magnitude'] <= mag2]
    print('Second star group:', len(c))
    s = Star(ra_hours=c.ra_hours, dec_degrees=c.dec_degrees)
    spos = o.observe(s)
    x, y = project(spos)
    m = (mag2 - c['magnitude']) / (mag2 - mag1)
    # Note that "gray_r" is white for 0.0 and black for 1.0
    art.append(ax.scatter(
        x, y, s=1.0,
        c=0.1 + 0.8 * m, cmap='gray_r', vmin=0.0, vmax=1.0,
    ))
    return art

X = []

class _Animation(FuncAnimation):
    patcher = None

    def _draw_next_frame(self, framedata, blit):
        #patcher.__exit__()

        #print('------------ _draw_next_frame()')
        #print('====>', self.AX.artists)
        #print('====>', self.AX.axes)
        #print('====>', self.AX.axes.lines)
        #print('====>', self.AX.axes.axes)  just itself
        #print('====>', self.AX.axes.get_children())
        #printout(self._fig)
        # if X:
        #     self.AX.artists[:] = X
        #     del X[:]

        blit = True    # override matplotlib's refusal to blit in save()
        super(_Animation, self)._draw_next_frame(framedata, blit)

        # X[:] = self.AX.artists
        # self.AX.artists[:] = []
        #print('------------ DONE')

        if self.patcher is None:
            self.patcher = patch('matplotlib.backends.backend_agg'
                            '.FigureCanvasAgg.draw', draw)
            self.patcher.__enter__()

def draw(self):
    pass

from unittest.mock import patch

def _animate(h, e, project, p):
    #from numpy import sin
    import matplotlib.pyplot as plt

    x, y = project(p)

    fig, ax = plt.subplots()
    print([name for name in dir(fig) if 'size' in name])
    fig.set_size_inches(9, 9)
    #(artist,) = ax.plot(x, sin(x))

    blue_art, = ax.plot(x, y, color='b', alpha=0.0)

    # print('blue_art:', id(blue_art))
    # print('planet art:', id(planet))

    _plot_stars(h, e, project, ax, 6.0, 8.0, 0.8)

    planet = None

    def init():
        nonlocal planet
        planet, = ax.plot(x[100], y[100], 'ro')
        #return ()
        #return a + [blue_art]
        #return a
        return (planet,)

    from time import time

    def update(i):
        # print(a)
        # print(fig.artists)

        # for artist in a:
        #     if a in fig.artists:
        #         print('=================== REMOVING')
        #         fig.artists.remove(a)

        # if t0[0] is not None:
        #     print(' {} seconds'.format(time() - t0[0]))
        t0[0] = time()
        # print('Frame', i, end=" ")
        #imshow(copy)
        #fig.canvas.restore_region(copy)
        planet.set_xdata(x[i])
        planet.set_ydata(y[i])
        #fig.canvas.blit(ax.bbox)
        #return ()
        return (planet,)

    #fig.canvas.draw()
    #fig.draw()
    #artist.axes.draw_artist(artist)

    # print(planet.figure is fig)

    #copy = artist.figure.canvas.copy_from_bbox(artist.axes.bbox)

    #anim = FuncAnimation(fig, update_frame, frames=10)
    t0 = [time()]

    anim = _Animation(fig, update, frames=len(x),
                      blit=True, init_func=init,
                      interval=50)
    plt.close()
    #anim.AX = ax
    html = anim.to_html5_video()

    anim.patcher.__exit__()
    return HTML(html)

class T(object):
    from time import time
    def __enter__(self):
        self.t0 = self.time()
    def __exit__(self, *args):
        print('<<< %s >>>' % (self.time() - self.t0))

def printout(thing):
    print()
    for line in debug(thing):
        print(line)
    print()

def debug(thing, indent=0):
    a = 'A' if thing.get_animated() else '-'
    yield '{:<{}}{} {} {}'.format('', indent, hex(id(thing)), a, repr(thing))
    indent += 2
    for child in thing.get_children():
        yield from debug(child, indent)
