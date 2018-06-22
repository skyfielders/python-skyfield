"""Routines to help draw star charts."""

import numpy as np
from IPython.display import HTML
from matplotlib import rc_context, rcParams
from matplotlib.animation import FuncAnimation, FFMpegWriter, writers
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

    # Dim stars: points of with varying gray levels.

    c = catalog
    c = c[c['magnitude'] > mag1]
    c = c[c['magnitude'] <= mag2]
    print('Second star group:', len(c))
    c = c.sort_values('magnitude', ascending=False)
    s = Star(ra_hours=c.ra_hours, dec_degrees=c.dec_degrees)
    spos = o.observe(s)
    x, y = project(spos)
    m = (mag2 - c['magnitude']) / (mag2 - mag1)
    # Note that "gray_r" is white for 0.0 and black for 1.0
    art.append(ax.scatter(
        x, y, s=1.0,
        c=0.8 * m, cmap='gray_r', vmin=0.0, vmax=1.0,
    ))

    # Bright stars: black circles of varying radius, surrounded by a
    # white gap in case stars are touching.  Draw the brightest stars
    # first to stop them from completely occluding smaller companions.

    c = catalog
    c = c[c['magnitude'] <= mag1]
    c = c.sort_values('magnitude', ascending=True)
    print('First star group:', len(c))
    s = Star(ra_hours=c.ra_hours, dec_degrees=c.dec_degrees)
    spos = o.observe(s)
    x, y = project(spos)
    scale = 1.5
    radius = (mag1 - c['magnitude']) * scale + 1.0

    x2 = np.repeat(x, 2)
    y2 = np.repeat(y, 2)
    radius2 = (radius[:,None] + (3.0,0.0)).flatten()
    #print(size2[:20])
    c2 = ('w', 'k')

    art.append(ax.scatter(x2, y2, s=radius2 ** 2.0, c=c2))

    return art

@writers.register('ffmpeg_small')
class _FFMpegWriter(FFMpegWriter):
    def _args(self):
        args = super(_FFMpegWriter, self)._args()
        i = args.index('h264') + 1
        args[i:i] = ['-x264-params', 'keyint=99999999:scenecut=0']
        #args.append('keyint=123:min-keyint=20')
        return args

X = []

class _Animation(FuncAnimation):
    patcher = None

    def _draw_next_frame(self, framedata, blit):
        blit = True    # override matplotlib's refusal to blit in save()
        super(_Animation, self)._draw_next_frame(framedata, blit)

        if self.patcher is None:
            self.patcher = patch('matplotlib.backends.backend_agg'
                            '.FigureCanvasAgg.draw', draw)
            self.patcher.__enter__()

def draw(self):
    pass

from unittest.mock import patch

def _animate(projection, t, stars, observer, planet):
    import matplotlib.pyplot as plt

    # print(t)
    # print(t[0].utc_strftime('%Y %B'))

    o = observer.at(t)
    p = o.observe(planet)

    p0 = planet.at(t)
    mag = -8.88 + 5.0 * np.log10(p.distance().au * p0.distance().au)

    project = projection(p)
    x, y = project(p)

    fig, ax = plt.subplots()

    plt.tick_params(axis='both', which='both',
                    bottom=False, labelbottom=False,
                    left=False, labelleft=False)

    fig.set_size_inches(9, 9)

    ax.scatter(x, y, color='b', alpha=0.0)

    _plot_stars(stars, o, project, ax, 6.0, 8.5, 0.8)

    # http://aa.usno.navy.mil/software/mica/USNO-AA-TN-2003-04.pdf
    # m - M = 5 (log10 d - 1)
    # m = 5 (log10 d - 1) + M

    x_left, x_right = ax.get_xlim()
    y_bottom, y_top = ax.get_ylim()
    text_x = x_left
    text_y = y_bottom - (y_top - y_bottom) * 0.01

    # date_text.set_clip_on(False)
    # date_text = fig.text(0.1, 0.1, 'date',
    #                  ha='center',
    #                  va='center',
    #                  transform = ax.transAxes)
    # date_text.set_clip_on(False)
    # date_text = fig.text(0.2, 0.2, 'date',
    #                  ha='center',
    #                  va='center',
    #                  transform = ax.transAxes)

    saturn_color = '#a69276'  # chroma('#d8c2a5').darken().hex()

    planet_art = None
    date_text = None

    def init():
        nonlocal planet_art, date_text
        planet_art, = ax.plot(x[100], y[100], color=saturn_color, marker='o')
        date_text = ax.text(text_x, text_y, '', ha='left', va='top')
        return planet_art, #date_text

    def update(i):
        planet_art.set_xdata(x[i])
        planet_art.set_ydata(y[i])

        scale = 2.0
        radius = (6.0 - mag[i]) * scale

        planet_art.set_ms(radius)

        #date_text = ax.text(text_x, text_y, 'date', ha='left', va='top')
        date_text.set_text(t[i].utc_strftime('%Y %B'))

        if i == 0:
            return planet_art,
        return planet_art, date_text

    frames = len(x)
    frames = 30

    anim = _Animation(fig, update, frames=frames,
                      blit=True, init_func=init,
                      interval=16)
    plt.close()
    rcParams.validate['animation.writer'].valid['ffmpeg_small'] = 'ffmpeg_small'
    with rc_context({'animation.writer': 'ffmpeg_small'}):
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
