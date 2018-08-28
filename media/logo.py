# -*- coding: utf-8 -*-

"""Build the Skyfield logo."""

# To build the logo and make a PNG:
#
# wget http://tdc-www.harvard.edu/catalogs/bsc5.dat.gz
# gunzip bsc5.dat.gz
#
# wget http://www.impallari.com/media/releases/dosis-v1.7.zip
# [extract Dosis-Medium.ttf]
#
# python logo.py
#
# convert -density 480 logo.pdf logo.png

from math import cos, pi, sin
try:
    intern         # Python 2
except NameError:  # Python 3
    from sys import intern
try:
    from reportlab.pdfgen.canvas import Canvas
    from reportlab.pdfbase import pdfmetrics
    from reportlab.pdfbase.ttfonts import TTFont
except ImportError:
    pass  # for when py.test discovers us but ReportLab is not installed

tau = pi * 2.0
quarter_tau = pi / 2.0

alpha_star_color = '#BC243C'
bright_star_color = '#002664'
dim_star_color = '#aaaaaa'
text_color = '#002664'
greek_color = 'black'

greeks = {
    'Alp': (u'α', 1.0),
    'Bet': (u'β', 1.0),
    'Gam': (u'γ', 1.3),
    'Del': (u'δ', 0.9),
    'Eps': (u'ε', 1.0),
    'Zet': (u'ζ', 1.0),
    'Eta': (u'η', 1.0),
    }

def main():
    pdfmetrics.registerFont(TTFont('Dosis', 'Dosis-Medium.ttf'))

    stars = []
    with open('bsc5.dat', 'rb') as f:
        for line in f:
            line = '.' + line  # switch to 1-based indexing
            if not line[62].strip():
                continue  # skip coordinate-less novas
            letter = intern(line[8:11])
            if letter == '   ':
                letter = None
            h, m, s = float(line[61:63]), float(line[63:65]), float(line[65:69])
            ra = (h + (m + s / 60.0) / 60.0) * tau / 24.0
            d, m, s = float(line[69:72]), float(line[72:74]), float(line[76:78])
            dec = (d + (m + s / 60.0) / 60.0) * tau / 360.0
            mag = float(line[103:108])
            stars.append((letter, ra, dec, mag))

    h, w = 48, 96
    c = Canvas('logo.pdf', pagesize=(w, h))

    c.setFillColor('white')
    c.rect(0, 0, w, h, stroke=0, fill=1)

    c.setFillColor(bright_star_color)

    rotation = 10.0 * tau / 360.0
    # magscale = 0.1
    # For 10 degrees:
    x_offset = 96 -33.5
    y_offset = h  +37.5
    # For 15 degrees:
    # x_offset = 96 -28.5
    # y_offset = 96 +0.5
    # for 45 degrees:
    # x_offset = 96 -13.5
    # y_offset = 96 -10

    small_glyphs = []
    c.setFont('Helvetica', 2)

    for letter, ra, dec, mag in stars:
        # if mag > 4.0:
        #     continue
        d = - (dec - quarter_tau) * 100
        ra += rotation
        x = d * sin(ra)
        y = d * cos(ra)

        if y < -63.0 or y > -39.0:
            continue
        if x < -43.0 or x > 19.0:
            continue

        x += x_offset
        y += y_offset

        r = ((13.0 - mag) / 10.0) ** 4.0 #* magscale
        r = min(r, 1.0)

        if r < 0.5:
            small_glyphs.append((x, y, r))
        else:
            if letter is not None:
                c.saveState()
                greek_letter, offset = greeks[letter]
                c.setFillColor(greek_color)
                c.drawString(x+offset, y+0.5, greek_letter)
                if letter == 'Alp':
                    c.setFillColor(alpha_star_color)
                    c.circle(x, y, r, stroke=0, fill=1)
                c.restoreState()
                if letter != 'Alp':
                    c.circle(x, y, r, stroke=0, fill=1)
            else:
                c.circle(x, y, r, stroke=0, fill=1)


    c.setFillColor(dim_star_color)
    for x, y, r in small_glyphs:
        c.circle(x, y, r, stroke=0, fill=1)

    c.setFillColor(text_color) #, alpha=0.5)
    c.setFont('Dosis', 24)
    sw = c.stringWidth('Skyfield')
    c.drawString(w // 2 - sw // 2, h - 40, 'Skyfield')

    c.showPage()
    with open('logo.pdf', 'wb') as f:
        f.write(c.getpdfdata())

if __name__ == '__main__':
    main()
