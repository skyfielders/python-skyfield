"""Build the Skyfield logo."""

# To make a PNG:
# convert -density 480 logo.pdf logo.png

from math import cos, pi, sin
from reportlab.pdfgen.canvas import Canvas
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont

tau = pi * 2.0
quarter_tau = pi / 2.0

def main():
    pdfmetrics.registerFont(TTFont('Dosis', 'Dosis.ttf'))

    stars = []
    with open('bsc5.dat', 'rb') as f:
        for line in f:
            line = '.' + line  # switch to 1-based indexing
            # print line
            if not line[62].strip():
                continue  # skip coordinate-less novas
            h, m, s = float(line[61:63]), float(line[63:65]), float(line[65:69])
            ra = (h + (m + s / 60.0) / 60.0) * tau / 24.0
            d, m, s = float(line[69:72]), float(line[72:74]), float(line[76:78])
            dec = (d + (m + s / 60.0) / 60.0) * tau / 360.0
            mag = float(line[103:108])
            stars.append((ra, dec, mag))

    h, w = 72, 96
    c = Canvas('logo.pdf', pagesize=(w, h))

    c.setFillColor('white')
    c.rect(0, 0, w, h, stroke=0, fill=1)

    c.setFillColor('black')
    c.setFont('Dosis', 24)
    sw = c.stringWidth('Skyfield')
    c.drawString(w // 2 - sw // 2, 20, 'Skyfield')

    rotation = 10.0 * tau / 360.0
    # magscale = 0.1
    # For 10 degrees:
    x_offset = 96 -33.76
    y_offset = 96 +1.63
    # For 15 degrees:
    # x_offset = 96 -28.5
    # y_offset = 96 +0.5
    # for 45 degrees:
    # x_offset = 96 -13.5
    # y_offset = 96 -10

    for ra, dec, mag in stars:
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

        r = ((13.0 - mag) / 10.0) ** 4.0 #* magscale
        r = min(r, 1.0)
        #print r
        c.circle(x_offset + x, y_offset + y, r, stroke=0, fill=1)
    c.showPage()
    with open('logo.pdf', 'wb') as f:
        f.write(c.getpdfdata())

if __name__ == '__main__':
    main()
