#!/usr/bin/python2.7
#
# Adapted from:
# http://stackoverflow.com/questions/14557944/downsizing-an-otf-font-by-removing-glyphs

import sys
try:
    import fontforge
except ImportError:
    print('ImportError, try: sudo apt install python-fontforge')
    exit(2)

if len(sys.argv) >= 4:
    font = fontforge.open(sys.argv[1])

    for path in sys.argv[3:]:
        with open(path, "rb") as f:
            data = f.read()
        for i in data.decode("UTF-8"):
            font.selection[ord(i)] = True
    f.close()

    font.selection.invert()

    for i in font.selection.byGlyphs:
        font.removeGlyph(i)

    font.generate(sys.argv[2])
else:
    print("Usage: {} source.woff output.woff source1.txt ...".format(sys.argv[0]))
    print("Example: {} /path/to/ukai.ttc chineseTranslation.txt ukaiStripped.ttf".format(sys.argv[0]))
    exit(2)
