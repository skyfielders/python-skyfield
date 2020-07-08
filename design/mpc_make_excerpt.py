# mpc_make_excerpt.py

"""Search the MPCORB file for minor planets, given their packed designations."""

import argparse
import re
import sys
import zlib

from skyfield.api import load
from skyfield.data import mpc

def main(argv):
    parser = argparse.ArgumentParser(description='Grep MPCORB.DAT.gz')
    parser.add_argument('designations', nargs='+', help='packed designations'
                        ' of the minor planets whose orbits you need')
    args = parser.parse_args(argv)

    designations = [re.escape(d.encode('ascii')) for d in args.designations]
    pattern = rb'^((?:%s) .*\n)' % rb'|'.join(designations)
    r = re.compile(pattern, re.M)

    data = load.open(mpc.MPCORB_URL).read()
    data = zlib.decompress(data, wbits = zlib.MAX_WBITS | 16)
    lines = r.findall(data)

    sys.stdout.buffer.write(b''.join(lines))

if __name__ == '__main__':
    main(sys.argv[1:])
