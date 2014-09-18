#!/usr/bin/env python2.7

#import argparse
from telnetlib import Telnet

def main(path):
    with open(path) as f:
        lines = f.read().split('\n')
    tn = Telnet('horizons.jpl.nasa.gov', 6775)
    out = open('output.txt', 'w')
    for line in lines:
        tn.write(line + '\r\n')
        data = tn.read_until('DUMMY PATTERN', 5.0)
        out.write(data)
        out.flush()

if __name__ == '__main__':
    import sys
    main(sys.argv[1])
