#!/usr/bin/env python2.7

#import argparse
from telnetlib import Telnet

def main(in_path, out_path):
    lines = read_lines(open(in_path))
    tn = Telnet('horizons.jpl.nasa.gov', 6775)
    out = open(out_path, 'w')
    for line in lines:
        print(repr(line))
        tn.write(line.encode('ascii') + b'\r\n')
        data = tn.read_until(b'DUMMY PATTERN', 5.0).decode('ascii')
        print(repr(data))
        out.write(data)
        out.flush()

def read_lines(f):
    for line in f:
        line = line.strip()
        if (not line) or line.startswith('#'):
            continue
        yield line

if __name__ == '__main__':
    try:
        main('horizons_input.txt', 'horizons_output.txt')
    except EOFError:
        print
        print('EOF')
