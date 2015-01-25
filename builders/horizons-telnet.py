#!/usr/bin/env python2.7

#import argparse
from telnetlib import Telnet

def main(in_path, out_path):
    with open(in_path) as f:
        lines = f.read().split('\n')
    tn = Telnet('horizons.jpl.nasa.gov', 6775)
    out = open(out_path, 'w')
    for line in lines:
        print(repr(line))
        tn.write(line + '\r\n')
        data = tn.read_until('DUMMY PATTERN', 2.0)
        print(data)
        out.write(data)
        out.flush()

if __name__ == '__main__':
    try:
        main('horizons-input.txt', 'horizons-output.txt')
    except EOFError:
        print
        print('EOF')
