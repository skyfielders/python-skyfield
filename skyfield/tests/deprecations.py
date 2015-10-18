"""Stand-alone Python script to test deprecated APIs."""

import sys

def main():
    args = sys.argv[1:]
    if not args:
        print('usage: deprecations.py (a|b...)')
        sys.exit(2)
    arg = args[0]
    if arg == 'a':
        from skyfield.api import earth, mars, now
        earth(now()).observe(mars).radec()
    elif arg == 'b':
        from skyfield.api import earth
        earth.topos('42.3583 N', '71.0636 W')
    elif arg == 'c':
        from skyfield.api import load
        eph = load('de421.bsp')
        earth = eph['earth']
        earth(100)

if __name__ == '__main__':
    main()
