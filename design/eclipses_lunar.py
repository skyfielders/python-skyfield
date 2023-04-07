"""Check Skyfield lunar eclipses against the huge NASA table of eclipses."""

from collections import Counter
from skyfield import eclipselib
from skyfield.api import load, GREGORIAN_START
import datetime as dt

with load.open('https://eclipse.gsfc.nasa.gov/5MCLE/5MKLEcatalog.txt') as f:
    data = f.read()

ts = load.timescale()
ts.julian_calendar_cutoff = GREGORIAN_START

table = []
lines = data.decode('ascii').splitlines()

for line in lines[14:]:
    year = int(line[7:12])
    #if year < 1 or year >= 3000:
    #if year < 1000 or year > 2500:  # Computes the claim made in `almanac.rst`.
    #if year < 1900 or year > 1920:
    if year >= 3000:  # Minimum necessary to work: DE422 stops at year 3000.
       continue
    datestr = line[13:29]
    d = dt.datetime.strptime('2024 ' + datestr, '%Y %b %d  %H:%M:%S')
    t = ts.tt(year, d.month, d.day, d.hour, d.minute, d.second)
    table.append((t, line[51]))

time0 = table[0][0]
timeN = table[-1][0]

start_time = ts.tt_jd(time0.whole - 1, time0.tt_fraction)
end_time = ts.tt_jd(timeN.whole + 1, timeN.tt_fraction)

print('Start time:', start_time.utc_iso(' '))
print('End time:  ', end_time.utc_iso(' '))

start_time.whole -= 1.0
end_time.whole += 1.0

eph = load('de422.bsp')
#print(eph)  # To double-check which dates it covers, in case of doubt.
t, y, details = eclipselib.lunar_eclipses(start_time, end_time, eph)

max_diff = 0.0
total_diff = 0.0
problems = Counter()

i = 0

for ti, yi in zip(t, y):
    letter2 = 'NPT'[yi]
    t5, letter = table[i]
    diff = (ti - t5) * 24.0 * 3600.0
    #diff = (ti.utc_datetime() - t5.utc_datetime()).total_seconds()
    if diff < -100.0:
        problems['Extra'] += 1
        print(ti.utc_iso(), yi, 'EXTRA:not in catalog', letter2)
        continue
    while diff > 100.0:
        problems['Missing'] += 1
        print(' ' * len(t5.utc_iso()), '-', t5.utc_iso(), letter,
              'MISSING eclipse!')
        i += 1
        t5, letter = table[i]
        diff = (ti - t5) * 24.0 * 3600.0
    max_diff = max(max_diff, abs(diff))
    total_diff += abs(diff)
    judge = '' if letter == letter2 else f'their {letter} != our {letter2}'
    if judge:
        problems['Mismatch'] += 1
    print(ti.utc_iso(), yi, t5.utc_iso(), letter, diff, judge)
    i += 1

print('Largest difference in time of an eclipse (seconds):', max_diff)
print('Total difference in seconds between eclipse times:', total_diff)
for name, count in sorted(problems.items()):
    print(f'{name}: {count} / {len(table)}')
