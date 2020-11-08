"""Check Skyfield lunar eclipses against the huge NASA table of eclipses."""

from skyfield import almanac
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
    if year < 1 or year > 1000:
    #if year < 1900 or year > 1920:
        continue
    datestr = line[8:29]
    d = dt.datetime.strptime(datestr, '%Y %b %d  %H:%M:%S')
    t = ts.tt(d.year, d.month, d.day, d.hour, d.minute, d.second)
    table.append((t, line[51]))

time0 = table[0][0]
timeN = table[-1][0]

start_time = ts.tt_jd(time0.whole - 1, time0.tt_fraction)
end_time = ts.tt_jd(timeN.whole + 1, timeN.tt_fraction)

start_time.whole -= 1.0
end_time.whole += 1.0

#eph = load('de421.bsp')
eph = load('de406.bsp')
t, y = almanac.lunar_eclipses(eph, start_time, end_time)

max_diff = 0.0
total_diff = 0.0
total_judge = 0

i = 0

for ti, yi in zip(t, y):
    t5, letter = table[i]
    diff = (ti.utc_datetime() - t5.utc_datetime()).total_seconds()
    if diff < -100.0:
        print('SKIPPING eclipse we found but is not in their list')
        continue
    i += 1
    while diff > 100.0:
        print('MISSING eclipse!')
        t5, letter = table[i]
        i += 1
        diff = (ti.utc_datetime() - t5.utc_datetime()).total_seconds()
    max_diff = max(max_diff, abs(diff))
    total_diff += abs(diff)
    letter2 = 'NPT'[yi]
    judge = '' if letter == letter2 else f'their {letter} != our {letter2}'
    if judge:
        total_judge += 1
    print(ti.tt_strftime(), yi, t5.tt_strftime(), letter, diff, judge)

print('Largest difference in time of an eclipse (seconds):', max_diff)
print('Total difference in seconds between eclipse times:', total_diff)
print('Number of eclipses with different descriptions:',
      total_judge, '/', len(table))
