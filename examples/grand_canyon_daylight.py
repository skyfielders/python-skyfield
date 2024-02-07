# Runtime was around 1 minute 4 seconds under the old find-discrete
# regime; with the new find_risings() and find_settings() routines, this
# script now takes around 5 seconds.  Nice.

from skyfield.api import load, wgs84
from skyfield import almanac

ts = load.timescale()
eph = load('de421.bsp')
grand_canyon_village = wgs84.latlon(+36.0544, -112.1401)
observer = eph['earth'] + grand_canyon_village
tz_hours = -7
start = ts.utc(1986, 1, 17, -tz_hours)
#end = ts.utc(1987, 7, 14, -tz_hours)
end = ts.utc(2022, 7, 14, -tz_hours)
sunrises, _ = almanac.find_risings(observer, eph['Sun'], start, end)
sunsets, _ = almanac.find_settings(observer, eph['Sun'], start, end)
lengths = sunsets - sunrises
print('date,hours_sunlight')
for ti, length in zip(sunrises, lengths):
    print('{},{}'.format(ti.utc_strftime('%Y-%m-%d'), length * 24.0))
