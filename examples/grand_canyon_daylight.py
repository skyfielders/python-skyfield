from skyfield.api import load, wgs84
from skyfield import almanac

ts = load.timescale()
eph = load('de421.bsp')
grand_canyon_village = wgs84.latlon(+36.0544, -112.1401)
tz_hours = -7
start = ts.utc(1986, 1, 17, -tz_hours)
#end = ts.utc(1987, 7, 14, -tz_hours)
end = ts.utc(2022, 7, 14, -tz_hours)
ss = almanac.sunrise_sunset(eph, grand_canyon_village)
t, y = almanac.find_discrete(start, end, ss)
sunrises = t[::2]
sunsets = t[1::2]
lengths = sunsets - sunrises
print('date,hours_sunlight')
for ti, length in zip(sunrises, lengths):
    print('{},{}'.format(ti.utc_strftime('%Y-%m-%d'), length * 24.0))
