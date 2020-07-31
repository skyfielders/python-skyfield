from skyfield import api, almanac

# Compare with USNO:
# http://aa.usno.navy.mil/cgi-bin/aa_moonill2.pl?form=1&year=2018&task=00&tz=-05

def test_fraction_illuminated():
    ts = api.load.timescale()
    t0 = ts.utc(2018, 9, range(9, 19), 5)
    e = api.load('de421.bsp')
    i = almanac.fraction_illuminated(e, 'moon', t0[-1]).round(2)
    assert i == 0.62
    i = almanac.fraction_illuminated(e, 'moon', t0).round(2)
    assert (i == (0, 0, 0.03, 0.08, 0.15, 0.24, 0.33, 0.43, 0.52, 0.62)).all()

# Compare with USNO:
# http://aa.usno.navy.mil/seasons?year=2018&tz=+0

def test_seasons():
    ts = api.load.timescale()
    t0 = ts.utc(2018, 9, 20)
    t1 = ts.utc(2018, 12, 31)
    e = api.load('de421.bsp')
    t, y = almanac.find_discrete(t0, t1, almanac.seasons(e))
    strings = t.utc_strftime('%Y-%m-%d %H:%M')
    print(strings)
    assert strings == ['2018-09-23 01:54', '2018-12-21 22:23']
    assert (y == (2, 3)).all()

# Compare with USNO:
# http://aa.usno.navy.mil/cgi-bin/aa_phases.pl?year=2018&month=9&day=11&nump=50&format=p

def test_moon_phases():
    ts = api.load.timescale()
    t0 = ts.utc(2018, 9, 11)
    t1 = ts.utc(2018, 9, 30)
    e = api.load('de421.bsp')
    t, y = almanac.find_discrete(t0, t1, almanac.moon_phases(e))
    strings = t.utc_strftime('%Y-%m-%d %H:%M')
    assert strings == ['2018-09-16 23:15', '2018-09-25 02:52']
    assert (y == (1, 2)).all()

def test_oppositions_conjunctions():
    ts = api.load.timescale()
    t0 = ts.utc(2019, 1, 1)
    t1 = ts.utc(2021, 1, 1)
    e = api.load('de421.bsp')
    f = almanac.oppositions_conjunctions(e, e['mars'])
    t, y = almanac.find_discrete(t0, t1, f)
    strings = t.utc_strftime('%Y-%m-%d %H:%M')
    assert strings == ['2019-09-02 10:42', '2020-10-13 23:26']
    assert (y == (0, 1)).all()

# Compare with USNO:
# http://aa.usno.navy.mil/rstt/onedaytable?ID=AA&year=2018&month=9&day=12&state=OH&place=Bluffton

def test_sunrise_sunset():
    ts = api.load.timescale()
    t0 = ts.utc(2018, 9, 12, 4)
    t1 = ts.utc(2018, 9, 13, 4)
    e = api.load('de421.bsp')
    bluffton = api.Topos('40.8939 N', '83.8917 W')
    t, y = almanac.find_discrete(t0, t1, almanac.sunrise_sunset(e, bluffton))
    strings = t.utc_strftime('%Y-%m-%d %H:%M')
    assert strings == ['2018-09-12 11:13', '2018-09-12 23:50']
    assert (y == (1, 0)).all()

def test_dark_twilight_day():
    ts = api.load.timescale()
    t0 = ts.utc(2019, 11, 8, 4)
    t1 = ts.utc(2019, 11, 9, 4)
    e = api.load('de421.bsp')
    defiance = api.Topos('41.281944 N', '84.362778 W')
    t, y = almanac.find_discrete(t0, t1, almanac.dark_twilight_day(e, defiance))
    strings = t.utc_strftime('%Y-%m-%d %H:%M')
    assert strings == [
        '2019-11-08 10:42', '2019-11-08 11:15', '2019-11-08 11:48',
        '2019-11-08 12:17', '2019-11-08 22:25', '2019-11-08 22:54',
        '2019-11-08 23:27', '2019-11-08 23:59',
    ]
    assert (y == (1, 2, 3, 4, 3, 2, 1, 0)).all()

# Logic.

def test_close_start_and_end():
    ts = api.load.timescale()
    t0 = ts.utc(2018, 9, 23, 1)
    t1 = ts.utc(2018, 9, 23, 2)
    e = api.load('de421.bsp')
    t, y = almanac.find_discrete(t0, t1, almanac.seasons(e))
    strings = t.utc_strftime('%Y-%m-%d %H:%M')
    assert strings == ['2018-09-23 01:54']
