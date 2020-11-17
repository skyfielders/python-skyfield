from __future__ import print_function, division
from skyfield import api

def test_sat_almanac_LEO():
    # Testcase from
    # https://github.com/skyfielders/astronomy-notebooks/blob/master/Solvers/Earth-Satellite-Passes.ipynb
    # with elevated horizon

    tle = ["TIANGONG 1",
           "1 37820U 11053A   14314.79851609  .00064249  00000-0  44961-3 0  5637",
           "2 37820  42.7687 147.7173 0010686 283.6368 148.1694 15.73279710179072"]
    sat = api.EarthSatellite(*tle[1:3], name=tle[0])
    topos = api.Topos('42.3581 N', '71.0636 W')
    timescale = api.load.timescale()
    t0 = timescale.tai(2014, 11, 10)
    t1 = timescale.tai(2014, 11, 11)
    horizon = 20
    nexpected = 12

    times, yis = sat.find_events(topos, t0, t1, 20.0)
    assert(verify_sat_almanac(times, yis, sat, topos, horizon, nexpected))


def test_sat_almanac_tricky():
    # Various tricky satellites
    # Integral: 3 days high eccentricity.
    # ANIK-F1R  Geo always visible from Boston
    # PALAPA D  Geo never visible from Boston
    # Ariane 5B GTO
    # Swift Low-inclination LEO never visible from Boston
    # Grace-FO 2  Low polar orbit
    tles = """\
        INTEGRAL
        1 27540U 02048A   20007.25125384  .00001047  00000-0  00000+0 0  9992
        2 27540  51.8988 127.5680 8897013 285.8757   2.8911  0.37604578 17780
        ANIK F-1R
        1 28868U 05036A   20011.46493281 -.00000066  00000-0  00000+0 0  9999
        2 28868   0.0175  50.4632 0002403 284.1276 195.8977  1.00270824 52609
        PALAPA D
        1 35812U 09046A   20008.38785173 -.00000341 +00000-0 +00000-0 0  9999
        2 35812 000.0518 095.9882 0002721 218.8296 045.1595 01.00269700038098
        Ariane 5B
        1 44802U 19080C   20010.68544515  .00001373  00000-0  27860-3 0  9997
        2 44802   5.1041 192.7327 7266711 217.6622  57.0965  2.30416801  1028
        Swift
        1 28485U 04047A   20010.76403232 +.00000826 +00000-0 +25992-4 0  9999
        2 28485 020.5579 055.7027 0010957 208.9479 151.0347 15.04516653829549
        GRACE-FO 2
        1 43477U 18047B   20011.66650462 +.00000719  00000-0 +29559-4 0    08
        2 43477  88.9974 159.0391 0019438 141.4770 316.8932 15.23958285 91199
    """.splitlines()
    expected = {'INTEGRAL': 36,
                'ANIK F-1R': 14,
                'PALAPA D': 0,
                'Ariane 5B': 37,
                'Swift': 0,
                'GRACE-FO 2': 90}
    for iline0 in range(0, len(tles)-3, 3):
        #print('---')
        sat = api.EarthSatellite(tles[1+iline0].strip(),
                                 tles[2+iline0].strip(),
                                 name=tles[iline0].strip())
        topos = api.Topos('42.3581 N', '71.0636 W')  # Boston
        timescale = api.load.timescale()
        t0 = timescale.tai(2020, 1, 1)
        t1 = timescale.tai(2020, 1, 15)
        horizon = 20
        nexpected = expected[sat.name]

        times, yis = sat.find_events(topos, t0, t1, 20.0)
        assert(verify_sat_almanac(times, yis, sat, topos, horizon, nexpected))



# Helper function to verify satellite events
def verify_sat_almanac(times, yis, sat, topos, horizon, nexpected):
    if nexpected is None:
        print("Number of satellite events expected for", sat.name," not specified.  Got ", len(times))
    else:
        if len(times) != nexpected or len(yis) != nexpected:
            raise RuntimeError("Expected {} events for satellite {}"
                               " but got {} times and {} y values"
                               .format(nexpected, sat.model.satnum,
                                       len(times), len(yis)))
    if len(times) == 0:
        return True     # Nothing to check
    # Verify 1) rises/sets cross the horizon in the right direction
    # 2) culminations are local maxima above horizon
    # 3) No double rises or double sets
    time_tolerance = 5/(24 * 60 * 60)   # Times must be accurate within 5 seconds
    ts = times[0].ts     # Use the timescale passed by
    altitudes_neartimes = [
        (sat - topos).at(ts.tai(jd=times.tai + dt)).altaz()[0].degrees
        for dt in (-time_tolerance, 0, time_tolerance)
    ]
    lastevent = None
    event_dict = {0: 'rise', 1: 'culminate', 2: 'set'}
    for time, yi, (alt_before, alt_at, alt_after) in zip(
            times, yis, zip(*altitudes_neartimes)):
        eventname = event_dict[yi]
        if eventname == 'rise':
            assert(alt_before < alt_at < alt_after)  # It is going up
            assert(alt_before < horizon < alt_after) # It is crossing horizon altitude
            assert(lastevent != 'rise')     # Can't have two consecutive rises
        elif eventname == 'culminate':
            assert(alt_before < alt_at)
            assert(alt_at > alt_after)
            assert(alt_at > horizon)
            # You can have two consecutive culminations (esp. for GEOs).
        elif eventname == 'set':
            assert(alt_before > alt_at > alt_after)  # It is going down
            assert(alt_before > horizon > alt_after) # It is crossing horizon altitude
            assert(lastevent != 'set')     # Can't have two consecutive sets
        else:
            raise RuntimeError("Unexpected satellite event type")
        lastevent = eventname
    return True
