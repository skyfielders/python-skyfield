from skyfield.api import EarthSatellite, load, wgs84

second = 1.0 / 24.0 / 3600.0

def test_satellite_events_on_several_satellites():
    ts = load.timescale()
    boston = wgs84.latlon(42.3581, -71.0636)
    horizon = +20.0

    # Verify 1) rises/sets cross the horizon in the right direction
    # 2) culminations are local maxima above horizon
    # 3) No double rises or double sets

    def run_sat(name, line1, line2, number_events_expected):
        print(name)
        sat = EarthSatellite(line1, line2, name)
        t, y = sat.find_events(topos, t0, t1, horizon)

        assert len(t) == len(y)
        assert len(t) == number_events_expected

        if len(t) == 0:
            return

        geometric = sat - topos
        t3 = ts.tt_jd((t.tt[:,None] + [[-second, 0, +second]]).flatten())
        alt = geometric.at(t3).altaz()[0].degrees
        last_event = None

        for i, (ti, yi) in enumerate(zip(t, y)):
            j = i*3
            alt1, alt2, alt3 = alt[j], alt[j+1], alt[j+2]
            print(i, ti, yi, alt1, alt2, alt3)
            event = ('rise', 'culminate', 'set')[yi]
            if event == 'rise':
                assert alt1 < alt2 < alt3
                assert alt1 < horizon < alt3
                assert last_event in (None, 'set')
            elif event == 'culminate':
                assert alt1 < alt2 > alt3
                assert alt2 > horizon
                assert last_event in (None, 'rise', 'culminate')
            elif event == 'set':
                assert alt1 > alt2 > alt3
                assert alt1 > horizon > alt3
                assert last_event in (None, 'culminate')
            last_event = event

    # Start easy: typical LEO.

    t0 = ts.tai(2014, 11, 10)
    t1 = ts.tai(2014, 11, 11)
    topos = boston

    run_sat(
        'TIANGONG 1',
        '1 37820U 11053A   14314.79851609  .00064249  00000-0  44961-3 0  5637',
        '2 37820  42.7687 147.7173 0010686 283.6368 148.1694 15.73279710179072',
        12,
    )

    # Integral: 3 days high eccentricity.

    t0 = ts.tai(2020, 1, 1)
    t1 = ts.tai(2020, 1, 15)

    run_sat(
        'INTEGRAL',
        '1 27540U 02048A   20007.25125384  .00001047  00000-0  00000+0 0  9992',
        '2 27540  51.8988 127.5680 8897013 285.8757   2.8911  0.37604578 17780',
        36,
    )

    # ANIK-F1R  Geo always visible from Boston

    run_sat(
        'ANIK F-1R',
        '1 28868U 05036A   20011.46493281 -.00000066  00000-0  00000+0 0  9999',
        '2 28868   0.0175  50.4632 0002403 284.1276 195.8977  1.00270824 52609',
        14,
    )

    # PALAPA D  Geo never visible from Boston

    run_sat(
        'PALAPA D',
        '1 35812U 09046A   20008.38785173 -.00000341 +00000-0 +00000-0 0  9999',
        '2 35812 000.0518 095.9882 0002721 218.8296 045.1595 01.00269700038098',
        0,
    )

    # Ariane 5B GTO

    run_sat(
        'Ariane 5B',
        '1 44802U 19080C   20010.68544515  .00001373  00000-0  27860-3 0  9997',
        '2 44802   5.1041 192.7327 7266711 217.6622  57.0965  2.30416801  1028',
        37,
    )

    # Swift Low-inclination LEO never visible from Boston

    run_sat(
        'Swift',
        '1 28485U 04047A   20010.76403232 +.00000826 +00000-0 +25992-4 0  9999',
        '2 28485 020.5579 055.7027 0010957 208.9479 151.0347 15.04516653829549',
        0,
    )

    # Grace-FO 2  Low polar orbit

    run_sat(
        'GRACE-FO 2',
        '1 43477U 18047B   20011.66650462 +.00000719  00000-0 +29559-4 0    08',
        '2 43477  88.9974 159.0391 0019438 141.4770 316.8932 15.23958285 91199',
        90,
    )

    # Issue #996: should detect rising without needing a culmination

    t0 = ts.utc(2024, 8, 26, 8, 38)
    t1 = t0 + 5.0/24.0
    topos = wgs84.latlon(55.671429, 37.62539, 180.0)
    horizon = 14.0
    run_sat(
        'ELEKTRO-L',
        '1 41105U 15074A   24238.84268576 -.00000128  00000+0  00000+0 0  9993',
        '2 41105   5.0153  78.6193 0002491 153.4123  31.4253  1.00270890 31881',
        1,
    )
