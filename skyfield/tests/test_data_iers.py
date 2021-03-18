from numpy import array, inf
from skyfield.data.iers import _build_timescale_arrays
from skyfield.timelib import Timescale

def test_build_timescale_arrays():
    mjd = array([42046.00, 42047.00, 42048.00, 42049.00])
    dut1 = array([-0.2942581, -0.2971424, 0.6999438, 0.6970539])
    delta_t, leap_dates, leap_offsets = _build_timescale_arrays(mjd, dut1)
    column_1, column_2 = delta_t

    assert list(column_1) == [2442046.5005113888, 2442047.5005113888,
                              2442048.500522963, 2442049.500522963]
    assert list(column_2) == [44.4782581, 44.481142399999996,
                              44.4840562, 44.4869461]
    assert list(leap_dates) == [-inf, 2441317.5, 2441499.5, 2441683.5,
                                2442048.5, inf]
    assert list(leap_offsets) == [10, 10, 10, 11, 12, 13]

    ts = Timescale(delta_t, leap_dates, leap_offsets)
    t = ts.tai(1973, 1, 1, 0, 0, [11, 12])
    assert t.utc_strftime() == ['1972-12-31 23:59:60 UTC',
                                '1973-01-01 00:00:00 UTC']
