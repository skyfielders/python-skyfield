"""Routines to download Earth orientation data."""

import numpy as np
from skyfield.iokit import load
from skyfield.timelib import (delta_t_formula_morrison_and_stephenson_2004,
                              julian_date)

def morrison_and_stephenson_2004_table():
    """Table of smoothed Delta T values from Morrison and Stephenson, 2004."""
    import pandas as pd
    f = load('http://eclipse.gsfc.nasa.gov/SEcat5/deltat.html')
    tables = pd.read_html(f.read())
    df = tables[0]
    return pd.DataFrame({'year': df[0], 'delta_t': df[1]})

def usno_historic_delta_t():
    import pandas as pd
    f = load('http://maia.usno.navy.mil/ser7/historic_deltat.data')
    df = pd.read_table(f, sep=r' +', engine='python', skiprows=[1])
    return pd.DataFrame({'year': df['Year'], 'delta_t': df['TDT-UT1']})

def usno_monthly_delta_t():
    import pandas as pd
    f = load('http://maia.usno.navy.mil/ser7/deltat.data')
    return pd.read_table(f, sep=r' +', engine='python',
                         names=['year', 'month', 'day', 'delta_t'])

def usno_predicted_delta_t():
    import pandas as pd
    f = load('http://maia.usno.navy.mil/ser7/deltat.preds')
    df = pd.read_table(f, sep=r'  +', engine='python')
    return pd.DataFrame({'year': df['YEAR'],
                         'delta_t': df['TT-UT PREDICTION']})

def build_delta_t_table():
    # TODO: julian_date() is proleptic Gregorian - but is that the
    # calendar used by Morrison and Stephenson?

    s = morrison_and_stephenson_2004_table()
    s['tt'] = julian_date(s.pop('year'))

    h = usno_historic_delta_t()
    h['tt'] = julian_date(h.pop('year'))

    m = usno_monthly_delta_t()
    m['tt'] = julian_date(m.pop('year'), m.pop('month'), m.pop('day'))

    p = usno_predicted_delta_t()
    p['tt'] = julian_date(p.pop('year'))

    s = s[s.tt < h.tt.min()]
    h = h[h.tt < m.tt.min()]
    p = p[p.tt > m.tt.max()]

    # Generate an initial and final data point from the long-term
    # formula to avoid a discontinuity between numbers we interpolate
    # from the table and numbers we generate from the formula.

    f = delta_t_formula_morrison_and_stephenson_2004

    step = 36525.0  # Julian century

    import pandas as pd

    tt = s.tt.iloc[0] - step
    start = pd.DataFrame({'tt': [tt], 'delta_t': [f(tt)]})

    tt = p.tt.iloc[-1] + step
    end = pd.DataFrame({'tt': [tt], 'delta_t': [f(tt)]})

    return pd.concat([start, s, h, m, p, end])

def delta_t():
    table = build_delta_t_table()
    return np.array([table.tt.values, table.delta_t.values])
