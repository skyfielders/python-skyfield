#!/usr/bin/env python3

import numpy as np
from skyfield.api import Time, load

ATTRIBUTES = (
    'J', 'delta_t', 'dut1', 'gmst',
    # (lambda t: t.toordinal()),
    # (lambda t: t._utc_float(0.0)),
    'tai_fraction', 'tdb_fraction', 'ut1_fraction',
)

def main():
    ts = load.timescale()
    t = ts.utc(2020, 10, 24)
    for attribute in ATTRIBUTES:
        step_width = measure_step_width(ts, t.whole, attribute)
        step_width_s = step_width * 24 * 60 * 60
        print('{:20}  {:.4g} seconds'.format(attribute, step_width_s))

def measure_step_width(ts, whole, attribute):
    samples = np.arange(-500, 501)
    widths = []

    for tenth in 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9:
        exp = -30
        while True:
            tt_fraction = tenth + samples * 10 ** exp
            t = Time(ts, whole, tt_fraction)
            output = getattr(t, attribute)
            diff = np.diff(output)
            steps_up = (diff > 0.0).sum()
            steps_down = (diff < 0.0).sum()
            steps = steps_up + steps_down
            if steps >= 10:
                break
            exp += 1
        step_width = (tt_fraction[-1] - tt_fraction[0]) / steps
        widths.append(step_width)

    return max(widths)

if __name__ == '__main__':
    main()
