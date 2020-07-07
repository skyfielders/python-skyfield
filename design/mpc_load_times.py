#!/usr/bin/env python

from time import time

from skyfield.api import load
from skyfield.data import mpc

def main():
    f = load.open(mpc.COMET_URL)

    t0 = time()
    c = mpc.load_comets_dataframe(f)
    print(time() - t0, 'seconds for load_comets_dataframe()')
    assert len(c) == 864

    f.seek(0)

    t0 = time()
    c = mpc.load_comets_dataframe_slow(f)
    print(time() - t0, 'seconds for load_comets_dataframe_slow()')
    assert len(c) == 864

    # f = load.open(mpc.MPCORB_URL)

    # t0 = time()
    # mpc.load_mpcorb_dataframe(f, slow=False)
    # print(time() - t0, 'seconds for load_mpcorb_dataframe() with slow=False')

    # f.seek(0)

    # t0 = time()
    # mpc.load_mpcorb_dataframe(f, slow=)
    # print(time() - t0, 'seconds for load_mpcorb_dataframe() with slow=')

if __name__ == '__main__':
    main()
