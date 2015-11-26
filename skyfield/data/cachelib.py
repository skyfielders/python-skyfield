"""Experiment with a local file cache; will probably delete soon."""

import os
from datetime import date
from numpy import save
from skyfield import iokit
from skyfield import timelib

_dirname = os.path.dirname(__file__)

cache = iokit.Cache('.')
cache.npy_dirname = _dirname
functions = set([
    timelib.usno_leapseconds,
    ])

def rebuild(remove_old_files=True):
    """Rebuild the data files in this ``skyfield.data`` sub-package.

    You can invoke this routine from the command line with::

        python -m skyfield.data

    """
    from skyfield import iokit
    from skyfield import timelib

    cache = iokit.Cache('.')
    cache.npy_dirname = _dirname
    functions = set([
        timelib.usno_leapseconds,
        ])

    if remove_old_files:
        for filename in os.listdir(_dirname):
            if filename.endswith('.npy'):
                os.unlink(os.path.join(_dirname, filename))

    for function in functions:
        return_value = function(cache)
        filename = function.__name__ + '.npy'
        path = os.path.join(_dirname, filename)
        save(path, return_value)

    with open(os.path.join(_dirname, 'rebuild_date.txt'), 'w') as f:
        f.write('%s\n' % date.today())
