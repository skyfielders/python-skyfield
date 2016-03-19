"""Test whether Skyfield handles file download and file age correctly."""

import shutil
import tempfile
from contextlib import contextmanager
from datetime import date
from mock import patch

from skyfield import api

#http://maia.usno.navy.mil/ser7/tai-utc.dat

'''
need files that contain internal state
too great a chance user will move files etc and ruin date

http://maia.usno.navy.mil/ser7/deltat.data - 1 year?
http://maia.usno.navy.mil/ser7/deltat.preds - 1 year?

https://hpiers.obspm.fr/iers/bul/bulc/Leap_Second.dat - "expires" message
'''

old_content = (b' 2015 10  1  67.9546\n'
               b' 2015 11  1  68.0055\n'
               b' 2015 12  1  68.0514\n'
               b' 2016  1  1  68.1024\n')
new_content = (old_content +
               b' 2016  2  1  68.1577\n')

def cache():
    path = tempfile.mkdtemp()
    try:
        yield api.Loader(path)
    finally:
        shutil.rmtree(path)

def save_file(cache, content=old_content):
    with open(cache.path_of('deltat.data'), 'wb') as f:
        f.write(content)

def file_contents(cache):
    with open(cache.path_of('deltat.data'), 'rb') as f:
        return f.read()

@contextmanager
def on(cache, year, month, day):
    fake_date = date(year, month, day)
    download = lambda *args, **kw: save_file(cache, new_content)
    with patch('skyfield.iokit.download', download), \
         patch('skyfield.iokit.today', lambda *args: fake_date):
        yield

# The tests.

def test_missing_file_gets_downloaded(cache):
    with on(cache, 2016, 1, 15):
        data = cache.load('deltat.data')
        assert file_contents(cache).endswith(' 68.1577\n')
    assert data[1][-1] == 68.1577

def test_11_month_old_file_gets_reused(cache):
    save_file(cache)
    with on(cache, 2016, 12, 15):
        data = cache.load('deltat.data')
        assert file_contents(cache).endswith(' 68.1024\n')
    assert data[1][-1] == 68.1024

def test_12_month_old_file_gets_redownloaded(cache):
    save_file(cache)
    with on(cache, 2017, 1, 15):
        data = cache.load('deltat.data')
        assert file_contents(cache).endswith(' 68.1577\n')
    assert data[1][-1] == 68.1577
