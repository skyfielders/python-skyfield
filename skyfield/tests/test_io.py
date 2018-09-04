"""Test whether Skyfield handles file download and file age correctly."""

import os
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

def load():
    path = tempfile.mkdtemp()
    try:
        yield api.Loader(path)
    finally:
        shutil.rmtree(path)

def save_file(load, content=old_content):
    with open(load.path_to('deltat.data'), 'wb') as f:
        f.write(content)

def file_contents(load):
    with open(load.path_to('deltat.data'), 'rb') as f:
        return f.read()

@contextmanager
def on(load, year, month, day):
    fake_date = date(year, month, day)
    download = lambda *args, **kw: save_file(load, new_content)
    with patch('skyfield.iokit.download', download):
        # Python 2.6 does not support the comma "with" statement, so:
        with patch('skyfield.iokit.today', lambda *args: fake_date):
            yield

# The tests.

def test_open_in_main_directory(load):
    with open(os.path.join(load.directory, 'file.tle'), 'wb') as f:
        f.write(b'example text\n')
    data = load.open('file.tle').read()
    print(repr(data))
    assert data == b'example text\n'

def test_open_in_subdirectory(load):
    os.mkdir(os.path.join(load.directory, 'folder'))
    with open(os.path.join(load.directory, 'folder', 'file.tle'), 'wb') as f:
        f.write(b'example text\n')
    data = load.open('folder/file.tle').read()
    assert data == b'example text\n'

def test_missing_file_gets_downloaded(load):
    with on(load, 2016, 1, 15):
        data = load('deltat.data')
        assert file_contents(load).endswith(b' 68.1577\n')
    assert data[1][-1] == 68.1577

def test_11_month_old_file_gets_reused(load):
    save_file(load)
    with on(load, 2016, 12, 15):
        data = load('deltat.data')
        assert file_contents(load).endswith(b' 68.1024\n')
    assert data[1][-1] == 68.1024

def test_12_month_old_file_gets_redownloaded(load):
    save_file(load)
    with on(load, 2017, 1, 15):
        data = load('deltat.data')
        assert file_contents(load).endswith(b' 68.1577\n')
    assert data[1][-1] == 68.1577
