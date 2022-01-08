"""Test whether Skyfield handles file download and file age correctly."""

import os
import shutil
import tempfile
from contextlib import contextmanager
from threading import Thread
try:
    from Queue import Queue
except ImportError:
    from queue import Queue
try:
    from unittest.mock import patch
except ImportError:
    from mock import patch

from assay import assert_raises
from skyfield import api

leap_second_text = (
    b'#  File expires on 28 June 2021\n'
    b'    41317.0    1  1 1972       10\n'
    b'    41499.0    1  7 1972       11\n'
)

def load():
    path = tempfile.mkdtemp()
    try:
        yield api.Loader(path)
    finally:
        shutil.rmtree(path)

def save_file(load, path, content):
    with open(load.path_to(path), 'wb') as f:
        f.write(content)

def file_contents(load, path):
    with open(load.path_to(path), 'rb') as f:
        return f.read()

@contextmanager
def fake_download(load, filename, content):
    download = lambda *args, **kw: save_file(load, filename, content)
    with patch('skyfield.iokit.download', download):
        yield

# Simple tests.

def test_build_url(load):
    url = 'https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/de421.bsp'
    assert load.build_url('de421.bsp') == url
    with assert_raises(ValueError, 'know the URL'):
        load.build_url('unknown.kind.of.file')

def test_open_in_main_directory(load):
    with open(os.path.join(load.directory, 'file.tle'), 'wb') as f:
        f.write(b'example text\n')
    data = load.open('file.tle').read()
    assert data == b'example text\n'

def test_open_in_subdirectory(load):
    os.mkdir(os.path.join(load.directory, 'folder'))
    with open(os.path.join(load.directory, 'folder', 'file.tle'), 'wb') as f:
        f.write(b'example text\n')
    data = load.open('folder/file.tle').read()
    assert data == b'example text\n'

def test_missing_file_gets_downloaded(load):
    with fake_download(load, 'Leap_Second.dat', leap_second_text):
        data = load('Leap_Second.dat')
    assert file_contents(load, 'Leap_Second.dat') == leap_second_text
    assert list(data[1][1]) == [10, 10, 10, 11]

def test_builtin_timescale_uses_recent_IERS_data(load):
    ts = load.timescale()
    # DUT1 cut and pasted from "20 1 1" line of "finals2000A.all":
    assert abs(ts.utc(2020, 1, 1).dut1 - (-0.1771554)) < 1e-8

def test_non_builtin_timescale_prefers_USNO_files(load):
    save_file(load, 'deltat.data',
              b' 1973  2  1  111.1\n'
              b' 1973  3  1  222.2\n'
              b' 1973  4  1  333.3\n')
    save_file(load, 'deltat.preds',
              b'   MJD        YEAR    TT-UT Pred  UT1-UTC Pred  ERROR\n'
              b'   58484.000  2019.00   69.34      -0.152       0.117\n'
              b'   58575.000  2019.25   69.48      -0.295       0.162\n')
    save_file(load, 'Leap_Second.dat', leap_second_text)
    save_file(load, 'finals2000A.all', b'invalid data')

    ts = load.timescale(builtin=False)
    assert abs(ts.utc(1973, 3, 1).delta_t - 222.2) < 1e-2

    # Did we correctly convert the old leap second table to the new
    # format? The offset prior to the first leap second should be +10.
    t = ts.tai(1970, 1, 1)
    assert t.utc == (1969, 12, 31, 23, 59, 50)

def test_non_builtin_timescale_tries_to_load_finals2000A_all(load):
    save_file(load, 'finals2000A.all', b'invalid data')
    with assert_raises(IndexError):
        load.timescale(builtin=False)

def test_close_of_ephemeris(load):
    eph = api.load('de421.bsp')
    assert eph.segments
    assert eph.codes
    eph.close()
    assert eph.segments == []
    assert eph.codes == set()

# Impressive tests: synchronize threads to reproduce concurrency bugs.

class FakeConnection():
    def __init__(self):
        self.control = Queue()
        self.blocks = Queue()
        self.headers = {}

    def __call__(self, *args, **kw):  # Intercept urlopen()
        return self

    def read(self, blocksize):  # Intercept connection.read()
        self.control.put('Read requested')
        return self.blocks.get()

def test_concurrent_downloads(load):
    # While downloads in a single process are protecte

    eof = b''
    data = b' 1973  2  1  43.4724\n'

    c1 = FakeConnection()
    c2 = FakeConnection()
    t1 = Thread(target=load, args=('Leap_Second.dat',))
    t2 = Thread(target=load, args=('Leap_Second.dat',))
    t1.daemon = True
    t2.daemon = True

    with patch('skyfield.iokit.urlopen', c1):
        t1.start()
        c1.control.get()  # has reached its first read()
    with patch('skyfield.iokit.urlopen', c2):
        t2.start()
        c2.control.get()  # has reached its first read()

    assert sorted(os.listdir(load.directory)) == [
        'Leap_Second.dat.download', 'Leap_Second.dat.download2',
    ]
    c1.blocks.put(data)
    c1.blocks.put(eof)
    t1.join()
    assert sorted(os.listdir(load.directory)) == [
        'Leap_Second.dat', 'Leap_Second.dat.download2',
    ]
    c2.blocks.put(data)
    c2.blocks.put(eof)
    t2.join()
    assert sorted(os.listdir(load.directory)) == ['Leap_Second.dat']
