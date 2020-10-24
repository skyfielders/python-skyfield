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

data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '../data'))
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

def save_file(load, path, content=old_content):
    with open(load.path_to(path), 'wb') as f:
        f.write(content)

def file_contents(load, path):
    with open(load.path_to(path), 'rb') as f:
        return f.read()

@contextmanager
def fake_download(load):
    download = lambda *args, **kw: save_file(load, 'deltat.data', new_content)
    with patch('skyfield.iokit.download', download):
        yield

# Simple tests.

def test_build_url(load):
    url = 'ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/de421.bsp'
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
    with fake_download(load):
        data = load('deltat.data')
    assert file_contents(load, 'deltat.data').endswith(b' 68.1577\n')
    assert data[1][-1] == 68.1577

def test_builtin_timescale_uses_recent_IERS_data(load):
    ts = load.timescale()
    # DUT1 cut and pasted from "20 1 1" line of "finals2000A.all":
    assert abs(ts.utc(2020, 1, 1).dut1 - (-0.1771547)) < 1e-8

def test_non_builtin_timescale_prefers_USNO_files(load):
    with open(os.path.join(data_dir, 'deltat.preds'), 'rb') as f:
        preds = f.read()
    with open(os.path.join(data_dir, 'Leap_Second.dat'), 'rb') as f:
        leaps = f.read()

    save_file(load, 'deltat.data',
              b' 1973  2  1  111.1\n'
              b' 1973  3  1  222.2\n'
              b' 1973  4  1  333.3\n')
    save_file(load, 'deltat.preds', preds)
    save_file(load, 'Leap_Second.dat', leaps)
    save_file(load, 'finals2000A.all', b'invalid data')

    ts = load.timescale(builtin=False)
    assert abs(ts.utc(1973, 3, 1).delta_t - 222.2) < 1e-2

def test_non_builtin_timescale_tries_to_load_finals2000A_all(load):
    save_file(load, 'finals2000A.all', b'invalid data')
    with assert_raises(IndexError):
        load.timescale(builtin=False)

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
    t1 = Thread(target=load, args=('deltat.data',))
    t2 = Thread(target=load, args=('deltat.data',))
    t1.daemon = True
    t2.daemon = True

    with patch('skyfield.iokit.urlopen', c1):
        t1.start()
        c1.control.get()  # has reached its first read()
    with patch('skyfield.iokit.urlopen', c2):
        t2.start()
        c2.control.get()  # has reached its first read()

    assert sorted(os.listdir(load.directory)) == [
        'deltat.data.download', 'deltat.data.download2',
    ]
    c1.blocks.put(data)
    c1.blocks.put(eof)
    t1.join()
    assert sorted(os.listdir(load.directory)) == [
        'deltat.data', 'deltat.data.download2',
    ]
    c2.blocks.put(data)
    c2.blocks.put(eof)
    t2.join()
    assert sorted(os.listdir(load.directory)) == ['deltat.data']
