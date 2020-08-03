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

from skyfield import api

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
def fake_download(load):
    download = lambda *args, **kw: save_file(load, new_content)
    with patch('skyfield.iokit.download', download):
        yield

# Simple tests.

def test_build_url(load):
    url = 'ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/de421.bsp'
    assert load.build_url('de421.bsp') == url
    assert load.build_url('unknown.kind.of.file') is None

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
    with fake_download(load):
        data = load('deltat.data')
    print(repr(file_contents(load)[:-20]))
    assert file_contents(load).endswith(b' 68.1577\n')
    assert data[1][-1] == 68.1577

def test_builtin_timescale(load):
    ts = load.timescale()
    ts.utc(2019, 7, 21, 11, 11)

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
