from __future__ import print_function
import os
import numpy as np
import sys
from datetime import datetime, timedelta
from time import time

from .jpllib import SpiceKernel

try:
    from fcntl import LOCK_EX, LOCK_UN, lockf
except:
    lockf = None

try:
    from urllib.request import urlopen, url2pathname
except:
    from urllib2 import urlopen, url2pathname

_missing = object()


def load(filename, autodownload=True, verbose=True):
    """Load the given file, possibly downloading it if it is not present."""
    if filename.endswith('.bsp'):
        url = url_for(filename)
        cls = SpiceKernel
    else:
        raise ValueError('Skyfield does not recognize that file extension')
    if not os.path.exists(filename):
        if not autodownload:
            raise IOError('you specified autodownload=False but {!r} cannot'
                          ' be found in the current directory'
                          .format(filename))
        url = url_for(filename)
        download(url, verbose)
    return cls(filename)


def url_for(filename):
    """Given a recognized filename, return its URL."""
    if filename.endswith('.bsp'):
        if filename.startswith('de'):
            return 'ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/' + filename
        elif filename.startswith('jup'):
            return ('http://naif.jpl.nasa.gov/pub/naif/generic_kernels'
                    '/spk/satellites/' + filename)
        raise ValueError('Skyfield does not know where to download {!r} from'
                         .format(filename))


def download(url, verbose=True, blocksize=128*1024):
    """Download a file from a URL, possibly displaying a progress bar."""
    filename = os.path.basename(url2pathname(url))
    if os.path.exists(filename):
        return filename
    tempname = filename + '.download'
    try:
        connection = urlopen(url)
    except Exception as e:
        raise IOError('cannot fetch file {0!r} from URL {1} because {2}'
                      .format(filename, url, e))
    content_length = int(connection.headers.get('content-length', -1))
    report = ProgressBar(filename).report if verbose else tuple
    with open(tempname, 'ab') as w:
        try:
            if lockf is not None:
                fd = w.fileno()
                lockf(fd, LOCK_EX)
                if os.fstat(fd).st_size:
                    return  # someone else wrote the file contents
            length = 0
            while True:
                data = connection.read(blocksize)
                if not data:
                    break
                w.write(data)
                length += len(data)
                report(length, content_length)
            os.rename(tempname, filename)
        except KeyboardInterrupt:# Exception as e:
            raise IOError('error getting {0} - {1}'.format(url, e))
        finally:
            if lockf is not None:
                lockf(fd, LOCK_UN)
    return filename


class ProgressBar(object):
    def __init__(self, filename):
        self.filename = filename
        self.t0 = 0

    def report(self, bytes_so_far, bytes_total):
        if bytes_total < 0:
            return
        percent = 100 * bytes_so_far // bytes_total
        if (percent != 100) and (time() - self.t0 < 0.5):
            return
        self.t0 = time()
        bar = '#' * (percent // 3)
        print('\r[{0:33}] {1:3}% {2}'.format(bar, percent, self.filename),
              end='\n' if (percent == 100) else '')
        sys.stdout.flush()


class Cache(object):
    """Early experiment in being sensitive to file dates and age."""
    def __init__(self, cache_path, days_old=0):
        self.cache_path = cache_path
        self.days_old = days_old
        self.ram_cache = {}
        self.npy_dirname = None

    def open_url(self, url, days_old=None):
        filename = url[url.rindex('/') + 1:]
        path = os.path.join(self.cache_path, filename)
        if days_old is None:
            days_old = self.days_old
        download_file(url, path, days_old)
        return open(path, 'rb')

    def run(self, function):
        """Return the result of running `function(this_cache)` one time only.

        If this cache has already been asked to run `function`, then the
        return value of its first run is returned without re-running it.

        """
        result = self.ram_cache.get(function, _missing)
        if result is not _missing:
            return result

        if self.npy_dirname:
            path = os.path.join(self.npy_dirname, function.__name__ + '.npy')
            if os.path.exists(path):
                # TODO: check whether data is recent enough
                result = np.load(path)
                self.ram_cache[function] = result
                return result

        result = function(self)
        self.ram_cache[function] = result
        return result


def download_file(url, filename, days_old=0):
    """Early experiment in what download logic might look like."""
    if os.path.exists(filename):
        if not is_days_old(filename, days_old):
            return

    response = urlopen(url)
    f = open(filename, 'wb')
    while True:
        block = response.read(4096)
        if not block:
            break
        f.write(block)

    f.close()

def is_days_old(filename, days_old):
    """Early experiment in being sensitive to file dates."""
    min_old = datetime.now()-timedelta(days=days_old)
    modified = datetime.fromtimestamp(os.path.getmtime(filename))
    return modified < min_old
