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
    from urllib.request import urlopen
except:
    from urllib2 import urlopen

_missing = object()


def load(filename, directory='.', autodownload=True, verbose=True):
    """Load the given file, possibly downloading it if it is not present."""
    if '/' in filename:
        url = filename
        filename = url.split('/')[-1]
        cls = None
    elif filename.endswith('.bsp'):
        url = url_for(filename)
        cls = SpiceKernel
    else:
        raise ValueError('Skyfield does not recognize that file extension')
    if directory == '.':
        path = filename
    else:
        directory = os.path.expanduser(directory)
        path = os.path.join(directory, filename)
    if not os.path.exists(path):
        if not autodownload:
            raise IOError('you specified autodownload=False but the file'
                          ' does not exist: {}'.format(path))
        download(url, path, verbose=verbose)
    return open(path) if (cls is None) else cls(path)


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


def download(url, path, verbose=None, blocksize=128*1024):
    """Download a file from a URL, possibly displaying a progress bar.

    Saves the output to the file named by `path`.  If the URL cannot be
    downloaded or the file cannot be written, an IOError is raised.

    Normally, if the standard error output is a terminal, then a
    progress bar is displayed to keep the user entertained.  Specify
    `verbose=True` or `verbose=False` to control this behavior.

    """
    tempname = path + '.download'
    try:
        connection = urlopen(url)
    except Exception as e:
        raise IOError('cannot get {1} because {2}'.format(url, e))
    if verbose is None:
        verbose = sys.stderr.isatty()
    if verbose:
        bar = ProgressBar(path)
        content_length = int(connection.headers.get('content-length', -1))

    # Python open() provides no way to achieve O_CREAT without also
    # truncating the file, which would ruin the work of another process
    # that is trying to download the same file at the same time.  So:

    flags = getattr(os, 'O_BINARY', 0) | os.O_CREAT | os.O_RDWR
    fd = os.open(tempname, flags)
    with os.fdopen(fd, 'wb') as w:
        try:
            if lockf is not None:
                fd = w.fileno()
                lockf(fd, LOCK_EX)           # only one download at a time
                if os.path.exists(path): # did someone else finish first?
                    if os.path.exists(tempname):
                        os.unlink(tempname)
                    return
            w.seek(0)
            length = 0
            while True:
                data = connection.read(blocksize)
                if not data:
                    break
                w.write(data)
                length += len(data)
                if verbose:
                    bar.report(length, content_length)
            w.flush()
            if lockf is not None:
                # On Unix, rename while still protected by the lock.
                try:
                    os.rename(tempname, path)
                except Exception as e:
                    raise IOError('error renaming {0} to {1} - {2}'.format(
                        tempname, path, e))
        except Exception as e:
            raise IOError('error getting {0} - {1}'.format(url, e))
        finally:
            if lockf is not None:
                lockf(fd, LOCK_UN)
    if lockf is None:
        # On Windows, rename here because the file needs to be closed first.
        try:
            os.rename(tempname, path)
        except Exception as e:
            raise IOError('error renaming {0} to {1} - {2}'.format(
                tempname, path, e))

class ProgressBar(object):
    def __init__(self, path):
        self.filename = os.path.basename(path)
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
              end='\n' if (percent == 100) else '', file=sys.stderr)
        sys.stderr.flush()


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
