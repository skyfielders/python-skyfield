from __future__ import print_function
import itertools
import os
import numpy as np
import sys
from datetime import date, datetime, timedelta
from pkgutil import get_data
from time import time

from .jpllib import SpiceKernel
from .timelib import julian_date

today = date.today

try:
    from fcntl import LOCK_EX, LOCK_UN, lockf
except:
    lockf = None

try:
    from io import BytesIO
except:
    from StringIO import StringIO as BytesIO

try:
    from urllib.parse import urlparse
    from urllib.request import urlopen
except:
    from urlparse import urlparse
    from urllib2 import urlopen

_missing = object()


class Cache(object):
    def __init__(self, directory):
        self.directory = directory

    def path_of(self, filename):
        return os.path.join(self.directory, filename)

    def load(self, filename):
        url, parser = _urls.get(filename, (None, None))
        if url is not None:
            expiration_date, data = parser(load(url, self.directory))
            if expiration_date < today():
                for n in itertools.count(1):
                    prefix, suffix = filename.rsplit('.', 1)
                    backup_name = '{}.old.{}'.format(prefix, n, suffix)
                    if not os.path.exists(backup_name):
                        break
                os.rename(self.path_of(filename), self.path_of(backup_name))
                expiration_date, data = parser(load(url, self.directory))
            return data
        return load(filename, self.directory)


def parse_deltat(text):
    array = np.loadtxt(text)
    year, month, day = array[-1,:3].astype(int)
    expiration_date = date(year + 1, month, day)
    year, month, day, delta_t = array.T
    data = np.array((julian_date(year, month, day), delta_t))
    return expiration_date, data


def _filename_of(url):
    return urlparse(url).path.split('/')[-1]

_urls = dict((_filename_of(url), (url, parser)) for url, parser in (
    ('http://maia.usno.navy.mil/ser7/deltat.data', parse_deltat),
    ('http://maia.usno.navy.mil/ser7/leapsec.dat', None),
    ))


def load_bundled_npy(filename):
    data = get_data('skyfield', 'data/{}.npy'.format(filename))
    return np.load(BytesIO(data))


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
    return open(path, 'rb') if (cls is None) else cls(path)


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
        raise IOError('cannot get {0} because {1}'.format(url, e))
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


def is_days_old(filename, days_old):
    """Early experiment in being sensitive to file dates."""
    min_old = datetime.now()-timedelta(days=days_old)
    modified = datetime.fromtimestamp(os.path.getmtime(filename))
    return modified < min_old
