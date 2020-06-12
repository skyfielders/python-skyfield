# -*- coding: utf-8 -*-
from __future__ import print_function
import itertools
import os
import errno
import numpy as np
import sys
from datetime import date
from fnmatch import fnmatch
from time import time

import certifi

from .io_timescale import (_build_builtin_timescale, parse_deltat_data,
                           parse_deltat_preds, parse_leap_seconds)
from .jpllib import SpiceKernel
from .sgp4lib import EarthSatellite
from .timelib import Timescale

today = date.today

try:
    from fcntl import LOCK_EX, LOCK_UN, lockf
except:
    lockf = None

if sys.version_info >= (3, 3):
    _replace = os.replace
else:
    _replace = os.rename  # Raises OSError on Windows if destination exists

try:
    from ssl import create_default_context
except ImportError:
    create_default_context = None

try:
    from urllib.parse import urlparse
    from urllib.request import urlopen
except:
    from urlparse import urlparse
    from urllib2 import urlopen

# If we are running under the built-in IDLE development environment, we
# cannot use '\r' to keep repainting the current line as a progress bar:
_running_IDLE = (sys.stderr.__class__.__name__ == 'PseudoOutputFile')

def _filename_of(url):
    """Return the last path component of a url."""
    return urlparse(url).path.split('/')[-1]

_IERS = 'https://hpiers.obspm.fr/iers/bul/bulc/'
_JPL = 'ftp://ssd.jpl.nasa.gov/pub/eph/planets/bsp/'
_NAIF_KERNELS = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/'
_NAIF = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/'
_CDDIS = 'ftp://cddis.nasa.gov/products/iers/'

def _open_binary(path):
    return open(path, mode='rb')

class Loader(object):
    """A tool for downloading and opening astronomical data files.

    A default `Loader` that saves data files to the current working
    directory can be imported directly from the Skyfield API::

        from skyfield.api import load

    But users can also create a `Loader` of their own, if there is
    another directory they want data files saved to, or if they want to
    specify different options.  The directory is created automatically
    if it does not yet exist::

        from skyfield.api import Loader
        load = Loader('~/skyfield-data')

    The options are:

    ``verbose``
      If set to ``False``, then the loader will not print a progress bar
      to the screen each time it downloads a file.  (If the standard
      output is not a TTY, then no progress bar is printed anyway.)

    ``expire``
      If set to ``False``, then Skyfield will always use an existing
      file on disk, instead of expiring files that are out of date and
      replacing them with newly downloaded copies.

    Once a `Loader` is created, it can be called like a function to
    open, or else to download and open, a file whose name it recognizes::

        planets = load('de405.bsp')

    Each loader also supports an attribute and a few methods.

    """
    def __init__(self, directory, verbose=True, expire=True):
        self.directory = os.path.expanduser(directory)
        self.verbose = verbose
        self.expire = expire
        self.events = []
        try:
            os.makedirs(self.directory)
        except OSError as e:
            if e.errno != errno.EEXIST and not os.path.isdir(self.directory):
                raise

        # Each instance gets its own copy of these data structures,
        # instead of sharing a single copy, so users can edit them
        # without changing the behavior of other Loader objects:

        self.urls = {
            'deltat.data': _CDDIS,
            'deltat.preds': _CDDIS,
            'Leap_Second.dat': _IERS,
            'moon_080317.tf': _NAIF_KERNELS + 'fk/satellites/',
            'moon_pa_de421_1900-2050.bpc': _NAIF_KERNELS + 'pck/',
            'pck00008.tpc': _NAIF_KERNELS + 'pck/a_old_versions/',
            '.bsp': [
                ('jup*.bsp', _NAIF),
                ('*.bsp', _JPL),
            ],
        }
        self.parsers = {
            'deltat.data': parse_deltat_data,
            'deltat.preds': parse_deltat_preds,
            'Leap_Second.dat': parse_leap_seconds,
        }
        self.openers = {
            # Old-fashioned: auto-create objects, leaving readers and
            # code tools guessing what kind of object we have returned.
            '.bsp': [
                ('*.bsp', SpiceKernel),
            ],
            # New approach: just return open files, which callers can
            # then pass to the right class, making the class visible in
            # the code to both human readers and their IDEs.
            '.bpc': [('*', _open_binary)],
            '.tpc': [('*', _open_binary)],
            '.tf': [('*', _open_binary)],
        }

    def path_to(self, filename):
        """Return the path to ``filename`` in this loader's directory."""
        return os.path.join(self.directory, filename)

    def __call__(self, filename):
        """Open the given file, downloading it first if necessary."""
        if '://' in filename:
            url = filename
            filename = urlparse(url).path.split('/')[-1]
        # Should this API accept full path names? It might look like:
        # elif os.sep in filename:
        #     os.path.expanduser(directory)
        #     path = filename
        #     filename = os.path.basename(path)
        #     url = _search(self.urls, filename)
        #     directory =
        else:
            url = _search(self.urls, filename)
            if url:
                url += filename

        path = self.path_to(filename)
        parser = _search(self.parsers, filename)
        opener = _search(self.openers, filename)
        if (parser is None) and (opener is None):
            raise ValueError('Skyfield does not know how to open a file'
                             ' named {0!r}'.format(filename))

        if os.path.exists(path):
            self._log('Already exists: {0}', path)
            if parser is not None:
                self._log('  Parsing with: {0}()', parser.__name__)
                with open(path, 'rb') as f:
                    expiration_date, data = parser(f)
                if not self.expire:
                    self._log('  Ignoring expiration: {0}', expiration_date)
                    return data
                if expiration_date is None:
                    self._log('  Does not specify an expiration date')
                    return data
                if today() <= expiration_date:
                    self._log('  Does not expire til: {0}', expiration_date)
                    return data
                self._log('  Expired on: {0}', expiration_date)
                for n in itertools.count(1):
                    prefix, suffix = filename.rsplit('.', 1)
                    backup_name = '{0}.old{1}.{2}'.format(prefix, n, suffix)
                    if not os.path.exists(self.path_to(backup_name)):
                        break
                self._log('  Renaming to: {0}', backup_name)
                os.rename(self.path_to(filename), self.path_to(backup_name))
            else:
                # Currently, openers have no concept of expiration.
                self._log('  Opening with: {0}', opener.__name__)
                return opener(path)

        if url is None:
            raise ValueError('Skyfield does not know where to download {!r}'
                             .format(filename))

        self._log('  Downloading: {0}', url)
        download(url, path, self.verbose)

        if parser is not None:
            self._log('  Parsing with: {0}()', parser.__name__)
            with open(path, 'rb') as f:
                expiration_date, data = parser(f)
            return data
        else:
            self._log('  Opening with: {0}', opener.__name__)
            return opener(path)

    def _log(self, message, *args):
        self.events.append(message.format(*args))

    def build_url(self, filename):
        """Return the URL Skyfield will try downloading for a given filename."""
        base = _search(self.urls, filename)
        return None if (base is None) else base + filename

    def tle(self, url, reload=False, filename=None):
        """Load and parse a satellite TLE file.

        DEPRECATED: in a misguided attempt to be overly convenient, this
        routine builds an unweildy dictionary of satellites with keys of
        two different Python types: integer keys for satellite numbers,
        and string keys for satellite names. It even lists satellites
        like ``ISS (ZARYA)`` twice, in case the user wants to look them
        up by a single name like ``ZARYA``.  What a mess.  Users should
        instead call the simple ``tle_file()`` method, and themselves
        build any dictionaries they need.

        See the :meth:`~skyfield.iokit.Loader.open()` documentation for
        the meaning of the ``reload`` and ``filename`` parameters.

        """
        d = {}
        with self.open(url, reload=reload, filename=filename) as f:
            for names, sat in parse_tle(f):
                d[sat.model.satnum] = sat
                for name in names:
                    d[name] = sat
        return d

    def tle_file(self, url, reload=False, filename=None, ts=None):
        """Load and parse a TLE file, returning a list of Earth satellites.

        Given a URL or local path to an ASCII text file, this loads a
        series of TLE “Two-Line Element” sets and returns a list of
        :class:`~skyfield.sgp4lib.EarthSatellite` objects for them.
        See :doc:`earth-satellites`.

        See the :meth:`~skyfield.iokit.Loader.open()` documentation for
        the meaning of the ``reload`` and ``filename`` parameters.

        """
        with self.open(url, reload=reload, filename=filename) as f:
            return list(parse_tle_file(f, ts))

    def open(self, url, mode='rb', reload=False, filename=None):
        """Open a file, downloading it first if it does not yet exist.

        Unlike when you call a loader directly like ``my_loader()``,
        this ``my_loader.open()`` method does not attempt to parse or
        interpret the file; it simply returns an open file object.

        The ``url`` can be either an external URL, or else the path to a
        file on the current filesystem.  A relative path will be assumed
        to be relative to the base directory of this loader object.

        If a URL was provided and the ``reload`` parameter is true, then
        any existing file will be removed before the download starts.

        The ``filename`` parameter lets you specify an alternative local
        filename instead of having the filename extracted from the final
        component of the URL.

        """
        if '://' not in url:
            path_that_might_be_relative = url
            path = os.path.join(self.directory, path_that_might_be_relative)
            return open(path, mode)
        if filename is None:
            filename = urlparse(url).path.split('/')[-1]
        path = self.path_to(filename)
        if reload and os.path.exists(path):
            os.remove(path)
        if not os.path.exists(path):
            download(url, path, self.verbose)
        return open(path, mode)

    def timescale(self, delta_t=None, builtin=False):
        """Open or download three time scale files, returning a `Timescale`.

        This method is how most Skyfield users build a `Timescale`
        object, which is necessary for building `Time` objects.  The
        safest approach is::

            ts = load.timescale(builtin=True)

        This avoids downloading any files by using built-in copies of
        them instead.  The problem is that the files distributed with
        any particular version of Skyfield will go gradually out of date
        and you will start missing leap seconds.  To instead download
        current files, omit the ``builtin`` option::

            ts = load.timescale()

        UT1 is tabulated by the United States Naval Observatory files
        ``deltat.data`` and ``deltat.preds``, while UTC is defined by
        ``Leap_Second.dat`` from the International Earth Rotation
        Service.

        """
        if builtin:
            return _build_builtin_timescale()

        if delta_t is not None:
            delta_t_recent = np.array(((-1e99, 1e99), (delta_t, delta_t)))
        else:
            try:
                data = self('deltat.data')
                preds = self('deltat.preds')
            except IOError as e:
                e.args = (e.args[0] + _TIMESCALE_IO_ADVICE,) + e.args[1:]
                raise
            data_end_time = data[0, -1]
            i = np.searchsorted(preds[0], data_end_time, side='right')
            delta_t_recent = np.concatenate([data, preds[:,i:]], axis=1)
        try:
            leap_dates, leap_offsets = self('Leap_Second.dat')
        except IOError as e:
            e.args = (e.args[0] + _TIMESCALE_IO_ADVICE,) + e.args[1:]
            raise
        return Timescale(delta_t_recent, leap_dates, leap_offsets)

    @property
    def log(self):
        return '\n'.join(self.events)

_TIMESCALE_IO_ADVICE = """

You can avoid this error by passing `timescale(builtin=True)` which
makes Skyfield use built-in copies of the timescale files instead of
downloading new ones.  The built-in leap second and Earth rotation files
will gradually go out of date unless you periodically upgrade Skyfield."""

def _search(mapping, filename):
    """Search a Loader data structure for a filename."""
    result = mapping.get(filename)
    if result is not None:
        return result
    name, ext = os.path.splitext(filename)
    result = mapping.get(ext)
    if result is not None:
        for pattern, result2 in result:
            if fnmatch(filename, pattern):
                return result2
    return None


def load_file(path):
    """Open a file on your local drive, using its extension to guess its type.

    This routine only works on ``.bsp`` ephemeris files right now, but
    will gain support for additional file types in the future. ::

        from skyfield.api import load_file
        planets = load_file('~/Downloads/de421.bsp')

    """
    path = os.path.expanduser(path)
    base, ext = os.path.splitext(path)
    if ext == '.bsp':
        return SpiceKernel(path)
    raise ValueError('unrecognized file extension: {}'.format(path))

def parse_tle(fileobj):
    """Parse a file of TLE satellite element sets.

    DEPRECATED: this routine is overly complicated, doing extra work to
    try to guess several ways in which the user might want to look up
    satellites by name.  Use ``parse_tle_file()`` instead.

    TODO: convert this into a wrapper around ``parse_tle_file()``.

    """
    b0 = b1 = b''
    for b2 in fileobj:
        if (b1.startswith(b'1 ') and len(b1) >= 69 and
            b2.startswith(b'2 ') and len(b2) >= 69):

            b0 = b0.rstrip(b'\n\r')
            if len(b0) == 24:   # Celestrak
                name = b0.decode('ascii').rstrip()
                names = [name]
            elif b0.startswith(b'0 '): # Spacetrack 3-line format
                name = b0[2:].decode('ascii').rstrip()
                names = [name]
            else:
                name = None
                names = ()

            line1 = b1.decode('ascii')
            line2 = b2.decode('ascii')
            sat = EarthSatellite(line1, line2, name)

            if name and ' (' in name:
                # Given a name like `ISS (ZARYA)` or `HTV-6 (KOUNOTORI
                # 6)`, also support lookup by the name inside or outside
                # the parentheses.
                short_name, secondary_name = name.split(' (')
                secondary_name = secondary_name.rstrip(')')
                names.append(short_name)
                names.append(secondary_name)

            yield names, sat

        b0 = b1
        b1 = b2

def parse_tle_file(lines, ts=None):
    """Parse TLE satellite elements, yielding a sequence of satellites.

    Given a sequence ``lines`` of byte strings that contain TLE
    “Two-Line Element” sets of satellite orbital elements, yields an
    :class:`~skyfield.sgp4lib.EarthSatellite` for each pair of adjacent
    lines that start with ``b"1 "`` and ``b"2 "`` and have 69 or more
    characters each.  If the preceding line looks like a name, it is set
    as the satellite’s ``.name``.  See :doc:`earth-satellites`.

    An exception is raised if the attempt to parse a pair of candidate
    lines as TLE elements fails.

    """
    b0 = b1 = b''
    for b2 in lines:
        if (b1.startswith(b'1 ') and len(b1) >= 69 and
            b2.startswith(b'2 ') and len(b2) >= 69):

            b0 = b0.rstrip(b'\n\r')
            if len(b0) == 24:   # Celestrak
                name = b0.decode('ascii').rstrip()
            elif b0.startswith(b'0 '): # Spacetrack 3-line format
                name = b0[2:].decode('ascii').rstrip()
            else:
                name = None

            line1 = b1.decode('ascii')
            line2 = b2.decode('ascii')
            yield EarthSatellite(line1, line2, name, ts)

        b0 = b1
        b1 = b2

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
        if create_default_context is None:
            connection = urlopen(url, cafile=certifi.where())
        else:
            ssl_context = create_default_context(cafile=certifi.where())
            connection = urlopen(url, context=ssl_context)
    except Exception as e:
        e2 = IOError('cannot get {0} because {1}'.format(url, e))
        e2.__cause__ = None
        raise e2
    if verbose is None:
        verbose = sys.stderr.isatty()

    bar = None
    if verbose:
        if _running_IDLE:
            print('Downloading {0} ...'.format(os.path.basename(path)),
                  file=sys.stderr)
        else:
            bar = ProgressBar(path)
            content_length = int(connection.headers.get('content-length', -1))

    # Python open() provides no way to achieve O_CREAT without also
    # truncating the file, which would ruin the work of another process
    # that is trying to download the same file at the same time.  So:

    flags = getattr(os, 'O_BINARY', 0) | os.O_CREAT | os.O_RDWR
    fd = os.open(tempname, flags, 0o666)
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
                if bar is not None:
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
            _replace(tempname, path)
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
