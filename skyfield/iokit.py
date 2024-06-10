# -*- coding: utf-8 -*-
from __future__ import print_function
import itertools
import os
import errno
import sys
from fnmatch import fnmatch
from pkgutil import get_data
from time import time

import certifi
import numpy as np

from .data import iers
from .curvelib import Splines
from .functions import load_bundled_npy
from .io_timescale import (
    _build_legacy_data,
    parse_deltat_data, parse_deltat_preds, parse_leap_seconds,
)
from .jpllib import SpiceKernel
from .sgp4lib import EarthSatellite
from .timelib import Timescale

try:
    from io import BytesIO
except:
    from StringIO import StringIO as BytesIO

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

try:
    urlopen('', cafile=None)
except TypeError:
    _supports_cafile_argument = False
except ValueError:  # Expected when the URL is an empty string.
    _supports_cafile_argument = True

# If we are running under the built-in IDLE development environment, we
# cannot use '\r' to keep repainting the current line as a progress bar:
_running_IDLE = (sys.stderr.__class__.__name__ == 'PseudoOutputFile')

def _filename_of(url):
    """Return the last path component of a url."""
    return urlparse(url).path.split('/')[-1]

_IERS = 'https://hpiers.obspm.fr/iers/bul/bulc/'
_IERS2 = 'ftp://ftp.iers.org/products/eop/rapid/standard/'
_JPL = 'https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/'
_NAIF_KERNELS = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/'
_NAIF = 'https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/'

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
      (This option is no longer supported.)

    Once a `Loader` is created, it can be called like a function to
    open, or else to download and open, a file whose name it recognizes::

        planets = load('de405.bsp')

    Each loader also supports an attribute and a few methods.

    """
    def __init__(self, directory, verbose=True, expire=False):
        self.directory = os.path.expanduser(directory)
        self.verbose = verbose
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
            'finals2000A.all': _IERS2,
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
        if self.directory == '.':
            return filename
        return os.path.join(self.directory, filename)

    def days_old(self, filename):
        """Return how recently ``filename`` was modified, measured in days."""
        mtime = os.stat(self.path_to(filename)).st_mtime
        seconds = time() - mtime
        return seconds / 86400.0

    def exists(self, filename):
        return os.path.exists(self.path_to(filename))

    def __call__(self, filename, reload=False, backup=False, builtin=False):
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

        parser = _search(self.parsers, filename)
        opener = _search(self.openers, filename)
        if (parser is None) and (opener is None):
            raise ValueError('Skyfield does not know how to open a file'
                             ' named {0!r}'.format(filename))

        if builtin:
            self._log('{0}\n  Parsing builtin file with {1}()',
                      filename, parser.__name__)
            f = BytesIO(get_data('skyfield.data', filename))
            return parser(f)

        path = self._assure(url, filename, reload, backup)

        if parser is not None:
            self._log('  Parsing with {0}()', parser.__name__)
            with open(path, 'rb') as f:
                return parser(f)
        else:
            self._log('  Opening with {0}', opener.__name__)
            return opener(path)

    def _assure(self, url, filename, reload, backup):
        path = self.path_to(filename)
        exists = os.path.exists(path)
        self._log(path)
        if exists:
            self._log('  File already exists')
        if (not exists) or reload:
            if url is None:
                raise ValueError('Skyfield does not know where to download {!r}'
                                 .format(filename))
            self._log('  Downloading {0}', url)
            download(url, path, self.verbose, backup=backup)
        return path

    def _log(self, message, *args):
        self.events.append(message.format(*args))

    def build_url(self, filename):
        """Return the URL Skyfield will try downloading for a given filename.

        Raises ``ValueError`` if Skyfield doesn't know where to get the
        file based on its name.

        """
        base = _search(self.urls, filename)
        if base:
            return base + filename
        raise ValueError("Skyfield doesn't know the URL of {0!r}"
                         .format(filename))

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

    def tle_file(self, url, reload=False, filename=None,
                 ts=None, skip_names=False):
        """Load and parse a TLE file, returning a list of Earth satellites.

        Given a URL or local path to an ASCII text file, this loads a
        series of TLE “Two-Line Element” sets and returns a list of
        :class:`~skyfield.sgp4lib.EarthSatellite` objects for them.
        See :doc:`earth-satellites`.

        See the :meth:`~skyfield.iokit.Loader.open()` method for the
        meaning of the ``reload`` and ``filename`` parameters.

        See the :meth:`parse_tle_file()` function for the meaning of the
        ``ts`` and ``skip_names`` parameters.

        """
        with self.open(url, reload=reload, filename=filename) as f:
            return list(parse_tle_file(f, ts, skip_names))

    def download(self, url, filename=None, backup=False):
        """Download a file, even if it’s already on disk; return its path.

        You can specify the local ``filename`` to which the file will be
        saved; the default is to use the final component of ``url``.
        Set ``backup`` to ``True`` if you want an already-existing file
        moved out of the way instead of overwritten.

        Your operating system may raise any of several errors during a
        download: hostname lookup failure (this is the usual symptom if
        you are disconnected from the Internet); the server refusing the
        connection; and the connection closing mid-download.  Skyfield
        makes no attempt to intercept or interpret these errors — which
        vary by operating system — so your application itself should
        catch network errors if it needs to avoid printing raw Python
        exceptions, or if you want to retry failed downloads.

        """
        if '://' not in url:
            url = self.build_url(url)
        if filename is None:
            filename = urlparse(url).path.split('/')[-1]
        path = self.path_to(filename)
        download(url, path, self.verbose, backup=backup)
        return path

    def open(self, url, mode='rb', reload=False, filename=None, backup=False):
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

        path = self._assure(url, filename, reload, backup)
        return open(path, mode)

    def timescale(self, delta_t=None, builtin=True):
        """Return a `Timescale` built using official Earth rotation data.

        ``delta_t`` — Lets you override the standard ∆T tables by
        providing your own ∆T offset in seconds.  For details, see
        :ref:`custom-delta-t`.

        ``builtin`` — By default, Skyfield uses ∆T and leap second
        tables that it carries internally; to instead load this data
        from files, set this option to ``False``.  For compatibility
        with Skyfield ≤ 1.30, if you have on disk the three files
        ``deltat.data``, ``deltat.preds``, and ``Leap_Second.dat``, then
        Skyfield will load them.  Otherwise, Skyfield will download and
        use ``finals2000A.all`` from the International Earth Rotation
        Service.  For details, see :ref:`downloading-timescale-files`.

        """
        e = self.exists
        if builtin:
            # See "build_arrays.py" for a notes on how these are stored.
            arrays = load_bundled_npy('iers.npz')
            daily_tt = arrays['tt_jd_minus_arange']
            daily_tt += np.arange(len(daily_tt))
            daily_delta_t = (arrays['delta_t_1e7'] / 1e7).round(7)
            delta_t_recent = daily_tt, daily_delta_t
            leap_dates = arrays['leap_dates']
            leap_offsets = arrays['leap_offsets']
        elif e('deltat.data') and e('deltat.preds') and e('Leap_Second.dat'):
            # Avoid changing the meaning of "builtin=False" and
            # surprising the user with a file download, if their
            # previous version of Skyfield already downloaded the three
            # old files we used to rely on.
            deltat_data = self('deltat.data')
            deltat_preds = self('deltat.preds')
            _, leap_second_dat = self('Leap_Second.dat')
            delta_t_recent, leap_dates, leap_offsets = _build_legacy_data(
                deltat_data, deltat_preds, leap_second_dat)
        else:
            url = self.build_url('finals2000A.all')
            with self.open(url) as f:
                utc_mjd, dut1 = iers.parse_dut1_from_finals_all(f)
            daily_tt, daily_delta_t, leap_dates, leap_offsets = (
                iers.build_timescale_arrays(utc_mjd, dut1))
            delta_t_recent = daily_tt, daily_delta_t

        if delta_t is not None:
            delta_t_recent = Splines([0, 1, delta_t])

        return Timescale(delta_t_recent, leap_dates, leap_offsets)

    @property
    def log(self):
        return '\n'.join(self.events)

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

def parse_tle_file(lines, ts=None, skip_names=False):
    """Parse lines of TLE satellite data, yielding a sequence of satellites.

    Given a sequence ``lines`` of byte strings (which can be an open
    binary file, which acts like a sequence of lines in Python), this
    routine yields an :class:`~skyfield.sgp4lib.EarthSatellite` for each
    pair of adjacent lines that start with ``"1 "`` and ``"2 "`` and
    have 69 or more characters each.  If the line preceding a TLE is not
    part of another TLE, it is used as the satellite’s ``.name``.

    If you pass a ``ts`` timescale, Skyfield will use it to build the
    ``.epoch`` date attribute on each satellite; otherwise a timescale
    derived from Skyfield’s built-in leap second files will be used.

    If for a particular file you see random lines of text being
    interpreted as satellite names, set ``skip_names`` to ``True`` and
    Skyfield will not try to store satellite names.

    See :doc:`earth-satellites` for details.  An exception is raised if
    the attempt to parse a pair of candidate lines as TLE lines fails.

    """
    b0 = b1 = b''
    for b2 in lines:
        if (b2.startswith(b'2 ') and len(b2) >= 69 and
            b1.startswith(b'1 ') and len(b1) >= 69):

            if not skip_names and b0:
                b0 = b0.rstrip(b' \n\r')
                if b0.startswith(b'0 '):
                    b0 = b0[2:]  # Spacetrack 3-line format
                name = b0.decode('ascii')
            else:
                name = None

            line1 = b1.decode('ascii')
            line2 = b2.decode('ascii')
            yield EarthSatellite(line1, line2, name, ts)

            b0 = b1 = b''
        else:
            b0 = b1
            b1 = b2

def download(url, path, verbose=None, blocksize=128*1024, backup=False):
    """Download a file from a URL, possibly displaying a progress bar.

    Saves the output to the file named by `path`.  If the URL cannot be
    downloaded or the file cannot be written, an ``IOError`` is raised.

    Normally, if the standard error output is a terminal, then a
    progress bar is displayed to keep the user entertained.  Specify
    `verbose=True` or `verbose=False` to override this behavior.

    """
    try:
        if create_default_context is not None:
            ssl_context = create_default_context(cafile=certifi.where())
            connection = urlopen(url, context=ssl_context)
        elif _supports_cafile_argument:
            connection = urlopen(url, cafile=certifi.where())
        else:
            connection = urlopen(url)  # Very old Python: no certificate check.
    except Exception as e:
        e2 = IOError('cannot download {0} because {1}'.format(url, e))
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

    # Claim our own unique download filename.

    tempbase = tempname = path + '.download'
    flags = getattr(os, 'O_BINARY', 0) | os.O_CREAT | os.O_EXCL | os.O_RDWR
    i = 1
    while True:
        try:
            fd = os.open(tempname, flags, 0o666)
        except OSError as e:  # "FileExistsError" is not supported by Python 2
            if e.errno != errno.EEXIST:
                raise
            i += 1
            tempname = '{0}{1}'.format(tempbase, i)
        else:
            break

    # Download to the temporary filename.

    with os.fdopen(fd, 'wb') as w:
        try:
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
        except Exception as e:
            raise IOError('error getting {0} - {1}'.format(url, e))

    # Move the original out of the way, if requested.

    if os.path.exists(path) and backup:
        _rename_original(path)

    # Rename the temporary file to the destination name.

    try:
        _replace(tempname, path)
    except Exception as e:
        raise IOError('error renaming {0} to {1} - {2}'.format(
            tempname, path, e))

def _rename_original(path):
    for n in itertools.count(1):
        prefix, suffix = path.rsplit('.', 1)
        backup_path = '{0}.old{1}.{2}'.format(prefix, n, suffix)
        if not os.path.exists(backup_path):
            break
    os.rename(path, backup_path)

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
