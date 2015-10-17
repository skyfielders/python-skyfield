import os
import numpy as np
from datetime import datetime, timedelta

from .jpllib import SpiceKernel

try:
    from urllib.request import urlopen
except:
    from urllib2 import urlopen

_missing = object()


def load(filename):
    if filename.endswith('.bsp'):
        if not os.path.exists(filename):
            url = url_for(filename)
            download(url)
        return SpiceKernel(filename)
    else:
        raise ValueError('Skyfield does not recognize that file extension')


def url_for(filename):
    if filename.endswith('.bsp'):
        url = 'http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/'
        if filename.startswith('de'):
            url += 'planets/'
            if filename < 'de430':
                url += 'a_old_versions/'
                if filename == 'de423.bsp':
                    url += 'de423_for_mercury_and_venus/'
        elif filename.startswith('jup'):
            url += 'satellites/'
        return url + filename


def download(url):
    filename = url.rstrip('/').rsplit('/', 1)[1]
    if not os.path.exists(filename):
        try:
            data = urlopen(url).read()
        except Exception as e:
            raise IOError('error getting {} - {}'.format(url, e))
        with open(filename, 'wb') as f:
            f.write(data)
    return filename


class Cache(object):
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
    min_old = datetime.now()-timedelta(days=days_old)
    modified = datetime.fromtimestamp(os.path.getmtime(filename))
    return modified < min_old
