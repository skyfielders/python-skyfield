import requests
import os
from datetime import datetime, timedelta


class Downloader(object):
    def __init__(self, cache_path, days_old=0):
        self.cache_path = cache_path
        self.days_old = days_old

    def get(self, url, days_old=None):
        filename = url[url.rindex('/') + 1:]
        path = os.path.join(self.cache_path, filename)
        if days_old is None:
            days_old = self.days_old
        download_file(url, path, days_old)
        return open(path)


def download_file(url, filename, days_old=0):
    if os.path.exists(filename):
        if not is_days_old(filename, days_old):
            return

    response = requests.get(url, stream=True)
    f = open(filename, 'wb')
    for chunk in response.iter_content(1024):
        f.write(chunk)

    f.close()

def is_days_old(filename, days_old):
    min_old = datetime.now()-timedelta(days=days_old)
    modified = datetime.fromtimestamp(os.path.getmtime(filename))
    return modified < min_old

