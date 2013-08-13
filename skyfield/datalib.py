import requests
import os
from datetime import datetime, timedelta

class DownloadFile(object):

    def __init__(self, url, filename, days_old=0):
        response = requests.get(url, stream=True)
        f = open(filename, 'w')
        for chunk in response.iter_content(1024):
            f.write(chunk)

def is_days_old(filename, days_old):
    min_old = datetime.now()-timedelta(days=days_old)
    modified = datetime.fromtimestamp(os.path.getmtime(filename))
    return modified < min_old

