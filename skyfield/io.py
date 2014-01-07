import requests
import os
from datetime import datetime, timedelta


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

