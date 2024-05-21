from skyfield.api import load

max_days = 7.0  # replace files that are one week old

name = 'stations.csv'
base = 'https://celestrak.org/NORAD/elements/gp.php'
url = base + '?GROUP=stations&FORMAT=csv'

if not load.exists(name) or load.days_old(name) >= max_days:
    load.download(url, filename=name)
