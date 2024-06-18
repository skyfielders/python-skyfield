from skyfield.api import load

max_days = 7.0         # download again once 7 days old
name = 'stations.csv'  # custom filename, not 'gp.php'

base = 'https://celestrak.org/NORAD/elements/gp.php'
url = base + '?GROUP=stations&FORMAT=csv'

if not load._exists(name) or load.days_old(name) >= max_days:
    load.download(url, filename=name)
