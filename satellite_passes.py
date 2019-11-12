# %% codecell
from datetime import datetime, timezone, timedelta
from skyfield.api import Loader, Topos
from skyfield.almanac import get_satellite_passes

# %% codecell
load = Loader("./", verbose=True, expire=True)

# Load timescale
ts = load.timescale()

# Load planets, get Earth and set the ground station
ephemeris = load("de421.bsp")
ground_station = Topos('44 N', '11 E', elevation_m=70)

# Load satellites
satellites = load.tle("http://celestrak.com/NORAD/elements/weather.txt")

# Get satellite
satellite = satellites["NOAA 18"]
print(satellite)

# %% codecell
get_satellite_passes(
    ephemeris,
    ground_station,
    satellite,
    ts.utc(datetime.now(timezone.utc)),
    ts.utc(datetime.now(timezone.utc) + timedelta(days=1)),
    alt_deg_thresh=15
)
