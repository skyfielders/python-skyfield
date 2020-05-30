"""Verdict: light-time delay from the Sun can be neglected in is_sunlit()!

For the Sun's position from Earth, does calling observe() really make
enough difference to justify the expense?  Here we take two approaches
to answering the question: we compare the difference every day over 40
years, and then we do a back-of-the-envelope estimate of how big we
might have expected the effect to be.  The two approaches agree!  The
maximum difference is around 10 mas.

What difference does that make for a satellite?  Let's take the ISS.
With its orbital period of 92 minutes, it sees the Earth swing in a full
circle around the sky in that amount of time.  That's 360/90 = 4 degrees
per minute (!) = 240 arcseconds per second.  At that speed, a difference
of 10 mas in the Sun's position would at most hasten or delay the moment
of sunrise for the ISS by 40 microseconds, which is far below the
accuracy of TLE position predictions and can thus be safely incurred.

"""
from skyfield import api
from skyfield.api import load

ts = load.timescale()
eph = load('de421.bsp')
sun = eph['sun']
earth = eph['earth']

t = ts.utc(2000, 1, range(40 * 365))
s1 = earth.at(t).observe(sun)
s2 = (sun - earth).at(t)
print('Milliarcseconds (mas) difference:', s1.separation_from(s2).mas().max())

print()
print('Does that make physical sense?')

# The Sun orbits around the Solar System barycenter which is usually
# inside the Sun's radius but occasionally a bit outside of it.  So we
# can very roughly imagine the Sun's orbit as its own circumference,
# give or take.
solar_radius_km = 696340
distance_sun_travels_in_one_orbit = solar_radius_km * api.tau

# It takes the Sun more than a decade to travel that path, as its orbit
# is roughly the opposite of Jupiter's (which takes 12 years to circle
# the Sun).  So it turns out that it travels a bit slowly.
sun_km_per_s = distance_sun_travels_in_one_orbit / 10 / 365.25 / 24 / 60 / 60
print('Sun km/s:', sun_km_per_s)

light_delay_seconds = s2[0].position.length().light_seconds()
print('Sample light delay from Sun to Earth (seconds):', light_delay_seconds)

print('How far does the Sun move in that time?')

travel_km = light_delay_seconds * sun_km_per_s
print('Sun km moved during that light travel time:', travel_km)

print('What angle does that many kilometers subtend from Earth?')

earth_sun_distance_km = 150e6
travel_angle = api.Angle(radians = travel_km / earth_sun_distance_km)
print('Angle traveled by sun in arcseconds:', travel_angle.arcseconds())
print('Angle traveled by sun in mas:', travel_angle.mas())

print()
print(__doc__.rstrip())
