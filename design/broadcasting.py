from numpy import array, max
from skyfield.api import load
from skyfield.constants import C_AUDAY
from skyfield.functions import length_of, _reconcile

def main():
    ts = load.timescale()

    cfltt = _correct_for_light_travel_time   # the original
    #cfltt = _correct_for_light_travel_time2   # possible replacement

    print('==== One time, one observer, one target ====')

    t = ts.utc(2021, 2, 19, 13, 47)
    planets = load('de421.bsp')
    observer = planets['earth'].at(t)
    target = planets['mars']
    r, v, t2, light_time = cfltt(observer, target)
    print(t.shape, observer.position.km.shape, r.shape)

    print('==== N times, one observer, one target ====')

    t = ts.utc(2021, 2, 19, 13, [46, 47, 48, 49])
    planets = load('de421.bsp')
    observer = planets['earth'].at(t)
    target = planets['mars']
    r, v, t2, light_time = cfltt(observer, target)
    print(t.shape, observer.position.km.shape, '->', r.shape)
    print('Here is where the planet wound up:')
    print(r)

    # The above maneuvers work fine even with the old version of the
    # routine.  But to proceed from here, we need to switch.

    print('==== N times, one observers, M targets ====')

    t = ts.utc(2021, 2, 19, 13, [46, 47, 48, 49])
    planets = load('de421.bsp')

    earth = planets['earth']
    if 0:
        # Turn Earth into two observation positions.
        earth = planets['earth']
        earth._at = build_multi_at(earth._at)
    observer = earth.at(t)
    print('observer', observer.position.au.shape)

    target = planets['mars']
    target._at = build_multi_at(target._at)  # Turn Mars into 2 planets.
    print('target', target.at(t).position.au.shape)

    t = ts.tt(t.tt[:,None])  # What if we add a dimension to t?
    print('t', t.shape)

    r, v, t2, light_time = _correct_for_light_travel_time2(observer, target)
    print(t.shape, observer.position.km.shape, '->', r.shape)

    print('Does it look like a second planet 1 AU away at the same 4 times?')
    print('First planet:')
    print(r[:,:,0])
    print('Second planet:')
    print(r[:,:,1])

offset = array([0, 1])[None,None,:]  # Dimensions: [xyz, time, offset]

def build_multi_at(_at):
    # Take a single planet position, and pretend that really there are
    # two targets returning their positions (like comets or asteroids).

    def wrapper(t):
        tposition, tvelocity, gcrs_position, message = _at(t)
        print('_at() t.shape', t.shape)
        if len(t.shape) < 2:
            # If the time lacks a second dimension, then let's expand
            # position and velocity by adding a new dimension (...,2) at
            # the bottom, as though there were two planets.
            tposition = tposition[:,:,None] + offset
            tvelocity = tvelocity[:,:,None] + offset
        else:
            # Otherwise, the time's extra dimension will already have us
            # producing what looks like two objects in the bottom
            # dimension!  We can simply apply our position offset.
            assert tposition.shape[-1] == tvelocity.shape[-1] == 2
            tposition = tposition + offset
            tvelocity = tvelocity + offset
        return tposition, tvelocity, gcrs_position, message
    return wrapper

# The original.

def _correct_for_light_travel_time(observer, target):
    """Return a light-time corrected astrometric position and velocity.

    Given an `observer` that is a `Barycentric` position somewhere in
    the solar system, compute where in the sky they will see the body
    `target`, by computing the light-time between them and figuring out
    where `target` was back when the light was leaving it that is now
    reaching the eyes or instruments of the `observer`.

    """
    t = observer.t
    ts = t.ts
    whole = t.whole
    tdb_fraction = t.tdb_fraction

    cposition = observer.position.au
    cvelocity = observer.velocity.au_per_d

    tposition, tvelocity, gcrs_position, message = target._at(t)

    distance = length_of(tposition - cposition)
    light_time0 = 0.0
    for i in range(10):
        light_time = distance / C_AUDAY
        delta = light_time - light_time0
        if abs(max(delta)) < 1e-12:
            break

        # We assume a light travel time of at most a couple of days.  A
        # longer light travel time would best be split into a whole and
        # fraction, for adding to the whole and fraction of TDB.
        t2 = ts.tdb_jd(whole, tdb_fraction - light_time)

        tposition, tvelocity, gcrs_position, message = target._at(t2)
        distance = length_of(tposition - cposition)
        light_time0 = light_time
    else:
        raise ValueError('light-travel time failed to converge')
    return tposition - cposition, tvelocity - cvelocity, t, light_time

# Try allowing vectors while keeping things "right side up" with x,y,z
# at the top dimension.

def sub(a, b):
    return (a.T - b.T).T

def _correct_for_light_travel_time2(observer, target):
    """Return a light-time corrected astrometric position and velocity.

    Given an `observer` that is a `Barycentric` position somewhere in
    the solar system, compute where in the sky they will see the body
    `target`, by computing the light-time between them and figuring out
    where `target` was back when the light was leaving it that is now
    reaching the eyes or instruments of the `observer`.

    """
    t = observer.t
    ts = t.ts
    whole = t.whole
    tdb_fraction = t.tdb_fraction

    cposition = observer.position.au
    cvelocity = observer.velocity.au_per_d

    tposition, tvelocity, gcrs_position, message = target._at(t)

    distance = length_of(sub(tposition, cposition))

    print('distance', distance.shape)  # (t, targets)

    light_time0 = 0.0
    for i in range(10):
        light_time = distance / C_AUDAY   # GOOD: scalar
        delta = light_time - light_time0  # GOOD first time: scalar; 2nd: ?
        if abs(max(delta)) < 1e-12:
            break

        # We assume a light travel time of at most a couple of days.  A
        # longer light travel time would best be split into a whole and
        # fraction, for adding to the whole and fraction of TDB.
        print('tdb_fraction', tdb_fraction.shape, '[ORIGINAL]')
        print('light_time', light_time.shape)
        diff = sub(tdb_fraction, light_time)
        print('sub()', diff.shape)  # (4,2)? YES!!!
        print('whole', whole.shape)  # (4,)? YES!!!
        # whole, diff = _reconcile(whole, diff)  # Winds up not needed?
        # print('sub()', diff.shape)  # (4,2)
        # print('whole', whole.shape)  # (4,1)
        t2 = ts.tdb_jd(whole, diff)
        print('t2', t2.shape)  # (4, 1)? Why not (4, 2)? Because it prints top.

        tposition, tvelocity, gcrs_position, message = target._at(t2)
        print('tposition', tposition.shape)  # Needs to be 3,4,2
        distance = length_of(sub(tposition, cposition))
        light_time0 = light_time
        #exit()
    else:
        raise ValueError('light-travel time failed to converge')
    return sub(tposition, cposition), sub(tvelocity, cvelocity), t, light_time

_reconcile  # So CI will think we used it, whether use above is commented or not

if __name__ == '__main__':
    main()
