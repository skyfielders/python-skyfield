from numpy import array, max
from skyfield.api import load
from skyfield.constants import C_AUDAY
from skyfield.functions import length_of, _reconcile

def main():
    ts = load.timescale()

    cfltt = _correct_for_light_travel_time   # the original
    #cfltt = _correct_for_light_travel_time2   # possible replacement
    fcfltt = _correct_for_light_travel_time4   # possible replacement

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

    print('==== N times, one observers, M targets [TRY: RIGHT SIDE UP] ====')

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

    # t = ts.tt(t.tt[:,None])  # What if we add a dimension to t?
    # print('t', t.shape)

    r, v, t2, light_time = _correct_for_light_travel_time2(observer, target)
    print(t.shape, observer.position.km.shape, '->', r.shape)

    print('Does it look like a second planet 1 AU away at the same 4 times?')
    print('First planet:')
    print(r[:,:,0])
    print('Second planet:')
    print(r[:,:,1])

    print('==== N times, one observers, M targets [TRY: UPSIDE DOWN] ====')

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

    # t = ts.tt(t.tt[:,None])  # What if we add a dimension to t?
    # print('t', t.shape)

    r, v, t2, light_time = _correct_for_light_travel_time3(observer, target)
    print(t.shape, observer.position.km.shape, '->', r.shape)

    print('Does it look like a second planet 1 AU away at the same 4 times?')
    print('First planet:')
    print(r[:,:,0])
    print('Second planet:')
    print(r[:,:,1])

    # Okay, so, the whole problem of broadcasting: is the only reason it
    # comes up because we have more targets than observers?  What if we
    # ask folks to expand the observers array dimensions as well, so
    # that NumPy is never faced with arrays of different numbers of
    # dimensions, which is what triggers its unfortunate decision to
    # start matching at the ends of the arrays instead of their
    # beginnings?

    print('==== N times, one observers, M targets'
          ' [TRY: EXTRA OBSERVER DIMENSION] ====')

    # So: time has its usual extra dimension to accommodate 4 times.
    t = ts.utc(2021, 2, 19, 13, [46, 47, 48, 49])
    planets = load('de421.bsp')

    earth = planets['earth']
    observer = earth.at(t)
    print('observer', observer.position.au.shape)
    print('Supplementing dimensions')
    observer.position.au = observer.position.au[:,:,None]
    observer.velocity.au_per_d = observer.velocity.au_per_d[:,:,None]
    print('observer', observer.position.au.shape)

    target = planets['mars']
    target._at = build_multi_at(target._at)  # Turn Mars into 2 planets.
    print('target', target.at(t).position.au.shape)

    # TODO: if we work out how we can expand the time's dimensions, as
    # in the following line, before earth.at(t) without raising an
    # exception, then the observer time would not be lacking a dimension
    # that would need to be expanded in

    t = ts.tt(t.tt[:,None])  # What if we add a dimension to t?
    # print('t !!!!!!!!!!!!!!', t.shape)
    # print('t !!!!!!!!!!!!!!', t.whole.shape)
    # print('t !!!!!!!!!!!!!!', t.tdb_fraction.shape)

    r, v, t2, light_time = _correct_for_light_travel_time4(observer, target)
    print(t.shape, observer.position.au.shape, '->', r.shape)

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
    #
    # This version makes subtraction work by replacing normal NumPy
    # broadcasting subtraction with a sub() function of our own that
    # broadcasts in the other direction.
    #
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

# What if we turn things upside down and then do normal subtraction?

def _correct_for_light_travel_time3(observer, target):
    #
    # This version uses normal NumPy subtraction with its default
    # broadcasting, which it makes work by always turning the position
    # vectors upside down when it receives them from other parts of
    # Skyfield, then turning them back over when its own computations
    # are done.
    #
    t = observer.t
    ts = t.ts
    whole = t.whole
    tdb_fraction = t.tdb_fraction

    cposition = observer.position.au
    cvelocity = observer.velocity.au_per_d

    cposition = cposition.T
    cvelocity = cvelocity.T

    tposition, tvelocity, gcrs_position, message = target._at(t)
    tposition = tposition.T
    tvelocity = tvelocity.T

    distance = length_of((tposition - cposition).T).T
    light_time0 = 0.0
    for i in range(10):
        light_time = distance / C_AUDAY
        delta = light_time - light_time0
        if abs(max(delta)) < 1e-12:
            break

        # We assume a light travel time of at most a couple of days.  A
        # longer light travel time would best be split into a whole and
        # fraction, for adding to the whole and fraction of TDB.
        t2 = ts.tdb_jd(whole, (tdb_fraction - light_time).T)

        tposition, tvelocity, gcrs_position, message = target._at(t2)
        tposition = tposition.T
        tvelocity = tvelocity.T
        distance = length_of((tposition - cposition).T).T
        light_time0 = light_time
    else:
        raise ValueError('light-travel time failed to converge')
    tposition = tposition - cposition
    tvelocity = tvelocity - cvelocity
    return tposition.T, tvelocity.T, t, light_time

# What if the observer had an extra dimension so it matches too?

def _correct_for_light_travel_time4(observer, target):
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

    print('======', cposition.shape)
    # ?
    print('tdb_fraction before:', tdb_fraction.shape)
    tdb_fraction = tdb_fraction[:,None]
    print('tdb_fraction after:', tdb_fraction.shape)

    tposition, tvelocity, gcrs_position, message = target._at(t)

    distance = length_of(tposition - cposition)
    light_time0 = 0.0
    for i in range(10):
        print('i =', i)
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

_reconcile  # So CI will think we used it, whether use above is commented or not

if __name__ == '__main__':
    main()
