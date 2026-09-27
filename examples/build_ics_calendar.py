from datetime import datetime
from skyfield import almanac
from skyfield.api import N, W, load, wgs84

def main():
    entries = []
    append = entries.append
    extend = entries.extend

    ts = load.timescale()
    t0 = ts.utc(2026, 1, 1)
    t1 = ts.utc(2028, 1, 1)
    eph = load('de421.bsp')
    sun = eph['sun']
    earth = eph['earth']
    observatory = earth + wgs84.latlon(35.2029, -111.6646)

    def find(f, *messages):
        t, y = almanac.find_discrete(t0, t1, f)
        for ti, yi in zip(t.utc_datetime(), y):
            message = messages[yi]
            yield format_ics_entry(
                message,
                ti,
                ti,
            )

    extend(find(
        almanac.seasons(eph),
        'March Equinox',
        'June Solstice',
        'September Equinox',
        'December Solstice',
    ))

    extend(find(
        almanac.moon_phases(eph),
        'New Moon',
        'Moon at First Quarter',
        'Full Moon',
        'Moon at Last Quarter',
    ))

    extend(find(
        almanac.moon_nodes(eph),
        'Moon at descending node',
        'Moon at ascending node',
    ))

    for name in ['Mars', 'Jupiter', 'Saturn']:
        planet = eph[name + ' Barycenter']
        extend(find(
            almanac.oppositions_conjunctions(eph, planet),
            '{} superior conjunction'.format(name),
            '{} at opposition'.format(name),
        ))

    s = ICS_WRAPPER.format(''.join(entries))
    print(s.replace('\n', '\r\n'), end='')

def format_ics_entry(summary, start_dt, end_dt):
    words = summary.lower().split()
    words.append(start_dt.strftime('%Y-%m-%d'))
    uid = '-'.join(words)

    fmt = "%Y%m%dT%H%M%SZ"
    dtstamp = datetime.utcnow().strftime(fmt)
    dtstart = start_dt.strftime(fmt)
    dtend = end_dt.strftime(fmt)

    ics_text = f"""BEGIN:VEVENT
UID:{uid}
DTSTAMP:{dtstamp}
DTSTART:{dtstart}
DTEND:{dtend}
SUMMARY:{summary}
END:VEVENT
"""

    return ics_text

ICS_WRAPPER = """BEGIN:VCALENDAR
VERSION:2.0
PRODID:-//Skyfield calendar example//EN
{}END:VCALENDAR
"""

if __name__ == '__main__':
    main()
