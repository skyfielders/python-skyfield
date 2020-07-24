"""Parse Stellarium data files."""

from collections import namedtuple

StarName = namedtuple('StarName', 'hip name')

def parse_constellations(lines):
    """Return a list of constellation outlines.

    Each constellation outline is a list of edges, each of which is
    drawn between a pair of specific stars::

        [
            (name, [(star1, star2), (star3, star4), ...]),
            (name, [(star1, star2), (star3, star4), ...]),
            ...
        ]

    Each name is a 3-letter constellation abbreviation; each star is an
    integer Hipparcos catalog number.  See :ref:`neowise-chart` for an
    example of how to combine this data with the Hipparcos star catalog
    to draw constellation lines on a chart.

    """
    constellations = []
    for line in lines:
        line = line.lstrip()
        if line.startswith(b'#'):
            continue
        fields = line.split()
        if not fields:
            continue
        name = fields[0]
        edges = [(int(fields[i]), int(fields[i+1]))
                 for i in range(2, len(fields), 2)]
        constellations.append((name.decode('utf-8'), edges))
    return constellations

def parse_star_names(lines):
    """Return the names in a Stellarium ``star_names.fab`` file.

    Returns a list of named tuples, each of which offers a ``.hip``
    attribute with a Hipparcos catalog number and a ``.name`` attribute
    with the star name.  Do not depend on the tuple having only length
    two; additional fields may be added in the future.

    """
    names = []
    for line in lines:
        line = line.strip()
        if line == b'' or line.startswith(b'#'):
            continue
        fields = line.split()
        hip, name = fields[0].split(b'|')
        names.append(StarName(
            int(hip),
            name.strip(b'_(")').decode('utf-8'),
        ))
    return names
