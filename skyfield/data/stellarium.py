"""Parse Stellarium data files."""

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
    integer Hipparcos catalog number.

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
