from skyfield.data import text_pck

TEXT = b"""JPL/PCK
This is ignored.
\\begindata
A = 3
B = (4)
C = ( 5 'Six' )
D_+=( 7,8 )
E = (9 1.0000D1)
E += (11)
E += 12
\\begintext
This is also ignored.
E += (13)
"""

def test_loading():
    lines = TEXT.splitlines(True)
    variables = {}
    text_pck.load(lines, variables)
    assert variables == {
        'A': 3,
        'B': 4,
        'C': [5, 'Six'],
        'D_': [7, 8],
        'E': [9, 10.0, 11, 12],
    }

def test_parsing():
    lines = TEXT.splitlines(True)
    assignments = list(text_pck.parse(lines))
    assert assignments == [
        ('A', b'=', [3]),
        ('B', b'=', [4]),
        ('C', b'=', [5, 'Six']),
        ('D_', b'+=', [7, 8]),
        ('E', b'=', [9, 10.0]),
        ('E', b'+=', [11]),
        ('E', b'+=', [12]),
    ]
