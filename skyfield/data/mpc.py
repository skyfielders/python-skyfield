"""Routines for interpreting data from the IAU Minor Planet Center."""

import pandas as pd

mpcorb_url = 'https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz'
_MPCORB_COLUMNS = [
    (0, 7),
    (8, 13),
    (14, 19),
    (20, 25),
    (26, 35),
    (37, 46),
    (48, 57),
    (59, 68),
    (70, 79),
    (80, 91),
    (92, 103),
    (105, 106),
    (107, 116),
    (117, 122),
    (123, 126),
    (127, 136),
    (137, 141),
    (142, 145),
    (146, 149),
    (150, 160),
]
_MPCORB_COLUMN_NAMES = [
    'designation',
    'H',
    'G',
    'epoch',
    'M',
    'Argument',
    'Node',
    'inclination',
    'e',
    'n',
    'a',
    'uncertainty',
    'reference',
    'observation_range',
    'oppositions',
    'observations',
    'rms_residual',
    'coarse_perturbers',
    'precise_perturbers',
    'computer_name',
]
_MPCORB_DTYPES = {
    'epoch': 'category',
    'uncertainty': 'category',
    'coarse_perturbers': 'category',
    'precise_perturbers': 'category',
    'computer_name': 'category',
}

def load_mpcorb_dataframe(fobj):
    """Parse a Minor Planet Center orbits file into a Pandas dataframe.

    The MPCORB file format is documented at:
    https://minorplanetcenter.net/iau/info/MPOrbitFormat.html

    """
    df = pd.read_fwf(fobj, _MPCORB_COLUMNS, names=_MPCORB_COLUMN_NAMES,
                     dtype=_MPCORB_DTYPES, skiprows=43, compression='gzip')
    return df

def mpc_comets(self, url, reload=False, filename=None):
    with self.open(url, reload=reload, filename=filename) as gzip_file:
        with gzip.open(gzip_file) as json_file:
            df = pd.read_json(json_file)

    df.Orbit_type = df.Orbit_type.astype('category')

    return df
