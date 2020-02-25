"""Routines for interpreting data from the IAU Minor Planet Center."""

import pandas as pd

MPCORB_URL = 'https://www.minorplanetcenter.net/iau/MPCORB/MPCORB.DAT.gz'

_MPCORB_COLUMNS = [
    ('designation', (0, 7)),
    ('magnitude_H', (8, 13)),
    ('magnitude_G', (14, 19)),
    ('epoch_packed', (20, 25)),
    ('mean_anomaly_degrees', (26, 35)),
    ('argument_of_perihelion_degrees', (37, 46)),
    ('longitude_of_ascending_node_degrees', (48, 57)),
    ('inclination_degrees', (59, 68)),
    ('eccentricity', (70, 79)),
    ('mean_daily_motion_degrees', (80, 91)),
    ('semimajor_axis_au', (92, 103)),
    ('uncertainty', (105, 106)),
    ('reference', (107, 116)),
    ('observations', (117, 122)),
    ('oppositions', (123, 126)),
    ('observation_period', (127, 136)),
    ('rms_residual_arcseconds', (137, 141)),
    ('coarse_perturbers', (142, 145)),
    ('precise_perturbers', (146, 149)),
    ('computer_name', (150, 160)),
]
_MPCORB_NECESSARY_COLUMNS = {
    'designation', 'epoch_packed', 'mean_anomaly_degrees',
    'argument_of_perihelion_degrees', 'longitude_of_ascending_node_degrees',
    'inclination_degrees', 'eccentricity', 'mean_daily_motion_degrees',
    'semimajor_axis_au',
}
_MPCORB_DTYPES = {
    'epoch_packed': 'category',
    'uncertainty': 'category',
    'coarse_perturbers': 'category',
    'precise_perturbers': 'category',
    'computer_name': 'category',
}

def load_mpcorb_dataframe(fobj, full=False):
    """Parse a Minor Planet Center orbits file into a Pandas dataframe.

    The MPCORB file format is documented at:
    https://minorplanetcenter.net/iau/info/MPOrbitFormat.html

    """
    columns = _MPCORB_COLUMNS
    if not full:
        columns = [tup for tup in columns if tup[0] in _MPCORB_NECESSARY_COLUMNS]
    names, colspecs = zip(*columns)
    df = pd.read_fwf(fobj, colspecs, names=names, dtypes=_MPCORB_DTYPES,
                     skiprows=43, compression='gzip')
    return df

def _mpc_comets(self, url, reload=False, filename=None):
    # TODO: switch from the expensive JSON format to parsing "CometEls.txt".
    with self.open(url, reload=reload, filename=filename) as gzip_file:
        with gzip.open(gzip_file) as json_file:
            df = pd.read_json(json_file)

    df.Orbit_type = df.Orbit_type.astype('category')

    return df
