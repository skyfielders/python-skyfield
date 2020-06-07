# -*- coding: utf-8 -*-
"""Routines specific to some calculation related with cultures in East Asia."""

from skyfield.constants import tau
from skyfield.nutationlib import iau2000b_radians

# Names of solar terms in Simplified Chinese
SOLAR_TERMS_ZHS = [
    '春分',
    '清明',
    '谷雨',
    '立夏',
    '小满',
    '芒种',
    '夏至',
    '小暑',
    '大暑',
    '立秋',
    '处暑',
    '白露',
    '秋分',
    '寒露',
    '霜降',
    '立冬',
    '小雪',
    '大雪',
    '冬至',
    '小寒',
    '大寒',
    '立春',
    '雨水',
    '惊蛰',
]

# Names of solar terms in Traditional Chinese
SOLAR_TERMS_ZHT = [
    '春分',
    '清明',
    '穀雨',
    '立夏',
    '小滿',
    '芒種',
    '夏至',
    '小暑',
    '大暑',
    '立秋',
    '處暑',
    '白露',
    '秋分',
    '寒露',
    '霜降',
    '立冬',
    '小雪',
    '大雪',
    '冬至',
    '小寒',
    '大寒',
    '立春',
    '雨水',
    '驚蟄',
]

# Names of solar terms in Japanese
SOLAR_TERMS_JP = [
    '春分',
    '清明',
    '穀雨',
    '立夏',
    '小満',
    '芒種',
    '夏至',
    '小暑',
    '大暑',
    '立秋',
    '処暑',
    '白露',
    '秋分',
    '寒露',
    '霜降',
    '立冬',
    '小雪',
    '大雪',
    '冬至',
    '小寒',
    '大寒',
    '立春',
    '雨水',
    '啓蟄',
]

# Names of solar terms in Vietnamese
SOLAR_TERMS_VN = [
    'Xuân phân',
    'Thanh minh',
    'Cốc vũ',
    'Lập hạ',
    'Tiểu mãn',
    'Mang chủng',
    'Hạ chí',
    'Tiểu thử',
    'Đại thử',
    'Lập thu',
    'Xử thử',
    'Bạch lộ',
    'Thu phân',
    'Hàn lộ',
    'Sương giáng',
    'Lập đông',
    'Tiểu tuyết',
    'Đại tuyết',
    'Đông chí',
    'Tiểu hàn',
    'Đại hàn',
    'Lập xuân',
    'Vũ thủy',
    'Kinh trập',
]

def solar_terms(ephemeris):
    """Build a function of time that returns the solar terms of the year.

    The function that this returns will expect a single argument that is
    a :class:`~skyfield.timelib.Time` and will return 0 through 23 for
    the solar terms.

    The name of the solar terms may vary in different cultures, so we
    have table of names in Simplified Chinese, Traditional Chinese,
    Japanese and Vietnamese

    Reference:

    https://en.wikipedia.org/wiki/Solar_term

    """
    earth = ephemeris['earth']
    sun = ephemeris['sun']

    def solar_term_at(t):
        """Return season 0 through 23 at time `t`."""
        t._nutation_angles_radians = iau2000b_radians(t)
        e = earth.at(t)
        _, slon, _ = e.observe(sun).apparent().ecliptic_latlon('date')
        return (slon.radians // (tau / 24) % 24).astype(int)

    solar_term_at.rough_period = 15.0
    return solar_term_at
