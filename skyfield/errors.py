"""Exceptions specific to the Skyfield library."""

class DeprecationError(Exception):
    """Explain that a Skyfield feature has been removed."""

class EphemerisRangeError(ValueError):
    """An ephemeris has been asked about positions outside its time range.

    Attributes:

    - `start_time`, `end_time`: the range of times supported by the segment
    - `time_mask`: Boolean array where ``True`` marks out-of-range times
    - `segment`: the ephemeris segment that was asked for positions

    """
    def __init__(self, message, start_time, end_time, time_mask, segment):
        self.args = message,
        self.start_time = start_time
        self.end_time = end_time
        self.time_mask = time_mask
        self.segment = segment
