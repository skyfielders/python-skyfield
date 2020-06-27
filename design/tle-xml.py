#!/usr/bin/env python3
from skyfield.api import load, EarthSatellite
from skyfield.timelib import julian_date
from skyfield.constants import AU_KM, AU_M
from sgp4.api import Satrec, WGS72

from math import pi, radians
import datetime
import xml.etree.ElementTree as ET

# Test for equivalence between TLE and XML data for the same epoch
#
# Load TLE  (local file "stations-tle.txt")
# Load XML  (local file "stations.xml")
# For each EarthSatellite in TLE
#     Locate in XML, construct 'twin' using from_satrec
#     Check for similar position, velocity at same time
#
# stations.xml is from:
#    https://celestrak.com/NORAD/elements/gp.php?GROUP=STATIONS&FORMAT=XML
# The XML structure is described in CCSDS 502.0-B-2 Orbit Data Messages
#    https://public.ccsds.org/Pubs/502x0b2c1.pdf
#
# stations-tle.txt is from:
#    https://celestrak.com/NORAD/elements/stations.txt
#
def satrec_from_XML(root, satellite_id):
    #iss_segment = root.find(f".//*[OBJECT_NAME='{satname}']/..")
    iss_segment = root.find(f".//*[NORAD_CAT_ID='{satellite_id}']/../..") # find the matching "segment" node
    meta = iss_segment.find("metadata")
    name = meta.find("OBJECT_NAME").text
    object_id = meta.find("OBJECT_ID").text
    ref_frame = meta.find("REF_FRAME").text
    mean_element_theory = meta.find("MEAN_ELEMENT_THEORY").text
    assert ref_frame == "TEME"
    assert mean_element_theory == "SGP4"
    # "data" node contains meanElements and tleParameters
    data = iss_segment.find("data")
    mean_elements = data.find("meanElements")
    tle_parameters = data.find("tleParameters")
    # Extract a subset of fields to create a Satrec
    epoch = mean_elements.find("EPOCH").text
    norad_cat_id = tle_parameters.find("NORAD_CAT_ID").text
    bstar = tle_parameters.find("BSTAR").text
    mm_dot = tle_parameters.find("MEAN_MOTION_DOT").text
    mm_ddot = tle_parameters.find("MEAN_MOTION_DDOT").text
    eccentricity = mean_elements.find("ECCENTRICITY").text
    arg_pericenter = mean_elements.find("ARG_OF_PERICENTER").text
    inclination = mean_elements.find("INCLINATION").text
    mean_anomaly = mean_elements.find("MEAN_ANOMALY").text
    mean_motion = mean_elements.find("MEAN_MOTION").text
    ra_asc_node = mean_elements.find("RA_OF_ASC_NODE").text
    # Per SGP calculate epoch from 1949 Dec 31
    base_epoch_julian = julian_date(1949, 12, 31)
    # Format of XML epoch is described in 6.5.9 of CCSDS Orbit Data Messages
    # (this does not cover all valid variations listed in the spec)
    dt = datetime.datetime.strptime(epoch, "%Y-%m-%dT%H:%M:%S.%f")
    epoch_julian = julian_date(dt.year, dt.month, dt.day, dt.hour, dt.minute,
                    dt.second + (dt.microsecond / 1E6))
    epoch_since_base = epoch_julian - base_epoch_julian
    # When MEAN_ELEMENT_THEORY is SGP4, convert revolutions per day -> radians/minute
    mean_motion_radians_per_minute = float(mean_motion) * 2 * pi / MINUTES_PER_DAY
    #
    # According to https://pypi.org/project/sgp4/
    # ndot and nddot are ignored by SGP4, so no further conversions are attempted here
    # sgp4.io processes them further, eg.
    # https://github.com/brandon-rhodes/python-sgp4/blob/4c13506a410f7fba17dbf38bd722a2ddb322f197/sgp4/io.py#L195
    #
    # "Build a satellite from orbital elements"
    # https://rhodesmill.org/skyfield/earth-satellites.html#from-satrec
    satrec = Satrec()
    satrec.sgp4init(
        WGS72,                  # gravity model
        "i",                    # 'a' = old AFSPC mode, 'i' = improved mode
        int(norad_cat_id),      # satnum: Satellite number
        epoch_since_base,       # epoch: days since 1949 December 31 00:00 UT
        float(bstar),           # bstar: drag coefficient (/earth radii)
        float(mm_dot),          # ndot: ballistic coefficient (revs/day)
        float(mm_ddot),         # nddot: second derivative of mean motion (revs/day^3)
        float(eccentricity),    # ecco: eccentricity
        radians(float(arg_pericenter)), # argpo: argument of perigee (radians)
        radians(float(inclination)),    # inclo: inclination (radians)
        radians(float(mean_anomaly)),   # mo: mean anomaly (radians)
        mean_motion_radians_per_minute, # no_kozai: mean motion (radians/minute)
        radians(float(ra_asc_node)),    # nodeo: right ascension of ascending node (radians)
    )
    return satrec

MINUTES_PER_DAY = 24 * 60
ONE_METER = 1.0 / AU_M
ONE_KM_PER_HOUR = 1.0 * 24.0 / AU_KM
ts = load.timescale(builtin=True)
# TLE
stations_tle = load.tle_file("stations-tle.txt", reload=False)
# XML
tree = ET.parse("stations.xml")
root = tree.getroot()

t0 = ts.utc(2020, 6, 27)
for sat_tle in stations_tle:
    name = sat_tle.name.strip()
    print(name, sat_tle.model.satnum)
    # Find and parse XML twin
    satrec = satrec_from_XML(root, sat_tle.model.satnum)
    sat_xml = EarthSatellite.from_satrec(satrec, ts)
    p_xml = sat_xml.at(t0)
    p_tle = sat_tle.at(t0)
    #print("XML Location", p_xml.position.au)
    #print("TLE location", p_tle.position.au)
    #print()
    #print("XML velocity", p_xml.velocity.au_per_d)
    #print("TLE velocity", p_tle.velocity.au_per_d)
    assert abs(p_xml.position.au - p_tle.position.au).max() < ONE_METER
    assert abs(p_xml.velocity.au_per_d - p_tle.velocity.au_per_d).max() < ONE_KM_PER_HOUR
