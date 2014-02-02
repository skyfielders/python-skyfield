from math import sin, cos
import math
import constants
from positionlib import ICRCoordinates

def semimajorAxisToOrbitalPeriod(axis):
    return (axis ** 3) ** 0.5

def orbitalPeriodToSemimajorAxis(period):
    return (period ** 2) ** (1.0 / 3.0)

def convergeEccentricAnomaly(mean_anomaly, eccentricity, precision):
    # calculate the delta
    delta = 10 ** -precision

    # normalize the mean anomaly
    m = mean_anomaly % constants.TAU

    # set up the first guess
    eccentric_anomaly = constants.TAU
    if eccentricity < 0.8:
        eccentric_anomaly = m

    # do the initial test
    test = eccentric_anomaly - eccentricity * sin(m) - m

    # while we're not happy with the result, and we haven't been dawdling too long
    max_iterations = 30
    count = 0
    while ((math.fabs(test) > delta) and (count < max_iterations)):
        # calculate the next guess for an eccentric anomaly
        eccentric_anomaly = (
            eccentric_anomaly - test /
            (1.0 - eccentricity * cos(eccentric_anomaly))
        )

        # try it
        test = eccentric_anomaly - eccentricity * sin(eccentric_anomaly) - m

        # count the runs, so we don't go forever
        count += 1

    # convert to degrees
    return eccentric_anomaly

def calculateMeanAnomaly(L, wb):
    return L - wb

class KeplerianOrbit:
    def __init__(
            self,
            semimajor_axis,
            eccentricity,
            inclination,
            longitude_ascending,
            argument_perihelion,
            mean_anomaly,
            epoch
        ):
        self.semimajor_axis = semimajor_axis
        self.eccentricity = eccentricity
        self.inclination = inclination
        self.longitude_ascending = longitude_ascending
        self.argument_perihelion = argument_perihelion
        self.mean_anomaly = mean_anomaly
        self.epoch = epoch

    def getECLCoordinatesOnJulianDate(self, date):
         # localize the orbital parameters
        a = self.semimajor_axis
        e = self.eccentricity
        I = self.inclination
        Om = self.longitude_ascending
        #n = 0.230605479
        n = 0.230652907

        w = self.argument_perihelion

        M = self.mean_anomaly
        d = date.tdb - self.epoch.tdb

        Om = Om / 360.0 * constants.TAU
        w = w / 360.0 * constants.TAU
        I = I / 360.0 * constants.TAU
        M = M / 360.0 * constants.TAU
        n = n / 360.0 * constants.TAU

        M += d * n

        # calculate the mean anomaly in rads
        E = convergeEccentricAnomaly(M, e, 30)

        # calculate the initial primes
        x_prime = a * (cos(E) - e)
        y_prime = a * (1 - e ** 2.0) ** (0.5) * sin(E)

        """
        http://ssd.jpl.nasa.gov/txt/aprx_pos_planets.pdf
        x_ecl = cos(w)cos(Om)-sin(w)sin(Om)cos(I) * x_prime +
            (-sin(w)cos(Om) - cos(w)sin(Om)cos(I)) * y_prime
        y_ecl = cos(w)sin(Om)-sin(w)cos(Om)cos(I) * x_prime +
            (-sin(w)cos(Om) - cos(w)sin(Om)cos(I)) * y_prime
        z_ecl = (sin(w)sin(I)) * x_prime +
            (cos(w)sin(I)) * y_prime
        """

        # calculate the ecliptic coordinates
        x_ecl = ((cos(w) * cos(Om) - sin(w) * sin(Om) * cos(I)) * x_prime +
                (-1 * sin(w) * cos(Om) - cos(w) * sin(Om) * cos(I)) * y_prime)
        y_ecl = ((cos(w) * sin(Om) + sin(w) * cos(Om) * cos(I)) * x_prime +
                (-1 * sin(w) * sin(Om) + cos(w) * cos(Om) * cos(I)) * y_prime)
        z_ecl = ((sin(w) * sin(I)) * x_prime + (cos(w) * sin(I)) * y_prime)

        return ICRCoordinates(x_ecl, y_ecl, z_ecl)

    def getICRSCoordinatesOnJulianDate(self, date):
        # J2000 obliquity
        e = 23.43928 * math.pi / 180.0

        # get the ecliptic coords
        ecliptic = self.getECLCoordinatesonJulianDate(date);

        # calculate the equatorial (ICRS) coordinates
        x_eq = ecliptic.x;
        y_eq = cos(e) * ecliptic.y - sin(e) * ecliptic.z
        z_eq = sin(e) * ecliptic.y + cos(e) * ecliptic.z

