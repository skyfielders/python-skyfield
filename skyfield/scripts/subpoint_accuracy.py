"""How accurate is Skyfield's subpoint computation routine?

Let's call it with a varying number of iterations, and see how
accurately it can turn a Topos-generated position back into a latitude,
longitude, and elevation.

"""
from numpy import einsum
from skyfield.api import Topos, load
from skyfield.constants import AU_M, DEG2RAD
from skyfield.earthlib import reverse_terra

def main():
    ts = load.timescale()
    t = ts.tt(2018, 1, 22, 9, 9, 20)
    trial_angles = 10, 20, 30, 40  # the error peaks around 20 degrees
    trial_elevations = 0, 6000, AU_M
    print(__doc__)
    for n in range(1, 5):
        print('=== {} iterations ==='.format(n))
        print('')
        for elevation_m in trial_elevations:
            for degrees in trial_angles:
                top = Topos(latitude_degrees=degrees, longitude_degrees=123,
                            elevation_m=elevation_m)
                xyz_au = top.at(t).position.au
                xyz_au = einsum('ij...,j...->i...', t.M, xyz_au)
                lat, lon, elev = reverse_terra(xyz_au, t.gast, n)
                lat = lat / DEG2RAD
                error_mas = 60.0 * 60.0 * 1000.0 * abs(degrees - lat)
                print('latitude {} degrees, elevation {} m'
                      ' -> error of {:.2f} mas'
                      .format(degrees, elevation_m, error_mas))
        print('')
    print("""\
Given that iterations=3 pushes the maximum error from tens of mas ("mas"
means " milli-arcsecond") down to hundredths of a mas, it is the value
we have chosen as a default.  A fourth iteration, if we ever chose to
perform one, pushes the error down to "0.00 mas".
""")

if __name__ == '__main__':
    main()
