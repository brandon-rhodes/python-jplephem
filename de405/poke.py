# First experiment in poking at the DE405 data.

import numpy as np

body_names = (None, 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter',
              'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Moon', 'Sun',
              'Solar System Barycenter', 'Earth-Moon Barycenter',
              'Nutations', 'Librations')
coordinate_names = ('x', 'y', 'z', 'xdot', 'ydot', 'zdot')

#

class Ephemeris(object):

    def __init__(self):
        # Load constants as instance attributes.
        self.__dict__.update(dict(np.load('constants.npy')))  # Ruby
        self.earth_share = 1.0 / (1.0 + self.EMRAT)
        self.moon_share = self.EMRAT / (1.0 + self.EMRAT)
        self.series = [None] * 14

    def load_series(self, n):
        s = self.series[n]
        if s is None:
            self.series[n] = s = np.load('series%02d.npy' % n)
        return s

    def compute(self, planet, jed):
        ja, jz, jd = self.jalpha, self.jomega, self.jdelta

        series = self.load_series(planet)
        step = (jz - ja) / series.shape[0]

        l, jremain = divmod(jed - ja, step)
        tc = 2.0 * jremain / step - 1.0
        ncf = series.shape[2]

        pc = np.zeros(ncf)
        pc[0] = 1.0
        pc[1] = tc
        twot = tc + tc
        for i in range(2, ncf):
            pc[i] = twot * pc[i-1] - pc[i-2]

        coords = np.sum(series[l] * pc, axis=1)
        # can return coords at this point if velocities not needed

        vfac = 2.0 / step

        vc = np.zeros(ncf)
        vc[1] = 1.0
        vc[2] = twot + twot
        for i in range(3, ncf):
            vc[i] = twot * vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2]

        velocities = np.sum(series[l] * vc, axis=1) * vfac
        return np.concatenate((coords, velocities))

#

def main():
    ephemeris = Ephemeris()

    testpo = open('ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/testpo.405')
    lines = iter(testpo)
    while next(lines).strip() != 'EOT':
        continue
    for line in lines:
        fields = line.split()
        jed = float(fields[2])
        target = int(fields[3])
        center = int(fields[4])
        coordinate_number = int(fields[5])
        coordinate = float(fields[6])
        print '%s %s %s(%d) -> %s(%d) field #%d' % (
            fields[1], jed, body_names[center], center,
            body_names[target], target, coordinate_number)
        if target == 14:
            r = ephemeris.compute(12, jed)
        elif target == 15:
            r = ephemeris.compute(13, jed)
        else:
            tpos = compute(ephemeris, jed, target)
            cpos = compute(ephemeris, jed, center)
            r = (tpos - cpos) / ephemeris.AU

        delta = r[coordinate_number - 1] - coordinate
        print '%.15f %.15f %.15f' % (
            r[coordinate_number - 1], coordinate, delta,
            )
        if abs(delta) >= 1e-13:
            print 'WARNING: difference =', delta
            break

def compute(ephemeris, jed, target):
    if target == 12:
        return np.zeros(6)  # solar system barycenter is our origin
    c = ephemeris.compute
    if target == 13:
        return c(3, jed)
    if target == 3:
        return c(3, jed) - c(10, jed) * ephemeris.earth_share
    if target == 10:
        return c(3, jed) + c(10, jed) * ephemeris.moon_share
    if target <= 11:
        return c(target, jed)
    raise ValueError('hmm %d' % target)

if __name__ == '__main__':
    main()
