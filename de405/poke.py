# First experiment in poking at the DE405 data.

import numpy as np

AU = 149597870.691  # TODO: read from file
EMRAT = 0.813005600000000044E+02 # TODO: read from file

body_names = (None, 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter',
              'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Moon', 'Sun',
              'Solar System Barycenter', 'Earth-Moon Barycenter',
              'Nutations', 'Librations')
coordinate_names = ('x', 'y', 'z', 'xdot', 'ydot', 'zdot')

#

class Ephemeris(object):

    def compute(self, planet, jed):
        constants = dict(np.load('constants.npy'))
        ja, jz, jd = constants['jalpha'], constants['jomega'], constants['jdelta']

        series = np.load('series%02d.npy' % planet)
        step = (jz - ja) / series.shape[0]  # TODO: isn't this in header file?

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
        r = pleph(ephemeris, jed, target, center)
        delta = r[coordinate_number - 1] - coordinate
        print '%.15f %.15f %.15f' % (
            r[coordinate_number - 1], coordinate, delta,
            )
        if abs(delta) >= 1e-13:
            print 'WARNING: difference =', delta
            break

def compute(ephemeris, jed, target):
    c = ephemeris.compute
    if target == 3:
        return c(3, jed) - c(10, jed) / (1.0 + EMRAT)
    if target <= 2 or 4 <= target <= 9:
        return c(target, jed)
    raise ValueError('hmm %d' % target)

def pleph(ephemeris, jed, target, center):
    # todo: nutations
    # todo: librations

    # Map between test-file target integers and ephemeris series.

    tpos = compute(ephemeris, jed, target)
    cpos = compute(ephemeris, jed, center)
    return (tpos - cpos) / AU

if __name__ == '__main__':
    main()
