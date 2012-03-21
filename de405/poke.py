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
        ja, jz = constants['JEDA'], constants['JEDZ']

        series = np.load('series%02d.npy' % planet)
        t2 = (jz - ja) / series.shape[0]  # TODO: isn't this in header file?

        l, jremain = divmod(jed - ja, t2)
        t1 = jremain / t2

        ncf = series.shape[2]
        ncm = 3
        tc = 2.0 * t1 - 1.0

        pc = [1.0, tc]
        twot = tc + tc

        for i in range(2, ncf):
            pc.append(twot * pc[-1] - pc[-2])

        answers = []
        for i in range(ncm):
            answers.append(sum(
                    pc[j] * series[l,i,j]
                    for j in reversed(range(ncf))
                    ))

        return np.array(answers)

#

def main():
    ephemeris = Ephemeris()

    testpo = open('ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/testpo.405')
    lines = iter(testpo)
    while next(lines).strip() != 'EOT':
        continue
    for line in lines:
        fields = line.split()
        print fields
        jed = float(fields[2])
        target = int(fields[3])
        center = int(fields[4])
        coordinate_number = int(fields[5])
        coordinate = float(fields[6])
        print jed, body_names[center], '->', body_names[target]
        r = pleph(ephemeris, jed, target, center)
        delta = r[coordinate_number - 1] - coordinate
        print '%.15f %.15f %.15f' % (
            r[coordinate_number - 1],
            coordinate,
            delta,
            )
        if abs(delta) >= 1e-13:
            print 'WARNING: difference =', delta
        break

def pleph(ephemeris, jed, target, center):
    # todo: nutations
    # todo: librations
    tpos = ephemeris.compute(target, jed)
    cpos = ephemeris.compute(center, jed)
    if center == 3:
        moonpos = ephemeris.compute(10, jed)
        cpos -= moonpos / (1.0 + EMRAT)
    return (tpos - cpos) / AU

if __name__ == '__main__':
    main()
