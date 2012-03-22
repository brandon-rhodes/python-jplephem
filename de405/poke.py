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
        self.sets = [None] * 14

    def load_sets(self, n):
        s = self.sets[n]
        if s is None:
            self.sets[n] = s = np.load('jpl-%02d.npy' % n)
        return s

    def compute(self, item, jed):
        """Given int `item` 1-13 compute its polynomials for date `jed`."""

        # Load the polynomial sets for this item.

        sets = self.load_sets(item)

        # How many days are covered by each polynomial set?

        interval = (self.jomega - self.jalpha) / sets.shape[0]

        # Select the paritcular polynomial set in which the date `jed`
        # falls, and determine the offset of `jed` into that date range.

        index, toffset = divmod(jed - self.jalpha, interval)
        coefficients = sets[index]

        # We make two passes for this set of coefficients, first
        # computing simple values, and then computing derivatives.  Each
        # time through we set up a list of `terms` then multiply by the
        # polynomical coefficients provided by JPL.

        length = sets.shape[2]
        pc = np.zeros(length)  # "p" = position
        vc = np.zeros(length)  # "v" = velocity

        pc[0] = 1.0
        pc[1] = t1 = 2.0 * toffset / interval - 1.0
        twot1 = t1 + t1
        for i in range(2, length):
            pc[i] = twot1 * pc[i-1] - pc[i-2]

        vc[1] = 1.0
        vc[2] = twot1 + twot1
        for i in range(3, length):
            vc[i] = twot1 * vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2]

        position = np.sum(coefficients * pc, axis=1)
        velocity = np.sum(coefficients * vc, axis=1) * (2.0 / interval)

        return np.concatenate((position, velocity))

#

def main():
    ephemeris = Ephemeris()
    verbose = False

    testpo = open('ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/testpo.405')
    lines = iter(testpo)
    while next(lines).strip() != 'EOT':
        continue
    successes = 0
    for line in lines:
        fields = line.split()
        jed = float(fields[2])
        target = int(fields[3])
        center = int(fields[4])
        coordinate_number = int(fields[5])
        coordinate = float(fields[6])
        if verbose:
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
        if verbose:
            print '%.15f %.15f %.15f' % (
                r[coordinate_number - 1], coordinate, delta,
                )
        if abs(delta) >= 1e-13:
            print 'WARNING: difference =', delta
            break
        successes += 1
    print '%d tests successful' % successes

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
    return c(target, jed)

if __name__ == '__main__':
    main()
