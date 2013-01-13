import os
import numpy as np


class Ephemeris(object):
    """Load and make computations with a JPL planetary ephemeris."""

    def __init__(self, module):
        self.dirpath = os.path.dirname(module.__file__)
        self.names = [ name.split('-')[-1].split('.')[0]
                       for name in os.listdir(self.dirpath)
                       if name.endswith('.npy') ]
        path = self.path('constants.npy')
        self.__dict__.update((k.decode('ascii'), v) for k, v in np.load(path))
        self.earth_share = 1.0 / (1.0 + self.EMRAT)
        self.moon_share = self.EMRAT / (1.0 + self.EMRAT)
        self.sets = {}

    def path(self, filename):
        """Compute the path to a particular file in the ephemeris."""
        return os.path.join(self.dirpath, filename)

    def load_set(self, name):
        """Load the polynomial series for `name`."""
        s = self.sets.get(name)
        if s is None:
            self.sets[name] = s = np.load(self.path('jpl-%s.npy' % name))
        return s

    def compute(self, name, jed):
        """Compute the position and velocity of `name` on date `jed`.

        Run the `names()` method on a given ephemeris to learn the
        values that it will accept for the `name` parameter.

        """
        # Load the polynomial sets for this item.

        sets = self.load_set(name)

        # How many days are covered by each polynomial set?

        interval = (self.jomega - self.jalpha) / sets.shape[0]

        # Select the paritcular polynomial set in which the date `jed`
        # falls, and determine the offset of `jed` into that date range.

        index, toffset = divmod(jed - self.jalpha, interval)

        # TODO: restore the ability to ask about the last date in the ephemeris
        # print index, toffset
        # at_end = index >= sets.shape[0]
        # index[at_end] -= 1
        # toffset[at_end] += interval
        # print index, toffset
        # if index == sets.shape[0]:
        #     index -= 1
        #     toffset += interval
        if hasattr(jed, 'shape'):
            index = index.astype(int)
        coefficients = sets[index]

        # We make two passes for this set of Chebyshev coefficients,
        # first computing simple values, and then computing derivatives.
        # Each time through we set up a list of `terms` then multiply by
        # the polynomical coefficients provided by JPL.

        length = sets.shape[2]
        # print length
        pc = [0.0] * length  # position
        vc = [0.0] * length  # velocity

        pc[0] = 1.0
        pc[1] = t1 = 2.0 * toffset / interval - 1.0
        twot1 = t1 + t1
        for i in range(2, length):
            pc[i] = twot1 * pc[i-1] - pc[i-2]

        vc[1] = 1.0
        vc[2] = twot1 + twot1
        for i in range(3, length):
            vc[i] = twot1 * vc[i-1] + pc[i-1] + pc[i-1] - vc[i-2]

        # print 'pc:', pc
        # print 'coeff:', coefficients

        # Generate return values of the same dimension as jed.

        x = jed * 0.0
        y = jed * 0.0
        z = jed * 0.0
        dx = jed * 0.0
        dy = jed * 0.0
        dz = jed * 0.0

        for i in range(length):
            # print 'x:', x
            # print 'pc[i]:', pc[i]
            # print 'coefficients[0, i]:', coefficients[0, i]
            if hasattr(jed, 'shape'):
                x += pc[i] * coefficients[:, 0, i]
                y += pc[i] * coefficients[:, 1, i]
                z += pc[i] * coefficients[:, 2, i]
                dx += vc[i] * coefficients[:, 0, i]
                dy += vc[i] * coefficients[:, 1, i]
                dz += vc[i] * coefficients[:, 2, i]
            else:
                x += pc[i] * coefficients[0, i]
                y += pc[i] * coefficients[1, i]
                z += pc[i] * coefficients[2, i]
                dx += vc[i] * coefficients[0, i]
                dy += vc[i] * coefficients[1, i]
                dz += vc[i] * coefficients[2, i]

        dd = (2.0 / interval)
        dx *= dd
        dy *= dd
        dz *= dd

        return x, y, z, dx, dy, dz
