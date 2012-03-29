import os
import numpy as np


class Ephemeris(object):
    """Load and make computations with a JPL planetary ephemeris."""

    def __init__(self, module):
        self.dirpath = os.path.dirname(module.__file__)
        self.__dict__.update(dict(np.load(self.path('constants.npy'))))
        self.earth_share = 1.0 / (1.0 + self.EMRAT)
        self.moon_share = self.EMRAT / (1.0 + self.EMRAT)
        self.sets = {}

    def path(self, filename):
        """Compute the path to a particular file in the ephemeris."""
        return os.path.join(self.dirpath, filename)

    def load_set(self, item):
        """Load the polynomial series for object `n`."""
        s = self.sets.get(item)
        if s is None:
            self.sets[item] = s = np.load(self.path('jpl-%s.npy' % item))
        return s

    def compute(self, item, jed):
        """Given int `item` 1-13 compute its polynomials for date `jed`."""

        # Load the polynomial sets for this item.

        sets = self.load_set(item)

        # How many days are covered by each polynomial set?

        interval = (self.jomega - self.jalpha) / sets.shape[0]

        # Select the paritcular polynomial set in which the date `jed`
        # falls, and determine the offset of `jed` into that date range.

        index, toffset = divmod(jed - self.jalpha, interval)
        if index == sets.shape[0]:
            index -= 1
            toffset += interval
        coefficients = sets[index]

        # We make two passes for this set of Chebyshev coefficients,
        # first computing simple values, and then computing derivatives.
        # Each time through we set up a list of `terms` then multiply by
        # the polynomical coefficients provided by JPL.

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
