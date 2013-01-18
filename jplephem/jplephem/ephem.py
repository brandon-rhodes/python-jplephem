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

    def compute(self, name, tdb, differentiate=False):
        """Compute the position and optionally velocity of `name` at `tdb`.

        The barycentric dynamical time `tdb` can be either a single
        value, or an array of many moments for which you want the
        computation run.  Run the `names()` method on a given ephemeris
        to learn the values that it will accept for the `name`
        parameter.

        """
        input_was_scalar = not hasattr(tdb, 'shape')
        if input_was_scalar:
            tdb = np.array((tdb,))

        coefficient_sets = self.load_set(name)
        number_of_sets, axis_count, coefficient_count = coefficient_sets.shape

        date_count = len(tdb)

        days_per_set = (self.jomega - self.jalpha) / number_of_sets
        index, offset = divmod(tdb - self.jalpha, days_per_set)
        index = index.astype(int)
        coefficients = np.rollaxis(coefficient_sets[index], 1)

        T = np.empty((coefficient_count, date_count))
        T[0] = 1.0
        T[1] = t1 = 2.0 * offset / days_per_set - 1.0
        twot1 = t1 + t1
        for i in range(2, coefficient_count):
            T[i] = twot1 * T[i-1] - T[i-2]

        if differentiate:
            dT = np.empty_like(T)
            dT[0] = 0.0
            dT[1] = 1.0
            dT[2] = twot1 + twot1
            for i in range(3, coefficient_count):
                dT[i] = twot1 * dT[i-1] - dT[i-2] + T[i-1] + T[i-1]
            dT *= 2.0 / days_per_set

        if differentiate:
            result = np.empty((2 * axis_count, date_count))
            result[:axis_count] = (T.T * coefficients).sum(axis=2)
            result[axis_count:] = (dT.T * coefficients).sum(axis=2)
        else:
            result = (T.T * coefficients).sum(axis=2)

        if input_was_scalar:
            return result[:,0]
        else:
            return result
