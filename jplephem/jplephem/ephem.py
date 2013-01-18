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
        # There are two main expenses involved in this computation.
        #
        # The first expense is fetching each coefficient set within
        # whose period at least one `tdb` value lies.  The most naive
        # implementation would do one lookup per `tdb`, either forcing
        # the creation of one view per `tdb` if we use a simple index
        # coefficients[n] while looping over `tdb`, or an actual data
        # copy per `tdb` if we try fancy indexing with coefficients[(n0,
        # n1, n2...)].  One lookup per `tdb` is, of course, inevitable
        # if `tdb` values are spaced far enough apart that each of them
        # falls under different coefficients; but since coefficients
        # tend to cover a couple of weeks, a few adjacent values of
        # `tdb` can often re-use the same coefficients (assuming an
        # ordered array `tdb`; we do not attempt to optimize for
        # unordered input).
        #
        # The second expense is computing the Chebyshev polynomial for
        # each `jremainder` by which the corresponding `tdb` exceeds the
        # start date for its polynomial.  Right now we do not attempt to
        # optimize this, but in the future we might try to detect if
        # successive runs of `tdb` fall at the same remainders, which
        # happens when the difference between uniformly spaced `tdb`
        # values is evenly divisible by the period covered by a
        # polynomial.

        coefficient_sets = self.load_set(name)
        number_of_sets, axis_count, coefficient_count = coefficient_sets.shape
        days_per_set = (self.jomega - self.jalpha) / number_of_sets

        input_was_scalar = not hasattr(tdb, 'shape')
        if input_was_scalar:
            tdb = np.array((tdb,))
        date_count = len(tdb)

        index, offset = divmod(tdb - self.jalpha, days_per_set)

        unique_index, put_index = np.unique(index, return_inverse=True)
        unique_offset, put_offset = np.unique(offset, return_inverse=True)

        unique_index = unique_index.astype(int)
        coefficients = [coefficient_sets[i] for i in unique_index]

        T = np.empty((coefficient_count, unique_offset.size))
        T[0] = 1.0
        T[1] = t1 = 2.0 * unique_offset / days_per_set - 1.0
        twot1 = t1 + t1
        for i in range(2, coefficient_count):
            T[i] = twot1 * T[i-1] - T[i-2]

        p_arrays = [T[:,i] for i in range(T.shape[1])]

        if differentiate:
            dT = np.zeros((coefficient_count, unique_offset.size))
            dT[1] = 1.0
            dT[2] = twot1 + twot1
            for i in range(3, coefficient_count):
                dT[i] = twot1 * dT[i-1] - dT[i-2] + T[i-1] + T[i-1]
            dT *= 2.0 / days_per_set

            v_arrays = [dT[:,i] for i in range(dT.shape[1])]

            result = np.empty((axis_count * 2, date_count))
        else:
            result = np.empty((axis_count, date_count))

        scratch = np.empty((axis_count, coefficient_count))

        for i in range(date_count):
            i1 = put_index[i]
            i2 = put_offset[i]

            c = coefficients[i1]
            scratch[:,:] = c
            scratch *= p_arrays[i2]
            result[:axis_count, i] = scratch.sum(axis=1)

            if differentiate:
                scratch[:,:] = c
                scratch *= v_arrays[i2]
                result[axis_count:,i] = scratch.sum(axis=1)

        if input_was_scalar:
            return result[:,0]
        else:
            return result
