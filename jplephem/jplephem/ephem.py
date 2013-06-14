import os
import numpy as np


class DateError(Exception):
    """Date input is outside the range covered by the ephemeris."""


class Ephemeris(object):
    """Load and make computations with a JPL planetary ephemeris."""

    def __init__(self, module):
        self.name = module.__name__.upper()
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

    def load(self, name):
        """Load the polynomial series for `name`."""
        s = self.sets.get(name)
        if s is None:
            self.sets[name] = s = np.load(self.path('jpl-%s.npy' % name))
        return s

    def position(self, name, tdb, tdb2=0.):
        """Compute the position of `name` at time `tdb` [+`tdb2`].

        Run the `names()` method on this ephemeris to learn the values
        it will accept for the `name` parameter, such as ``'mars'`` and
        ``'earthmoon'``.

        The barycentric dynamical time `tdb` can be either a normal number or
        a NumPy array of times, in which case each of the three return values
        ``(x, y, z)`` will be an array.

        For extra precision, one can give a two-part tdb; rounding errors are
        avoided if `tdb` is an (half-)integer part and `tdb2` a fraction.
        """
        return self._interpolate(name, tdb, tdb2, False)

    def compute(self, name, tdb, tdb2=0.):
        """Compute the position and velocity of `name` at time `tdb` [+`tdb2`].

        Run the `names()` method on this ephemeris to learn the values
        it will accept for the `name` parameter, such as ``'mars'`` and
        ``'earthmoon'``.

        The barycentric dynamical time `tdb` can be either a normal number or
        a NumPy array of times, in which case each of the three return values
        ``(x, y, z)`` will be an array.

        For extra precision, one can give a two-part tdb; rounding errors are
        avoided if `tdb` is an (half-)integer part and `tdb2` a fraction.
        """
        return self._interpolate(name, tdb, tdb2, True)

    def _interpolate(self, name, tdb, tdb2=0., differentiate=True):
        input_was_scalar = getattr(tdb, 'shape', ()) == ()
        if input_was_scalar:
            tdb = np.array((tdb,))
        # no need to deal with tdb2; numpy broadcast will add fine below.

        coefficient_sets = self.load(name)
        number_of_sets, axis_count, coefficient_count = coefficient_sets.shape

        jalpha, jomega = self.jalpha, self.jomega
        days_per_set = (jomega - jalpha) / number_of_sets
        # to keep precision, first subtract, then add
        index, offset = divmod((tdb - jalpha) + tdb2, days_per_set)
        index = index.astype(int)

        if (index < 0).any() or (number_of_sets < index).any():
            raise DateError('ephemeris %s only covers dates %.1f through %.1f'
                            % (self.name, jalpha, jomega))

        omegas = (index == number_of_sets)
        index[omegas] -= 1
        offset[omegas] += days_per_set

        coefficients = np.rollaxis(coefficient_sets[index], 1)

        # Chebyshev recurrence:

        T = np.empty((coefficient_count, len(index)))
        T[0] = 1.0
        T[1] = t1 = 2.0 * offset / days_per_set - 1.0
        twot1 = t1 + t1
        for i in range(2, coefficient_count):
            T[i] = twot1 * T[i-1] - T[i-2]

        if not differentiate:
            result = (T.T * coefficients).sum(axis=2)
        else:
            dT = np.empty_like(T)
            dT[0] = 0.0
            dT[1] = 1.0
            dT[2] = twot1 + twot1
            for i in range(3, coefficient_count):
                dT[i] = twot1 * dT[i-1] - dT[i-2] + T[i-1] + T[i-1]
            dT *= 2.0 / days_per_set

            result = np.empty((2 * axis_count, len(index)))
            result[:axis_count] = (T.T * coefficients).sum(axis=2)
            result[axis_count:] = (dT.T * coefficients).sum(axis=2)

        if input_was_scalar:
            return np.squeeze(result)

        return result
