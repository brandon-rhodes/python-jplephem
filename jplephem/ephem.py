"""Compute positions from an ephemeris installed as a Python package.

Note: This entire module is DEPRECATED.  The idea of distributing JPL
ephemerides as Python packages proved to be impractical (they were much
too large for the Python Package Index to easily store and distribute),
and it forced `jplephem` users get their ephemerides from a different
source than mainline astronomers, who use SPICE files.  This package's
documentation now recommends avoiding this old code, and using the
features now build directly into the `SPK` class instead.

"""
import os
import numpy as np


class DateError(ValueError):
    """Date input is outside the range covered by the ephemeris."""


class Ephemeris(object):
    """[DEPRECATED] JPL planetary ephemeris for computing positions on dates."""

    def __init__(self, module):
        self.name = module.__name__.upper()
        self.dirpath = os.path.dirname(module.__file__)
        self.names = tuple(sorted(
            name.split('-')[-1].split('.')[0]
            for name in os.listdir(self.dirpath)
            if not name.startswith('constants') and name.endswith('.npy')
            ))
        path = self.path('constants.npy')
        self.__dict__.update((k.decode('ascii'), v) for k, v in np.load(path))
        self.earth_share = 1.0 / (1.0 + self.EMRAT)
        self.moon_share = self.EMRAT / (1.0 + self.EMRAT)
        self.sets = {}

    def path(self, filename):
        """[DEPRECATED] Compute the path to a particular file in the ephemeris."""
        return os.path.join(self.dirpath, filename)

    def load(self, name):
        """[DEPRECATED] Load the polynomial series for `name` and return it."""
        s = self.sets.get(name)
        if s is None:
            self.sets[name] = s = np.load(self.path('jpl-%s.npy' % name))
        return s

    def position(self, name, tdb, tdb2=0.0):
        """[DEPRECATED] Compute the position of `name` at time ``tdb [+ tdb2]``.

        The position is returned as a NumPy array ``[x y z]``.

        The barycentric dynamical time `tdb` argument should be a float.
        If there are many dates you want computed, then make `tdb` an
        array, which is more efficient than calling this method multiple
        times; the return value will be a two-dimensional array giving a
        row of values for each coordinate.

        For extra precision, the time can be split into two floats; a
        popular choice is to use `tdb` for the integer or half-integer
        date, and `tdb2` to hold the remaining fraction.

        Consult the `names` attribute of this ephemeris for the values
        of `name` it supports, such as ``'mars'`` or ``'earthmoon'``.

        """
        bundle = self.compute_bundle(name, tdb, tdb2)
        return self.position_from_bundle(bundle)

    def position_and_velocity(self, name, tdb, tdb2=0.0):
        """[DEPRECATED] Compute the position and velocity of `name` at ``tdb [+ tdb2]``.

        The position and velocity are returned in a 2-tuple::

            ([x y z], [xdot ydot zdot])

        The barycentric dynamical time `tdb` argument should be a float.
        If there are many dates you want computed, then make `tdb` an
        array, which is more efficient than calling this method multiple
        times; the return values will be two-dimensional arrays giving a
        row of values for each coordinate.

        For extra precision, the time can be split into two floats; a
        popular choice is to use `tdb` for the integer or half-integer
        date, and `tdb2` to hold the remaining fraction.

        Consult the `names` attribute of this ephemeris for the values
        of `name` it supports, such as ``'mars'`` or ``'earthmoon'``.

        """
        bundle = self.compute_bundle(name, tdb, tdb2)
        position = self.position_from_bundle(bundle)
        velocity = self.velocity_from_bundle(bundle)
        return position, velocity

    def compute(self, name, tdb):
        """[DEPRECATED] Legacy routine that concatenates position and velocity vectors.

        This routine is deprecated.  Use the methods `position()` and
        `position_and_velocity()` instead.  This method follows the same
        calling convention, but incurs extra copy operations in order to
        return a single NumPy array::

            [x y z xdot ydot zdot]

        """
        bundle = self.compute_bundle(name, tdb, 0.0)
        position = self.position_from_bundle(bundle)
        velocity = self.velocity_from_bundle(bundle)
        return np.concatenate((position, velocity))

    def compute_bundle(self, name, tdb, tdb2=0.0):
        """[DEPRECATED] Return a tuple of coefficients and parameters for `tdb`.

        The return value is a tuple that bundles together the
        coefficients and other Chebyshev intermediate values that are
        needed for the computation of either the position or velocity.
        The bundle can then be passed to either `position_from_bundle()`
        or `velocity_from_bundle()` to finish the computation.  See the
        package-level documentation for details; most users will simply
        call `position()` or `position_and_velocity()` instead.

        The barycentric dynamical time `tdb` argument should be a float.
        If there are many dates you want computed, then make `tdb` an
        array, which is more efficient than calling this method multiple
        times; the return values will be arrays providing a value for
        each time in `tdb`.

        For extra precision, the time can be split into two floats; a
        popular choice is to use `tdb` for the integer or half-integer
        date, and `tdb2` to hold the remaining fraction.

        Consult the `names` attribute of this ephemeris for the values
        of `name` it supports, such as ``'mars'`` or ``'earthmoon'``.

        """
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

        bundle = coefficients, days_per_set, T, twot1
        return bundle

    def position_from_bundle(self, bundle):
        """[DEPRECATED] Return position, given the `coefficient_bundle()` return value."""

        coefficients, days_per_set, T, twot1 = bundle
        return (T.T * coefficients).sum(axis=2)

    def velocity_from_bundle(self, bundle):
        """[DEPRECATED] Return velocity, given the `coefficient_bundle()` return value."""

        coefficients, days_per_set, T, twot1 = bundle
        coefficient_count = coefficients.shape[2]

        # Chebyshev derivative:

        dT = np.empty_like(T)
        dT[0] = 0.0
        dT[1] = 1.0
        dT[2] = twot1 + twot1
        for i in range(3, coefficient_count):
            dT[i] = twot1 * dT[i-1] - dT[i-2] + T[i-1] + T[i-1]
        dT *= 2.0
        dT /= days_per_set

        return (dT.T * coefficients).sum(axis=2)
