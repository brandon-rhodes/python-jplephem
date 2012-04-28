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

    def compute_spherical(self, item, jed):
        """Cartesian 6D vector to spherical components.

        Returns r, alpha, delta, rdot, alphadot, deltadot.
        """
        # ref: TPM manual.
        # See http://phn.github.com/pytpm/_downloads/tpm.pdf
        # r = sqrt(x^2 + y^2 + z^2)
        # alpha = arctan(y/x)  # longitude/azimuth
        # delta = arcsin(z/r)  # latitude/elevation
        # rdot = (x*xdot + y*ydot + z*zdot) / r
        # alphadot = (x*xdot - y*ydot) / (r*cos(delta))**2
        # deltadot = (zdot - (rdot*sin(delta))) / (r*cos(delta))
        # Special conditions are handled as follows:
        # If r==0: rdot = xdot and alpha, delta, alphadot, deltadot = 0.
        # If x==0: alpha = -pi/2 if y < 0, 0 if y == 0 and pi/2 if y > 0.
        # If cos(delta) == 0: rdot = zdot/sin(delta) and
        #   if cos(alpha)==0: deltadot = -ydot /(r.sin(delta).cos(alpha))
        #   else: deltadot = -xdot/ (r.sin(delta).cos(alpha))
        CLOSE_TO_ZERO = 1e-15  # 2.2e16 is the lowest?
        abs = np.abs  # overwrites builtin abs
        sqrt = np.sqrt
        cos = np.cos
        sin = np.sin
        atan2 = np.arctan2
        pi = np.pi

        x, y, z, xdot, ydot, zdot = self.compute(item, jed)
        print x, y, z, xdot, ydot, zdot

        r = sqrt(x ** 2 + y ** 2 + z ** 2)

        # Check for potential singularites i.e, division by 0.
        if r <= CLOSE_TO_ZERO:  # r == 0, r is always +ve
            rdot = xdot
            alpha, alphadot, delta, deltadot = 0.0, 0.0, 0.0, 0.0
            return r, alpha, delta, rdot, alphadot, deltadot

        if (abs(x) <= CLOSE_TO_ZERO):
            if abs(y) <= CLOSE_TO_ZERO:
                alpha = 0.0
            elif y < 0:
                alpha = -pi / 2.0
            elif y > 0:
                alpha = pi / 2.0
        else:
            alpha = atan2(y, x)

        delta = atan2(z, sqrt(x ** 2 + y ** 2))

        if abs(cos(delta)) <= CLOSE_TO_ZERO:  # cos(delta) == 0; poles
            rdot = zdot / sin(delta)
            if abs(cos(alpha)) <= CLOSE_TO_ZERO:
                deltadot = -ydot / (r * sin(delta) * sin(alpha))
            else:
                deltadot = -xdot / (r * sin(delta) * cos(alpha))
            return r, alpha, delta, rdot, alphadot, deltadot

        # If we are here then no singularities should occur.
        rdot = (x * xdot + y * ydot + z * zdot) / r
        alphadot = (x * ydot - y * xdot) / (r * cos(delta)) ** 2
        deltadot = (zdot - rdot * sin(delta)) / (r * cos(delta))

        return r, alpha, delta, rdot, alphadot, deltadot
