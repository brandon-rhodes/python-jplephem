"""Tests for ``jplephem``.

See the accompanying ``jpltest`` module for a more intense numerical
test suite that can verify that ``jplephem`` delivers, in a large number
of cases, the same results as when the ephemerides are run at JPL.  This
smaller and more feature-oriented suite can be run with::

    python -m unittest discover jplephem

"""
import numpy as np
from functools import partial
from jplephem import Ephemeris
from unittest import TestCase

class Tests(TestCase):

    def check0(self, x, y, z, dx, dy, dz):
        eq = partial(self.assertAlmostEqual, delta=1.0)
        eq(x, -3.84258043e+07)
        eq(y, 1.29508901e+08)
        eq(z, 5.62536412e+07)
        eq(dx, -2.52115204e+06)
        eq(dy, -6.40613144e+05)
        eq(dz, -2.78642593e+05)

    def check1(self, x, y, z, dx, dy, dz):
        eq = partial(self.assertAlmostEqual, delta=1.0)
        eq(x, -1.07003862e+08)
        eq(y, 9.16969864e+07)
        eq(z, 3.98134914e+07)
        eq(dx, -1.79386097e+06)
        eq(dy, -1.73625791e+06)
        eq(dz, -7.54793163e+05)

    def test_scalar(self):
        import de405
        e = Ephemeris(de405)

        self.check0(*e.compute('earthmoon', 2305447.5))
        self.check1(*e.compute('earthmoon', 2305478.5))

    def test_array(self):
        import de405
        e = Ephemeris(de405)

        v = e.compute('earthmoon', np.array([2305447.5, 2305478.5]))
        v = np.array(v)
        self.check0(*v[:,0])
        self.check1(*v[:,1])

    def test_scalar_at_ephemeris_end(self):
        # TODO
        pass

    def test_array_at_ephemeris_end(self):
        # TODO
        pass
