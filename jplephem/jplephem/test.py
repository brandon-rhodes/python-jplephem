"""Tests for ``jplephem``.

See the accompanying ``jpltest`` module for a more intense numerical
test suite that can verify that ``jplephem`` delivers, in a large number
of cases, the same results as when the ephemerides are run at JPL.  This
smaller and more feature-oriented suite can be run with::

    python -m unittest discover jplephem

"""
import numpy as np
from functools import partial
from jplephem import Ephemeris, DateError
from unittest import TestCase

class Tests(TestCase):

    def check0(self, x, y, z, dx, dy, dz):
        eq = partial(self.assertAlmostEqual, delta=1.0)
        eq(x, 39705023.28)
        eq(y, 131195345.65)
        eq(z, 56898495.41)
        eq(dx, -2524248.19)
        eq(dy, 619970.11)
        eq(dz, 268928.26)

    def check1(self, x, y, z, dx, dy, dz):
        eq = partial(self.assertAlmostEqual, delta=1.0)
        eq(x, -144692624.00)
        eq(y, -32707965.14)
        eq(z, -14207167.26)
        eq(dx, 587334.38)
        eq(dy, -2297419.36)
        eq(dz, -996628.74)

    def test_scalar_input(self):
        import de421
        e = Ephemeris(de421)

        self.check0(*e.compute('earthmoon', 2414994.0))
        self.check1(*e.compute('earthmoon', 2415112.5))

    def test_array_input(self):
        import de421
        e = Ephemeris(de421)

        v = e.compute('earthmoon', np.array([2414994.0, 2415112.5]))

        v = np.array(v)
        self.check0(*v[:,0])
        self.check1(*v[:,1])

    def test_ephemeris_end_date(self):
        import de421
        e = Ephemeris(de421)
        x, y, z = e.position('earthmoon', e.jomega)
        self.assertAlmostEqual(x, -2.81196460e+07, delta=1.0)
        self.assertAlmostEqual(y, 1.32000379e+08, delta=1.0)
        self.assertAlmostEqual(z, 5.72139011e+07, delta=1.0)

    def test_too_early_date(self):
        import de421
        e = Ephemeris(de421)
        self.assertRaises(DateError, e.compute, 'earthmoon', e.jalpha - 0.01)

    def test_too_late_date(self):
        import de421
        e = Ephemeris(de421)
        self.assertRaises(DateError, e.compute, 'earthmoon', e.jomega + 16.01)
