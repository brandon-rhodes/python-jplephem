"""Tests for ``jplephem``.

See the accompanying ``jpltest`` module for a more intense numerical
test suite that can verify that ``jplephem`` delivers, over hundreds of
examples, the same results as when the ephemerides are run at JPL.  This
smaller and more feature-oriented suite can be run with::

    python -m unittest discover jplephem

"""
import numpy as np
from functools import partial
from jplephem import Ephemeris, DateError
try:
    from unittest import SkipTest, TestCase
except ImportError:
    from unittest2 import SkipTest, TestCase

epsilon_m = 0.01


class LegacyTests(TestCase):

    def setUp(self):
        try:
            import de421
        except ImportError:
            raise SkipTest('the "de421" ephemeris package is not installed')
        self.eph = Ephemeris(de421)

    def position(self, name, tdb, tdb2=0.0):
        return self.eph.position(name, tdb, tdb2)

    def position_and_velocity(self, name, tdb, tdb2=0.0):
        return self.eph.position_and_velocity(name, tdb, tdb2)

    def check0(self, xyz, xyzdot=None):
        eq = partial(self.assertAlmostEqual, delta=epsilon_m)
        x, y, z = xyz
        eq(x, 39705023.28)
        eq(y, 131195345.65)
        eq(z, 56898495.41)
        if xyzdot is None:
            return
        dx, dy, dz = xyzdot
        eq(dx, -2524248.19)
        eq(dy, 619970.11)
        eq(dz, 268928.26)

    def check1(self, xyz, xyzdot=None):
        eq = partial(self.assertAlmostEqual, delta=epsilon_m)
        x, y, z = xyz
        eq(x, -144692624.00)
        eq(y, -32707965.14)
        eq(z, -14207167.26)
        if xyzdot is None:
            return
        dx, dy, dz = xyzdot
        eq(dx, 587334.38)
        eq(dy, -2297419.36)
        eq(dz, -996628.74)

    def test_names(self):
        self.assertEqual(self.eph.names,  (
            'earthmoon', 'jupiter', 'librations', 'mars', 'mercury',
            'moon', 'neptune', 'nutations', 'pluto', 'saturn', 'sun',
            'uranus', 'venus',
            ))

    def test_scalar_tdb(self):
        self.check0(self.position('earthmoon', 2414994.0))
        self.check1(self.position('earthmoon', 2415112.5))

    def test_scalar_tdb2(self):
        self.check0(self.position('earthmoon', 2414990.0, 4.0))
        self.check1(self.position('earthmoon', 2415110.0, 2.5))

    def test_scalar_tdb_keyword(self):
        self.check0(self.position('earthmoon', tdb=2414994.0))
        self.check1(self.position('earthmoon', tdb=2415112.5))

    def test_scalar_tdb2_keyword(self):
        self.check0(self.position('earthmoon', tdb=2414990.0, tdb2=4.0))
        self.check1(self.position('earthmoon', tdb=2415110.0, tdb2=2.5))

    def check_2d_result(self, name, tdb, tdb2):
        p = self.position(name, tdb + tdb2)
        self.check0(p[:,0])
        self.check1(p[:,1])

        p = self.position(name, tdb, tdb2)
        self.check0(p[:,0])
        self.check1(p[:,1])

        p, v = self.position_and_velocity(name, tdb + tdb2)
        self.check0(p[:,0], v[:,0])
        self.check1(p[:,1], v[:,1])

        p, v = self.position_and_velocity(name, tdb, tdb2)
        self.check0(p[:,0], v[:,0])
        self.check1(p[:,1], v[:,1])

    def test_array_tdb(self):
        tdb = np.array([2414994.0, 2415112.5])
        tdb2 = 0.0
        self.check_2d_result('earthmoon', tdb, tdb2)

    def test_array_tdb_scalar_tdb2(self):
        tdb = np.array([2414991.5, 2415110.0])
        tdb2 = 2.5
        self.check_2d_result('earthmoon', tdb, tdb2)

    def test_scalar_tdb_array_tdb2(self):
        tdb = 2414990.0
        d = 2415112.5 - tdb
        tdb2 = np.array([4.0, d])
        self.check_2d_result('earthmoon', tdb, tdb2)

    def test_array_tdb_array_tdb2(self):
        tdb = np.array([2414990.0, 2415110.0])
        tdb2 = np.array([4.0, 2.5])
        self.check_2d_result('earthmoon', tdb, tdb2)

    def test_legacy_compute_method(self):
        pv = self.eph.compute('earthmoon', 2414994.0)
        self.check0(pv[:3], pv[3:])
        pv = self.eph.compute('earthmoon', np.array([2414994.0, 2415112.5]))
        self.check0(pv[:3,0], pv[3:,0])
        self.check1(pv[:3,1], pv[3:,1])

    def test_ephemeris_end_date(self):
        x, y, z = self.position('earthmoon', self.eph.jomega)
        self.assertAlmostEqual(x, -94189805.73967789, delta=epsilon_m)
        self.assertAlmostEqual(y, 1.05103857e+08, delta=1.0)
        self.assertAlmostEqual(z, 45550861.44383482, delta=epsilon_m)

    def test_too_early_date(self):
        tdb = self.eph.jalpha - 0.01
        self.assertRaises(DateError, self.position, 'earthmoon', tdb)

    def test_too_late_date(self):
        tdb = self.eph.jomega + 16.01
        self.assertRaises(DateError, self.position, 'earthmoon', tdb)
