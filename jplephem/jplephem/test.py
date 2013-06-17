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

    def setUp(self):
        import de421
        self.e = Ephemeris(de421)

    def check0(self, xyz, xyzdot=None):
        eq = partial(self.assertAlmostEqual, delta=1.0)
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
        eq = partial(self.assertAlmostEqual, delta=1.0)
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
        self.assertEqual(self.e.names,  (
            'earthmoon', 'jupiter', 'librations', 'mars', 'mercury',
            'moon', 'neptune', 'nutations', 'pluto', 'saturn', 'sun',
            'uranus', 'venus',
            ))

    def test_scalar_tdb(self):
        self.check0(self.e.position('earthmoon', 2414994.0))
        self.check1(self.e.position('earthmoon', 2415112.5))

    def test_scalar_tdb2(self):
        self.check0(self.e.position('earthmoon', 2414990.0, 4.0))
        self.check1(self.e.position('earthmoon', 2415110.0, 2.5))

    def test_scalar_tdb_keyword(self):
        self.check0(self.e.position('earthmoon', tdb=2414994.0))
        self.check1(self.e.position('earthmoon', tdb=2415112.5))

    def test_scalar_tdb2_keyword(self):
        self.check0(self.e.position('earthmoon', tdb=2414990.0, tdb2=4.0))
        self.check1(self.e.position('earthmoon', tdb=2415110.0, tdb2=2.5))

    def check_2d_result(self, name, tdb, tdb2):
        p = self.e.position(name, tdb, tdb2)
        self.check0(p[:,0])
        self.check1(p[:,1])

        p, v = self.e.position_and_velocity(name, tdb, tdb2)
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
        pv = self.e.compute('earthmoon', 2414994.0)
        self.check0(pv[:3], pv[3:])
        pv = self.e.compute('earthmoon', np.array([2414994.0, 2415112.5]))
        self.check0(pv[:3,0], pv[3:,0])
        self.check1(pv[:3,1], pv[3:,1])

    def test_ephemeris_end_date(self):
        x, y, z = self.e.position('earthmoon', self.e.jomega)
        self.assertAlmostEqual(x, -2.81196460e+07, delta=1.0)
        self.assertAlmostEqual(y, 1.32000379e+08, delta=1.0)
        self.assertAlmostEqual(z, 5.72139011e+07, delta=1.0)

    def test_too_early_date(self):
        tdb = self.e.jalpha - 0.01
        self.assertRaises(DateError, self.e.position, 'earthmoon', tdb)

    def test_too_late_date(self):
        tdb = self.e.jomega + 16.01
        self.assertRaises(DateError, self.e.position, 'earthmoon', tdb)
