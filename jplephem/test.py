"""Tests for ``jplephem``.

See the accompanying ``jpltest`` module for a more intense numerical
test suite that can verify that ``jplephem`` delivers, over hundreds of
examples, the same results as when the ephemerides are run at JPL.  This
smaller and more feature-oriented suite can be run with::

    python -m unittest discover jplephem

"""
import numpy as np
import subprocess
from functools import partial
from jplephem import Ephemeris, commandline
from jplephem.daf import NAIF_DAF
from jplephem.spk import SPK
try:
    from unittest import SkipTest, TestCase
except ImportError:
    from unittest2 import SkipTest, TestCase

epsilon_m = 0.01
target_names = {
    'mercury barycenter': 1, # BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
    'venus barycenter': 2, # BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
    'earthmoon': 3,        # BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
    'mars barycenter': 4,  # BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
    'jupiter': 5,          # BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
    'saturn': 6,           # BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
    'uranus': 7,           # BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
    'neptune': 8,          # BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
    'pluto': 9,            # BARYCENTER w.r.t. 0 SOLAR SYSTEM BARYCENTER
    'sun': 10,             # w.r.t. 0 SOLAR SYSTEM BARYCENTER
    'mercury': 199,        # w.r.t. 1 MERCURY BARYCENTER
    'venus': 299,          # w.r.t. 2 VENUS BARYCENTER
    'moon': 301,           # w.r.t. 3 EARTH BARYCENTER
    'earth': 399,          # w.r.t. 3 EARTH BARYCENTER
    'mars': 499,           # w.r.t. 4 MARS BARYCENTER
    }


class _CommonTests(object):

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

    def test_ephemeris_end_date(self):
        x, y, z = self.position('earthmoon', self.jomega)
        # These positions are actually from HORIZONS and thus DE431,
        # hence the low precision match:
        self.assertAlmostEqual(x, 1.442502234663646E+08, delta=1.0)
        self.assertAlmostEqual(y, 3.690043031712407E+07, delta=1.0)
        self.assertAlmostEqual(z, 1.599543968176661E+07, delta=1.0)

    def test_too_early_date(self):
        tdb = self.jalpha - 0.01
        self.assertRaises(ValueError, self.position, 'earthmoon', tdb)

    def test_too_late_date(self):
        tdb = self.jomega + 16.01
        self.assertRaises(ValueError, self.position, 'earthmoon', tdb)


class SPKTests(_CommonTests, TestCase):

    def setUp(self):
        try:
            self.spk = SPK.open('de421.bsp')
        except IOError:
            raise SkipTest('the "de421.bsp" SPK file is not available')
        segment = self.spk[0,1]
        self.jalpha = segment.start_jd
        self.jomega = segment.end_jd

    def tearDown(self):
        self.spk.close()

    def position(self, name, tdb, tdb2=0.0):
        segment = self.spk[0, target_names[name]]
        return segment.compute(tdb, tdb2)

    def position_and_velocity(self, name, tdb, tdb2=0.0):
        segment = self.spk[0, target_names[name]]
        return segment.compute_and_differentiate(tdb, tdb2)

    def test_segment_with_only_two_coefficients(self):
        tdb = 2414990.0
        tup = target_names['mercury barycenter'], target_names['mercury']
        segment = self.spk[tup]
        segment.compute_and_differentiate(tdb)

    def test_str(self):
        str(self.spk)  # just to confirm it does not raise an exception
        segment = self.spk[0,4]
        self.assertEqual(str(segment), segment.describe(verbose=False))
        self.assertEqual(segment.describe(verbose=False),
  '2414864.50..2471184.50  Solar System Barycenter (0) -> Mars Barycenter (4)')
        self.assertEqual(segment.describe(verbose=True),
  '2414864.50..2471184.50  Solar System Barycenter (0) -> Mars Barycenter (4)'
  '\n  frame=1 data_type=2 source=DE-0421LE-0421')


class LegacyTests(_CommonTests, TestCase):

    def setUp(self):
        try:
            import de421
        except ImportError:
            raise SkipTest('the "de421" ephemeris package has not been'
                           ' installed with "pip install de421"')
        self.eph = Ephemeris(de421)
        self.jalpha = self.eph.jalpha
        self.jomega = self.eph.jomega

    def position(self, name, tdb, tdb2=0.0):
        return self.eph.position(name, tdb, tdb2)

    def position_and_velocity(self, name, tdb, tdb2=0.0):
        return self.eph.position_and_velocity(name, tdb, tdb2)

    def test_names(self):
        self.assertEqual(self.eph.names,  (
            'earthmoon', 'jupiter', 'librations', 'mars', 'mercury',
            'moon', 'neptune', 'nutations', 'pluto', 'saturn', 'sun',
            'uranus', 'venus',
            ))

    def test_legacy_compute_method(self):
        pv = self.eph.compute('earthmoon', 2414994.0)
        self.check0(pv[:3], pv[3:])
        pv = self.eph.compute('earthmoon', np.array([2414994.0, 2415112.5]))
        self.check0(pv[:3,0], pv[3:,0])
        self.check1(pv[:3,1], pv[3:,1])

    def test_ephemeris_end_date(self):
        x, y, z = self.position('earthmoon', self.jomega)
        self.assertAlmostEqual(x, -94189805.73967789, delta=epsilon_m)
        self.assertAlmostEqual(y, 1.05103857e+08, delta=1.0)
        self.assertAlmostEqual(z, 45550861.44383482, delta=epsilon_m)


class NAIF_DAF_Tests(TestCase):

    def test_single_position(self):
        kernel = SPK(NAIF_DAF(open('de405.bsp', 'rb')))
        x, y, z = kernel[0,4].compute(2457061.5)
        # Expect rough agreement with a DE430 position from our README:
        self.assertAlmostEqual(x, 2.05700211e+08, delta=2.0)
        self.assertAlmostEqual(y, 4.25141646e+07, delta=2.0)
        self.assertAlmostEqual(z, 1.39379183e+07, delta=2.0)
        kernel.close()


class CommandLineTests(TestCase):
    maxDiff = 9999

    def test_command_line(self):
        self.assertEqual(commandline.main(['de405.bsp']), """\
File type NAIF/DAF and format BIG-IEEE with 15 segments:
2433282.50..2469807.50  Solar System Barycenter (0) -> Mercury Barycenter (1)
2433282.50..2469807.50  Solar System Barycenter (0) -> Venus Barycenter (2)
2433282.50..2469807.50  Solar System Barycenter (0) -> Earth Barycenter (3)
2433282.50..2469807.50  Solar System Barycenter (0) -> Mars Barycenter (4)
2433282.50..2469807.50  Solar System Barycenter (0) -> Jupiter Barycenter (5)
2433282.50..2469807.50  Solar System Barycenter (0) -> Saturn Barycenter (6)
2433282.50..2469807.50  Solar System Barycenter (0) -> Uranus Barycenter (7)
2433282.50..2469807.50  Solar System Barycenter (0) -> Neptune Barycenter (8)
2433282.50..2469807.50  Solar System Barycenter (0) -> Pluto Barycenter (9)
2433282.50..2469807.50  Solar System Barycenter (0) -> Sun (10)
2433282.50..2469807.50  Earth Barycenter (3) -> Moon (301)
2433282.50..2469807.50  Earth Barycenter (3) -> Earth (399)
2433282.50..2469807.50  Mercury Barycenter (1) -> Mercury (199)
2433282.50..2469807.50  Venus Barycenter (2) -> Venus (299)
2433282.50..2469807.50  Mars Barycenter (4) -> Mars (499)""")
