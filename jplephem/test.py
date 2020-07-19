"""Tests for ``jplephem``.

See the accompanying ``jpltest`` module for a more intense numerical
test suite that can verify that ``jplephem`` delivers, over hundreds of
examples, the same results as when the ephemerides are run at JPL.  This
smaller and more feature-oriented suite can be run with::

    python -m unittest discover jplephem

"""
import mmap
import numpy as np
import sys
import tempfile
from doctest import DocTestSuite, ELLIPSIS
from functools import partial
from io import BytesIO
from jplephem import Ephemeris, commandline
from jplephem.exceptions import OutOfRangeError
from jplephem.daf import DAF, FTPSTR, NAIF_DAF
from jplephem.spk import SPK
from struct import Struct
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


class TestDAFBytesIO(TestCase):
    def sample_daf(self):
        word = Struct('d').pack
        integer = Struct('i').pack
        return BytesIO(b''.join([
            # Record 1 - File Record
            b'DAF/SPK ',
            b'\x02\x00\x00\x00', # ND
            b'\x03\x00\x00\x00', # NI
            b'Internal Name'.ljust(60, b' '), # LOCIFN
            b'\x03\x00\x00\x00', # FWARD
            b'\x07\x00\x00\x00', # BWARD
            b'\x01\x04\x00\x00', # FREE
            b'LTL-IEEE', # LOCFMT
            b'\0' * 603, # PRENUL
            FTPSTR,
            b'\0' * 297, # PSTNUL

            # Record 2
            b'Comment Record'.ljust(1024, b'\0'),

            # Record 3 - first Summary Record
            b''.join([
                word(7), # next summary record
                word(0), # previous summary record
                word(1), # number of summaries
                word(101),
                word(202),
                integer(303),
                integer(1024 * 4 // 8 + 1), # Record 5 start
                integer(1024 * 5 // 8),     # Record 5 end
                integer(0),
            ]).ljust(1024, b'\0'),

            # Record 4 - first Name Record
            b'Summary Name 1'.ljust(1024, b' '),

            # Record 5
            word(1001) * 128,

            # Record 6
            word(2002) * 128,

            # Record 7 - second Summary Record
            b''.join([
                word(0), # next summary record
                word(3), # previous summary record
                word(1), # number of summaries
                word(111),
                word(222),
                integer(333),
                integer(1024 * 5 // 8 + 1), # Record 6 start
                integer(1024 * 6 // 8),     # Record 6 end
                integer(0),
            ]).ljust(1024, b'\0'),

            # Record 8 - second Name Record
            b'Summary Name 2'.ljust(1024, b' '),
        ]))

    def test_header(self):
        f = self.sample_daf()
        d = DAF(f)
        eq = self.assertEqual
        eq(d.locidw, b'DAF/SPK')
        eq(d.nd, 2)
        eq(d.ni, 3)
        eq(d.locifn_text, b'Internal Name')
        eq(d.fward, 3)
        eq(d.bward, 7)
        eq(d.free, 0x401)
        eq(d.locfmt, b'LTL-IEEE')

    def test_segments(self):
        f = self.sample_daf()
        d = DAF(f)

        summaries = list(d.summaries())
        eq = self.assertEqual
        eq(len(summaries), 2)
        eq(summaries[0], (b'Summary Name 1', (101.0, 202.0, 303, 513, 640)))
        eq(summaries[1], (b'Summary Name 2', (111.0, 222.0, 333, 641, 768)))

        eq = self.assertSequenceEqual
        eq(list(d.map(summaries[0][1])), [1001.0] * 128)
        eq(list(d.map(summaries[1][1])), [2002.0] * 128)

    def test_add_segment(self):
        f = self.sample_daf()
        d = DAF(f)

        d.add_array(b'Summary Name 3', (121.0, 232.0, 343), [3003.0] * 128)

        summaries = list(d.summaries())
        eq = self.assertEqual
        eq(len(summaries), 3)
        eq(summaries[0], (b'Summary Name 1', (101.0, 202.0, 303, 513, 640)))
        eq(summaries[1], (b'Summary Name 2', (111.0, 222.0, 333, 641, 768)))
        eq(summaries[2], (b'Summary Name 3', (121.0, 232.0, 343, 1025, 1152)))

        eq = self.assertSequenceEqual
        eq(list(d.map(summaries[0][1])), [1001.0] * 128)
        eq(list(d.map(summaries[1][1])), [2002.0] * 128)
        eq(list(d.map(summaries[2][1])), [3003.0] * 128)

    def test_add_segment_when_summary_block_is_full(self):
        f = self.sample_daf()
        d = DAF(f)

        # Update n_summaries of final summary block to full.
        d.file.seek(6 * 1024 + 16)
        d.file.write(Struct('d').pack(d.summaries_per_record))

        d.add_array(b'Summary Name 3', (121.0, 232.0, 343), [3003.0] * 200)

        # Reset n_summaries of that block back to its real value.
        d.file.seek(6 * 1024 + 16)
        d.file.write(Struct('d').pack(1))

        summaries = list(d.summaries())
        eq = self.assertEqual
        eq(len(summaries), 3)
        eq(summaries[0], (b'Summary Name 1', (101.0, 202.0, 303, 513, 640)))
        eq(summaries[1], (b'Summary Name 2', (111.0, 222.0, 333, 641, 768)))
        eq(summaries[2], (b'Summary Name 3', (121.0, 232.0, 343, 1281, 1480)))

        eq = self.assertSequenceEqual
        eq(list(d.map(summaries[0][1])), [1001.0] * 128)
        eq(list(d.map(summaries[1][1])), [2002.0] * 128)
        eq(list(d.map(summaries[2][1])), [3003.0] * 200)


class TestDAFRealFile(TestDAFBytesIO):
    # Where "Real" = "written to disk with a real file descriptor
    # instead of an in-memory BytesIO".

    def sample_daf(self):
        bytes_io = super(TestDAFRealFile, self).sample_daf()
        f = tempfile.NamedTemporaryFile(mode='w+b', prefix='jplephem_test')
        f.write(bytes_io.getvalue())
        f.seek(0)
        return f


def fake_mmap_that_raises_OSError(*args, **kw):
    raise OSError('mmap() not supported on this platform')

class TestDAFRealFileWithoutMMap(TestDAFRealFile):
    # Where "Real" = "written to disk with a real file descriptor
    # instead of an in-memory BytesIO".  And we turn off mmap() to
    # simulate platforms like pyodide.

    def setUp(self):
        self.mmap = mmap.mmap
        mmap.mmap = fake_mmap_that_raises_OSError

    def tearDown(self):
        mmap.mmap = self.mmap


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

    def test_jitter(self):
        usecond = 1e-6 / 24.0 / 3600.0
        tdb = np.ones(12) * 2414998.0
        tdb2 = np.linspace(1.0 * usecond, 2.0 * usecond, 12)
        x, y, z = self.position('earthmoon', tdb, tdb2)
        for component in x, y, z:
            size = component[0]
            relative_jitter = np.diff(np.diff(x)) / size
            self.assertLess(max(abs(relative_jitter)), 3e-16)

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
  '2414864.50..2471184.50  Type 2  Solar System Barycenter (0) -> Mars Barycenter (4)')
        self.assertEqual(segment.describe(verbose=True),
  '2414864.50..2471184.50  Type 2  Solar System Barycenter (0) -> Mars Barycenter (4)'
  '\n  frame=1 source=DE-0421LE-0421')

    def test_loading_array(self):
        segment = self.spk[0,4]
        initial_epoch, interval_length, coefficients = segment.load_array()
        self.assertEqual(coefficients.shape, (3, 1760, 11))

    def test_out_of_range_dates(self):
        segment = self.spk[0,4]
        tdb = np.array([-1e3, 0, +1e5]) + 2414990.0
        try:
            segment.compute_and_differentiate(tdb)
        except OutOfRangeError as e:
            self.assertEqual(str(e), 'segment only covers dates'
                             ' 2414864.5 through 2471184.5')
            self.assertIs(type(e.out_of_range_times), np.ndarray)
            self.assertEqual(list(e.out_of_range_times), [True, False, True])

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
        with SPK(NAIF_DAF(open('de405.bsp', 'rb'))) as kernel:
            x, y, z = kernel[0,4].compute(2457061.5)
            # Expect rough agreement with a DE430 position from our README:
            self.assertAlmostEqual(x, 2.05700211e+08, delta=2.0)
            self.assertAlmostEqual(y, 4.25141646e+07, delta=2.0)
            self.assertAlmostEqual(z, 1.39379183e+07, delta=2.0)


class CommandLineTests(TestCase):
    maxDiff = 9999

    def test_comment_command(self):
        output = commandline.main(['comment', 'de405.bsp'])
        self.assertEqual(output[:30], '; de405.bsp LOG FILE\n;\n; Creat')
        self.assertEqual(output[-30:], "rom Standish's DE405 memo <<<\n")

    def test_daf_command(self):
        self.assertEqual(commandline.main(['daf', 'de405.bsp']), """\
 1 DE-405 -1577879958.8160586 1577880064.1839132 1 0 1 2 1409 202316
 2 DE-405 -1577879958.8160586 1577880064.1839132 2 0 1 2 202317 275376
 3 DE-405 -1577879958.8160586 1577880064.1839132 3 0 1 2 275377 368983
 4 DE-405 -1577879958.8160586 1577880064.1839132 4 0 1 2 368984 408957
 5 DE-405 -1577879958.8160586 1577880064.1839132 5 0 1 2 408958 438653
 6 DE-405 -1577879958.8160586 1577880064.1839132 6 0 1 2 438654 464923
 7 DE-405 -1577879958.8160586 1577880064.1839132 7 0 1 2 464924 487767
 8 DE-405 -1577879958.8160586 1577880064.1839132 8 0 1 2 487768 510611
 9 DE-405 -1577879958.8160586 1577880064.1839132 9 0 1 2 510612 533455
10 DE-405 -1577879958.8160586 1577880064.1839132 10 0 1 2 533456 613364
11 DE-405 -1577879958.8160586 1577880064.1839132 301 3 1 2 613365 987780
12 DE-405 -1577879958.8160586 1577880064.1839132 399 3 1 2 987781 1362196
13 DE-405 -1577879958.8160586 1577880064.1839132 199 1 1 2 1362197 1362208
14 DE-405 -1577879958.8160586 1577880064.1839132 299 2 1 2 1362209 1362220
15 DE-405 -1577879958.8160586 1577880064.1839132 499 4 1 2 1362221 1362232
""")

    def test_spk_command(self):
        self.assertEqual(commandline.main(['spk', 'de405.bsp']), """\
File type NAIF/DAF and format BIG-IEEE with 15 segments:
2433282.50..2469807.50  Type 2  Solar System Barycenter (0) -> Mercury Barycenter (1)
2433282.50..2469807.50  Type 2  Solar System Barycenter (0) -> Venus Barycenter (2)
2433282.50..2469807.50  Type 2  Solar System Barycenter (0) -> Earth Barycenter (3)
2433282.50..2469807.50  Type 2  Solar System Barycenter (0) -> Mars Barycenter (4)
2433282.50..2469807.50  Type 2  Solar System Barycenter (0) -> Jupiter Barycenter (5)
2433282.50..2469807.50  Type 2  Solar System Barycenter (0) -> Saturn Barycenter (6)
2433282.50..2469807.50  Type 2  Solar System Barycenter (0) -> Uranus Barycenter (7)
2433282.50..2469807.50  Type 2  Solar System Barycenter (0) -> Neptune Barycenter (8)
2433282.50..2469807.50  Type 2  Solar System Barycenter (0) -> Pluto Barycenter (9)
2433282.50..2469807.50  Type 2  Solar System Barycenter (0) -> Sun (10)
2433282.50..2469807.50  Type 2  Earth Barycenter (3) -> Moon (301)
2433282.50..2469807.50  Type 2  Earth Barycenter (3) -> Earth (399)
2433282.50..2469807.50  Type 2  Mercury Barycenter (1) -> Mercury (199)
2433282.50..2469807.50  Type 2  Venus Barycenter (2) -> Venus (299)
2433282.50..2469807.50  Type 2  Mars Barycenter (4) -> Mars (499)
""")


def load_tests(loader, tests, ignore):
    """Run our main documentation as a test."""

    # If we are running in CI, where we test against an old version of
    # NumPy, skip the doctests since NumPy will print whitespace
    # differently (and worse).
    version = tuple(int(s) for s in np.__version__.split('.'))
    if version < (1, 17):
        return tests

    # Python 2.6 formats floating-point numbers a bit differently and
    # breaks the doctest.
    if sys.version_info <= (2, 6):
        return tests

    tests.addTests(DocTestSuite('jplephem', optionflags=ELLIPSIS))
    return tests
