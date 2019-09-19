"""Test harness for checking ``jplephem`` against actual JPL results.

This test can be invoked with a simple::

    python -m jplephem.jpltest

"""
import numpy as np
from .spk import SPK

AU = eval('0.149597870700000000D+09'.replace('D', 'E'))  # km
meter = 0.001 / AU
epsilon = 0.01 * meter
field_names = {1: 'x', 2: 'y', 3: 'z', 4: 'vx', 5: 'vy', 6: 'vz'}

def run_testpo(spk, testpo_file):
    """Compare the positions we compute against those computed by the JPL."""

    lines = iter(testpo_file)

    while next(lines).strip() != 'EOT':
        continue

    skips = 0
    successes = 0
    targets = set([segment.target for segment in spk.segments])

    # Old codes, special-cased in _position():
    targets.add(11)
    targets.add(12)
    targets.add(13)

    for line in lines:
        de, date, jed, target, center, number, value = [f(v) for f, v
            in zip((str, str, float, int, int, int, float), line.split())]

        if (target not in targets) or (center not in targets):
            assert target == 14 or target == 15
            assert center == 0
            skips += 1
            continue

        if 14 <= target <= 15:
            r = _position(spk, jed, target)
        else:
            tpos = _position(spk, jed, target)
            cpos = _position(spk, jed, center)
            r = (tpos - cpos) / AU

        delta = r[number - 1] - value
        # if (target == 15 and number == 3):
        #     delta = delta / (0.23 * (jed - 2451545.0))
        # elif (target == 15 and number == 6):
        #     delta = delta * 0.01 / (1.0 + (jed - 2451545.0) / 365.25)

        if abs(delta) >= epsilon:
            print('%s %s %s->%s field %d (%s)'
                  % (date, jed, center, target, number, field_names[number]))
            print('  JPL result: %.15f' % value)
            print('  Our result: %.15f' % r[number - 1])
            print('    ERROR: difference = %s' % (delta,))
            exit(1)

        successes += 1

    print('  {0} tests successful, {1} skipped'.format(successes, skips))


def _position(spk, jed, target):
    """Compute position given a JPL test file target integer identifier."""

    if target == 3:
        p1, v1 = spk[0,3].compute_and_differentiate(jed)
        p2, v2 = spk[3,399].compute_and_differentiate(jed)
        p = p1 + p2
        v = v1 + v2
    elif target == 10:
        p1, v1 = spk[0,3].compute_and_differentiate(jed)
        p2, v2 = spk[3,301].compute_and_differentiate(jed)
        p = p1 + p2
        v = v1 + v2
    elif target == 11:
        p, v = spk[0,10].compute_and_differentiate(jed)
    elif target == 12:
        return np.zeros((6,))  # solar system barycenter is the origin
    elif target == 13:
        p, v = spk[0,3].compute_and_differentiate(jed)
    else:
        p, v = spk[0,target].compute_and_differentiate(jed)

    return np.concatenate((p, v))


class MissingFile(Exception):
    pass


def test_all():
    number = 430
    spk_path = 'de%d_test_excerpt.bsp' % number
    testpo_path = 'testpo.%d' % number
    try:
        spk = SPK.open(spk_path)
    except IOError:
        raise MissingFile('cannot open: %s' % spk_path)
    try:
        testpo_file = open(testpo_path)
    except IOError:
        raise MissingFile('cannot open: %s' % testpo_path)
    run_testpo(spk, testpo_file)


if __name__ == '__main__':
    try:
        test_all()
    except MissingFile as e:
        print("""
Cannot find the JPL SPK files or test files against which this test
suite validates that the positions it generates are correct. To fetch
them, run these commands in your current working directory:

  wget http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp
  wget ftp://ssd.jpl.nasa.gov/pub/eph/planets/test-data/430/testpo.430

When you are done running the tests, you can remove the files.
""")
        print(str(e))
        exit(1)
