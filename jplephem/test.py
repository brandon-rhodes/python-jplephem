"""Test harness for checking jplephem against actual JPL computations."""

import numpy as np
from sys import exit
from .ephem import Ephemeris

def testpo(module, testpo_path):
    """Compare the positions we calculate against those computed by the JPL."""
    ephemeris = Ephemeris(module)
    lines = iter(open(testpo_path))

    while next(lines).strip() != 'EOT':
        continue

    successes = 0

    for line in lines:
        de, date, jed, target, center, number, value = [f(v) for f, v
            in zip((str, str, float, int, int, int, float), line.split())]

        if 14 <= target <= 15:
            r = ephemeris.compute(target - 2, jed)
        else:
            tpos = _position(ephemeris, jed, target)
            cpos = _position(ephemeris, jed, center)
            r = (tpos - cpos) / ephemeris.AU

        delta = r[number - 1] - value
        if abs(delta) >= 1e-13:
            print '%s %s %s->%s field %d' % (date, jed, center, target, number)
            print 'JPL result: %.15f' % value
            print 'Our result: %.15f' % r[number - 1]
            print 'ERROR: difference =', delta
            exit(1)

        successes += 1
    print '%d tests successful' % successes


def _position(ephemeris, jed, target):
    """Compute position given a JPL test file target integer identifier."""
    if target == 12:
        return np.zeros(6)  # solar system barycenter is the origin
    c = ephemeris.compute
    if target == 13:
        return c(3, jed)
    if target == 3:
        return c(3, jed) - c(10, jed) * ephemeris.earth_share
    if target == 10:
        return c(3, jed) + c(10, jed) * ephemeris.moon_share
    return c(target, jed)


if __name__ == '__main__':
    import de405
    testpo(de405, 'ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/testpo.405')
