#!/usr/bin/env python

import numpy as np
import de421
from time import time
from jplephem import Ephemeris
from jplephem.spk import SPK

def main():
    for size in 10, 1000, 100000:
        jd = np.linspace(2414992.5, 2471184.50, size)
        spk = SPK.open('de421.bsp')
        ephem = Ephemeris(de421)
        mars = spk.targets[4]

        print(size)
        print('-- old code (2 successive runs):')

        t0 = time()
        ephem.position('mars', jd)
        print(time() - t0)

        t0 = time()
        ephem.position('mars', jd)
        print(time() - t0)

        print('-- new SPK-powered code (2 successive runs):')

        t0 = time()
        mars.compute(jd)
        print(time() - t0)

        t0 = time()
        mars.compute(jd)
        print(time() - t0)

        print()

if __name__ == '__main__':
    main()
    print(' Warmed up, running again '.center(72, '-'))
    main()
