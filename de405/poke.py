# First experiment in poking at the DE405 data.

import os
from bisect import bisect
from math import floor

AU = 149597870.7

body_names = (None, 'Mercury', 'Venus', 'Earth', 'Mars', 'Jupiter',
              'Saturn', 'Uranus', 'Neptune', 'Pluto', 'Moon', 'Sun',
              'Solar System Barycenter', 'Earth-Moon Barycenter',
              'Nutations', 'Librations')
coordinate_names = ('x', 'y', 'z', 'xdot', 'ydot', 'zdot')

#

class Ephemeris(object):

    def __init__(self, dirpath):
        self.dirpath = dirpath
        filenames = os.listdir(dirpath)

        headername = [ n for n in filenames if n.startswith('header') ][0]
        with open(os.path.join(dirpath, headername)) as f:
            lines = iter(f)
            while next(lines).strip() != 'GROUP   1050':
                continue
            assert next(lines).strip() == ''
            self.starts = [ int(field) for field in next(lines).split() ]
            self.coeffs = [ int(field) for field in next(lines).split() ]
            self.cosets = [ int(field) for field in next(lines).split() ]
            print self.starts

        step = self.starts[-1] + self.coeffs[-1] * self.cosets[-1] * 3

        datanames = [ n for n in filenames if n.startswith('asc') ]
        ranges = []  # each item is ((start, end), [...])
        for dataname in datanames:
            if '2000' not in dataname:
                continue
            with open(os.path.join(dirpath, dataname)) as f:
                body = f.read()
            array = [ float(f) for f in body.replace('D', 'E').split() ]
            n = 0
            i = 0
            while i < len(array):
                n += 1
                assert array[i] == n
                entries = int(array[i + 1])
                jedrange = (array[i + 2], array[i + 3])
                ranges.append((jedrange, array[i + 4 : i + 4 + entries]))
                i += 4 + entries

        ranges.sort()
        self.ranges = ranges

    def compute(self):
        d = 2456007.0
        i = bisect(self.ranges, ((d,),))
        dates, coefficients = self.ranges[i - 1]
        t2 = dates[1] - dates[0]
        t1 = (d - dates[0]) / t2
        print dates, t2, t1

        EARTH = 3 - 1
        print 'starts:', self.starts[EARTH], len(body_names)

        ncf = self.coeffs[EARTH]
        ncm = 3
        na = self.cosets[EARTH]

        dna = float(na)
        dt1 = floor(t1)
        temp = dna * t1
        l = int(temp - dt1) # + 1 offset was because fortran arrays
        print 'ell', l

        tc = 2.0 * (temp % 1.0 + dt1) - 1.0
        print 'tc', tc

        # np = 2
        # nv = 3
        pc = [1.0, tc]
        twot = tc + tc

        for i in range(2, ncf):
            pc.append(twot * pc[-1] - pc[-2])

        answers = []
        for i in range(ncm):
            answers.append(sum(
                    pc[j] * coefficients[self.starts[EARTH]
                                         - 3
                                         # BUF(J,I,L)
                                         + j
                                         + i * ncf
                                         + l * ncf * ncm]
                    for j in range(ncf)
                    ))

        print answers
        print [a/AU for a in answers]

#

def main():
    ephemeris = Ephemeris('ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405')
    ephemeris.compute()
    return
    testpo = open('ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405/testpo.405')
    lines = iter(testpo)
    while next(lines).strip() != 'EOT':
        continue
    for line in lines:
        fields = line.split()
        print fields
        jed = float(fields[2])
        target = int(fields[3])
        center = int(fields[4])
        coordinate_number = int(fields[5])
        coordinate = float(fields[6])
        print jed, body_names[center], '->', body_names[target]
        r = pleph(jed, target, center)
        print r[coordinate_number], coordinate
        break

def pleph(jed, target, center):
    # todo: nutations
    # todo: librations
    bary = True
    LIST = [0.0] * 12
    IPT = [0.0] * 39

    for k in target, center:
        if k <= 10:
            LIST[k] = 2.0
        if k == 10:
            LIST[3] = 2.0
        if k == 3:
            LIST[10] = 2.0
        if k == 13:
            LIST[3] = 2.0

    return (0.0, 0.0, 0.0, 0.0, 0.0, 0.0)

if __name__ == '__main__':
    main()
