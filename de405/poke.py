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

    def compute(self, planet, d):
        i = bisect(self.ranges, ((d,),))
        dates, coefficients = self.ranges[i - 1]
        t2 = dates[1] - dates[0]
        t1 = (d - dates[0]) / t2

        ncf = self.coeffs[planet]
        ncm = 3
        na = self.cosets[planet]

        dna = float(na)
        dt1 = floor(t1)
        temp = dna * t1
        l = int(temp - dt1) # + 1 offset was because fortran arrays

        tc = 2.0 * (temp % 1.0 + dt1) - 1.0

        # np = 2
        # nv = 3
        pc = [1.0, tc]
        twot = tc + tc

        for i in range(2, ncf):
            pc.append(twot * pc[-1] - pc[-2])

        answers = []
        for i in range(ncm):
            answers.append(sum(
                    pc[j] * coefficients[self.starts[planet]
                                         - 3
                                         # BUF(J,I,L)
                                         + j
                                         + i * ncf
                                         + l * ncf * ncm]
                    for j in range(ncf)
                    ))

        return [a/AU for a in answers]

#

def main():
    ephemeris = Ephemeris('ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405')
    d0 = 2456007.0
    dates = [ d0 + i for i in range(365) ]

    for planet in range(4):
        xyz = [ ephemeris.compute(planet, d) for d in dates ]
        x = [ k[0] for k in xyz ]
        from math import copysign
        y = [ copysign((k[1]*k[1] + k[2]*k[2]) ** 0.5, k[1]) for k in xyz ]
        from matplotlib import pyplot
        pyplot.plot(x, y)

    pyplot.show()
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
