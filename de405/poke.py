# First experiment in poking at the DE405 data.

import os

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
            with open(os.path.join(dirpath, dataname)) as f:
                body = f.read()
            array = [ float(f) for f in body.replace('D', 'E').split() ]
            n = 0
            while n < len(array):
                count = int(array[n + 1])
                jedrange = (array[n + 2], array[n + 3])
                ranges.append((jedrange, array[n + 4:n + 4 + count]))
                n += 4 + count

        ranges.sort()
        print ranges[0][0]
        print ranges[-2][0]
        print ranges[-1][0]
#

def main():
    ephemeris = Ephemeris('ssd.jpl.nasa.gov/pub/eph/planets/ascii/de405')
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
