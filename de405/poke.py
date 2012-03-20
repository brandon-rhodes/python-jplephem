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
        for filename in os.listdir(dirpath):
            pass

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
