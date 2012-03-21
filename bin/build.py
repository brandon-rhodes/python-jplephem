#!/usr/bin/env python

# Read big NASA files full of ASCII ephemeris coefficients, and write
# small numpy files that make the data easy for Python to use.

import numpy as np
import os

def e(s):
    """Convert a string in 0.1D+01 FORTRAN notation into 0.1e+10."""
    return s.replace('D', 'e')

def main():
    topdir = 'ssd.jpl.nasa.gov/pub/eph/planets/ascii'
    dirs = os.listdir(topdir)
    for dirname in dirs:
        if not dirname.startswith('de'):
            continue
        dirpath = os.path.join(topdir, dirname)
        filenames = os.listdir(dirpath)
        headername = [ n for n in filenames if n.startswith('header') ][0]
        with open(os.path.join(dirpath, headername)) as f:
            lines = iter(f)

            while next(lines).strip() != 'GROUP   1030':
                continue
            assert next(lines).strip() == ''
            jeda, jedz, _ = (float(s) for s in e(next(lines)).split())

            while next(lines).strip() != 'GROUP   1040':
                continue
            assert next(lines).strip() == ''
            nconstants = int(next(lines))
            names = []
            while len(names) < nconstants:
                names.extend(next(lines).split())

            while next(lines).strip() != 'GROUP   1041':
                continue
            assert next(lines).strip() == ''
            assert int(next(lines)) == nconstants
            values = []
            while len(values) < nconstants:
                values.extend(float(s) for s in e(next(lines)).split())

            constants = np.zeros(nconstants + 2, dtype=[
                    ('name','a6'), ('value','f8')])
            constants['name'] = names + ['JEDA', 'JEDZ']
            constants['value'] = values + [jeda, jedz]

            while next(lines).strip() != 'GROUP   1050':
                continue
            assert next(lines).strip() == ''
            planet_offsets = np.array(next(lines).split(), np.int_)
            num_coefficients = np.array(next(lines).split(), np.int_)
            coefficient_sets = np.array(next(lines).split(), np.int_)

        for filename in filenames:
            if '1600' not in filename:
                continue
            print filename

            n = 0
            tdatasets = []
            with open(os.path.join(dirpath, filename)) as f:
                lines = iter(f)
                for line in lines:

                    # Read a line reading "3 1018" or whatever,
                    # describing a big block of data.

                    n += 1
                    nn, datalen = (int(s) for s in line.split())
                    assert n == nn

                    # Read the first data line, which starts with two
                    # dates before actually starting the coefficients.

                    j0, j1, datum = (float(s) for s in e(next(lines)).split())

                    # Read the rest of the coefficients in this block.

                    datalen -= 2
                    data = [ datum ]
                    while len(data) < datalen:
                        data.extend(float(s) for s in e(next(lines)).split())

                    # Store the date range and coefficients.

                    tdatasets.append((j0, j1, data))

        # Verify that all time periods are equal, and that the datasets
        # in sequence cover adjacent time periods.

        tdatasets.sort()
        step = tdatasets[0][1] - tdatasets[0][0]
        j = tdatasets[0][0]
        for tds in tdatasets:
            assert j == tds[0]
            j += step
            assert j == tds[1]

        # Remove the (j0, j1) tuples and create a pure list of datasets.

        datasets = [ td[2] for td in tdatasets ]

        # Save each planet's coefficients as its own array.

        planet_offsets -= 2  # We handle the two Julian dates separately
        planet_offsets -= 1  # Python arrays index from zero

        for planet, offset in enumerate(planet_offsets):
            nc = num_coefficients[planet]
            cs = coefficient_sets[planet]
            a = np.array([[
                        dataset[
                            offset + 3 * csi * nc :
                            offset + 3 * csi * nc + nc
                            ]
                        for j in range(3)
                        ] for dataset in datasets for csi in range(cs) ])
            print a.shape
            np.save('series%d' % (planet + 1), a)

        np.save('constants', constants)

if __name__ == '__main__':
    main()
