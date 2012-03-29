#!/usr/bin/env python

# Read big NASA files full of ASCII ephemeris coefficients, and write
# small numpy files that make the data easy for Python to use.

import numpy as np
import os
from operator import itemgetter

def e(s):
    """Convert a string in 0.1D+01 FORTRAN notation into 0.1e+10."""
    return s.replace('D', 'e')

def main():
    topdir = 'ssd.jpl.nasa.gov/pub/eph/planets/ascii'
    dirs = os.listdir(topdir)
    for dirname in sorted(dirs):
        if not dirname.startswith('de'):
            continue
        print '=' * 16, dirname
        if not os.path.isdir(dirname):
            os.mkdir(dirname)
        with open(os.path.join(dirname, '__init__.py'), 'w') as f:
            pass
        dirpath = os.path.join(topdir, dirname)
        filenames = os.listdir(dirpath)
        headername = [ n for n in filenames if n.startswith('header') ][0]
        with open(os.path.join(dirpath, headername)) as f:
            lines = iter(f)

            while next(lines).strip() != 'GROUP   1030':
                continue
            assert next(lines).strip() == ''
            jalpha, jomega, jdelta = (float(s) for s in e(next(lines)).split())

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

            constants = np.zeros(nconstants + 3, dtype=[
                    ('name','a6'), ('value','f8')])
            constants['name'] = names + ['jalpha', 'jomega', 'jdelta']
            constants['value'] = values + [jalpha, jomega, jdelta]

            while next(lines).strip() != 'GROUP   1050':
                continue
            assert next(lines).strip() == ''
            planet_offsets = np.array(next(lines).split(), np.int_)
            num_coefficients = np.array(next(lines).split(), np.int_)
            coefficient_sets = np.array(next(lines).split(), np.int_)

        tdatasets = []

        for filename in sorted(filenames):
            if not filename.startswith('asc'):
                continue
            print filename

            n = 0
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
                    dataset = [ datum ]
                    while len(dataset) < datalen:
                        dataset.extend(
                            float(s) for s in e(next(lines)).split()
                            )

                    # Store the date range and coefficients.

                    dataset = np.array(dataset)  # array uses half the memory
                    tdatasets.append((j0, j1, dataset))

        # Sort all coefficient sets and remove adacent duplicates, which
        # occur because the last record of one file is sometimes (!) the
        # first record of the next file.

        tdatasets.sort(key=itemgetter(0))
        for i in reversed(range(0, len(tdatasets) - 1)):
            if (tdatasets[i][0] == tdatasets[i+1][0] and
                tdatasets[i][1] == tdatasets[i+1][1] and
                (tdatasets[i][2] == tdatasets[i+1][2]).all()):
                del tdatasets[i+1]

        # Verify that all time periods are equal, and that the datasets
        # in sequence cover adjacent time periods.

        step = tdatasets[0][1] - tdatasets[0][0]
        j = tdatasets[0][0]
        for tds in tdatasets:
            assert j == tds[0], (j, tds[0])
            j += step
            assert j == tds[1], (j, tds[1])

        # Remove the (j0, j1) tuples and create a pure list of datasets.

        datasets = [ td[2] for td in tdatasets ]

        # Save each planet's coefficients as its own array.

        planet_offsets -= 2  # We handle the two Julian dates separately
        planet_offsets -= 1  # Python arrays index from zero

        for planet, offset in enumerate(planet_offsets):
            nc = num_coefficients[planet]
            cs = coefficient_sets[planet]
            ncm = 2 if planet == 11 else 3
            a = np.array([[
                        dataset[
                            offset + j * nc + ncm * csi * nc :
                            offset + j * nc + ncm * csi * nc + nc
                            ]
                        for j in range(ncm)
                        ] for dataset in datasets for csi in range(cs) ])

            print 'polynomials', planet + 1, a.shape
            np.save(os.path.join(dirname, 'jpl-%02d' % (planet + 1)), a)

        np.save(os.path.join(dirname, 'constants'), constants)

if __name__ == '__main__':
    main()
