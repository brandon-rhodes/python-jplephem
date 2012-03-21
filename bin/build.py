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

            constants = np.zeros(nconstants, dtype=[
                    ('name','a6'), ('value','f8')])
            constants['name'] = names
            constants['value'] = values

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
        #print starts, coeffs, cosets
        #print np.array(starts) - 3

        np.save('constants', constants)

if __name__ == '__main__':
    main()
