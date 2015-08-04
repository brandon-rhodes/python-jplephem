"""Support function for parsing JPL ephemeris text files.

This is for parsing a NASA ephemeris text header file, like:

ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de421/header.421

You can use this routine like this::

    from jplephem.ascii import parse_header
    d = parse_header(open('header.421'))

    from pprint import pprint
    pprint(d)
    pprint(dict(zip(d['names'], d['values'])))

"""
import numpy as np


def parse_header(lines):
    lines = iter(lines)

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

    while next(lines).strip() != 'GROUP   1050':
        continue
    assert next(lines).strip() == ''
    planet_offsets = np.array(next(lines).split(), np.int_)
    num_coefficients = np.array(next(lines).split(), np.int_)
    coefficient_sets = np.array(next(lines).split(), np.int_)

    del lines
    return(locals())


def e(s):
    """Convert a string in 0.1D+01 FORTRAN notation into 0.1e+10."""
    return s.replace('D', 'e')
