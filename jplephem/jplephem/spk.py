"""Compute positions from a NASA SPICE SPK ephemeris kernel file.

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

"""
from collections import namedtuple
from numpy import ndarray
from .ephem import Ephemeris
from .daf import DAF

Summary = namedtuple('Summary', 'source start_second stop_second target'
                     ' center frame data_type start_index stop_index')

S_PER_DAY = 86400.0
T0 = 2451545.0

class SPK(Ephemeris):
    """A JPL SPK-format ephemeris that computes positions and velocities."""

    def __init__(self, path):
        self.daf = DAF(path)
        self.sets = {}
        g = self.daf.summaries()
        g = (Summary(source, *values) for source, values in g)
        self.summaries = {s.target: s for s in g}

    def array(self, start, stop):
        """Return the array of floats at words `start` to `stop`, inclusive."""
        data = self.daf.bytes(start, stop)
        return ndarray(stop - start + 1, self.daf.endian + 'd', data)

    def load(self, target):
        s = self.sets.get(target)
        if s is None:
            summary = self.summaries[target]
            if summary.data_type == 2:
                component_count = 3
            elif summary.data_type == 3:
                component_count = 6
            else:
                raise ValueError('only SPK data types 2 and 3 are supported')
            stop = summary.stop_index
            init, intlen, rsize, n = self.array(stop - 3, stop)
            coefficient_count = (rsize - 2) // component_count
            # TODO: use intlen directly to create days_per_set
            self.jalpha = T0 + init / S_PER_DAY
            self.jomega = self.jalpha + intlen * n / S_PER_DAY
            print('omega:', init + intlen * n)
            print(init, intlen, rsize, n)
            s = self.array(summary.start_index, stop-4)
            s.shape = (n, rsize)
            s = s[:,2:]
            s.shape = (n, component_count, coefficient_count)
        return s


def main2():
    T0 = 2451545.0
    s = SPK('jup310.bsp')
    p = s.position(3, T0)
    print(p)
    p = s.position(502, T0)
    print(p)


main2()
