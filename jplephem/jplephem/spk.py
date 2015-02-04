"""Compute positions from a NASA SPICE SPK ephemeris kernel file.

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

"""
import numpy as np
from .ephem import Ephemeris
from .daf import DAF

S_PER_DAY = 86400.0
T0 = 2451545.0

class SPK(Ephemeris):
    """A JPL SPK-format ephemeris that computes positions and velocities."""

    def __init__(self, path):
        self.daf = DAF(path)
        self.sets = {}
        self.summaries = {s.target: s for s in self.daf.summaries()}

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
            init, intlen, rsize, n = self.daf[stop-3:stop]
            coefficient_count = (rsize - 2) // component_count
            print(summary)
            # TODO: use intlen directly to create days_per_set
            self.jalpha = T0 + init / S_PER_DAY
            self.jomega = self.jalpha + intlen * n / S_PER_DAY
            print('omega:', init + intlen * n)
            print(init, intlen, rsize, n)
            data = self.daf.bytes(summary.start_index, stop)
            s = np.ndarray((n, rsize), self.daf.endian + 'd', data)
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
