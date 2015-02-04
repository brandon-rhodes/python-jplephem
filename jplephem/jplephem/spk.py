"""Compute positions from a NASA SPICE SPK ephemeris kernel file.

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

"""
from collections import namedtuple
from numpy import array, empty, empty_like, ndarray, rollaxis
from .daf import DAF

Segment = namedtuple('Segment', 'source start_second stop_second target'
                     ' center frame data_type start_index stop_index')

S_PER_DAY = 86400.0
T0 = 2451545.0


class SPK(object):
    """A JPL SPK ephemeris kernel for computing positions and velocities."""

    def __init__(self, path):
        self.daf = DAF(path)
        g = self.daf.summaries()
        self.segment_list = [Segment(source, *values) for source, values in g]
        self.segments = {s.target: s for s in self.segment_list}
        self._coefficients = {}

    def array(self, start, stop):
        """Return the array of floats from `start` to `stop` inclusive."""
        data = self.daf.bytes(start, stop)
        return ndarray(stop - start + 1, self.daf.endian + 'd', data)

    def _load(self, segment):
        if segment.data_type == 2:
            component_count = 3
        elif segment.data_type == 3:
            component_count = 6
        else:
            raise ValueError('only SPK data types 2 and 3 are supported')
        stop = segment.stop_index
        init, intlen, rsize, n = self.array(stop - 3, stop)
        initial_epoch = T0 + init / S_PER_DAY
        interval_length = intlen / S_PER_DAY
        coefficient_count = (rsize - 2) // component_count
        coefficients = self.array(segment.start_index, stop-4)
        coefficients.shape = (n, rsize)
        coefficients = coefficients[:,2:]  # ignore MID and RADIUS elements
        coefficients.shape = (n, component_count, coefficient_count)
        coefficients = rollaxis(coefficients, 1)
        return initial_epoch, interval_length, coefficients

    def unload(self, segment):
        self._coefficients.pop(segment, None)

    def compute(self, segment, tdb, tdb2=0.0, differentiate=False):
        """Compute the component values for the time `tdb` plus `tdb2`.

        If `differentiate` is false, then an array of components is
        returned.

        If `differentiate` is true, then a tuple is returned whose first
        element is the array of components and whose second element is
        an array of rates at which the components are changing.

        """
        is_scalar = getattr(tdb, 'shape', ()) == ()
        if is_scalar:
            tdb = array((tdb,))

        info = self._coefficients.get(segment)
        if info is None:
            self._coefficients[segment] = info = self._load(segment)

        initial_epoch, interval_length, coefficients = info
        component_count, n, coefficient_count = coefficients.shape

        # Subtracting tdb before adding tdb2 affords greater precision.
        index, offset = divmod((tdb - initial_epoch) + tdb2, interval_length)
        index = index.astype(int)

        if (index < 0).any() or (index > n).any():
            final_epoch = initial_epoch + interval_length * n
            raise ValueError('segment only covers dates %.1f through %.1f'
                            % (initial_epoch, final_epoch))

        omegas = (index == n)
        index[omegas] -= 1
        offset[omegas] += interval_length

        coefficients = coefficients[:,index]

        # Chebyshev polynomial.

        T = empty((coefficient_count, len(index)))
        T[0] = 1.0
        T[1] = t1 = 2.0 * offset / interval_length - 1.0
        twot1 = t1 + t1
        for i in range(2, coefficient_count):
            T[i] = twot1 * T[i-1] - T[i-2]

        components = (T.T * coefficients).sum(axis=2)
        if not differentiate:
            return components

        # Chebyshev differentiation.

        dT = empty_like(T)
        dT[0] = 0.0
        dT[1] = 1.0
        dT[2] = twot1 + twot1
        for i in range(3, coefficient_count):
            dT[i] = twot1 * dT[i-1] - dT[i-2] + T[i-1] + T[i-1]
        dT *= 2.0
        dT /= interval_length

        rates = (dT.T * coefficients).sum(axis=2)
        return components, rates


def main2():
    T0 = 2451545.0
    spk = SPK('jup310.bsp')
    segment = spk.segments[3]
    p = spk.compute(segment, T0)
    print(p)
    segment = spk.segments[502]
    p = spk.compute(segment, T0)
    print(p)


main2()
