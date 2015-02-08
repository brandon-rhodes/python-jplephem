"""Compute positions from a NASA SPICE SPK ephemeris kernel file.

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

"""
from numpy import array, empty, empty_like, rollaxis
from .daf import DAF
from .names import target_names

T0 = 2451545.0
S_PER_DAY = 86400.0


def jd(seconds):
    """Convert a number of seconds since J2000 to a Julian Date."""
    return T0 + seconds / S_PER_DAY


class SPK(object):
    """A JPL SPK ephemeris kernel for computing positions and velocities."""

    def __init__(self, path):
        self.daf = DAF(path)
        self.segments = [Segment(self.daf, *t) for t in self.daf.summaries()]
        self.targets = dict((s.target, s) for s in self.segments)  # Python 2.6

    def __str__(self):
        daf = self.daf
        d = lambda b: b.decode('latin-1')
        lines = (str(segment) for segment in self.segments)
        return 'File type {0} and format {1} with {2} segments:\n{3}'.format(
            d(daf.locidw), d(daf.locfmt), len(self.segments), '\n'.join(lines))

    def comments(self):
        """Return the file comments, as a string."""
        return self.daf.comments()


class Segment(object):

    def __init__(self, daf, source, descriptor):
        self.daf = daf
        self.source = source
        (self.start_second, self.end_second, self.target, self.center,
         self.frame, self.data_type, self.start_i, self.end_i) = descriptor
        self.start_jd = jd(self.start_second)
        self.end_jd = jd(self.end_second)

    def __str__(self):
        return self.describe(verbose=False)

    def describe(self, verbose=True):
        center = target_names.get(self.center, 'Unknown center')
        target = target_names.get(self.target, 'Unknown target')
        text = ('{0.start_jd:.2f}..{0.end_jd:.2f}  {1} ({0.center})'
                ' -> {2} ({0.target})'.format(self, center, target))
        if verbose:
            text += ('\n  frame={0.frame} data_type={0.data_type} source={1}'
                     .format(self, self.source.decode('ascii')))
        return text

    def _load(self):
        """Map the coefficients into memory using a NumPy array."""

        if self.data_type == 2:
            component_count = 3
        elif self.data_type == 3:
            component_count = 6
        else:
            raise ValueError('only SPK data types 2 and 3 are supported')

        init, intlen, rsize, n = self.daf.array(self.end_i - 3, self.end_i)
        initial_epoch = jd(init)
        interval_length = intlen / S_PER_DAY
        coefficient_count = int(rsize - 2) // component_count
        coefficients = self.daf.array(self.start_i, self.end_i - 4)

        coefficients.shape = (int(n), rsize)
        coefficients = coefficients[:,2:]  # ignore MID and RADIUS elements
        coefficients.shape = (int(n), component_count, coefficient_count)
        coefficients = rollaxis(coefficients, 1)
        return initial_epoch, interval_length, coefficients

    def compute(self, tdb, tdb2=0.0, differentiate=False):
        """Compute the component values for the time `tdb` plus `tdb2`.

        If `differentiate` is false, then an array of components is
        returned.

        If `differentiate` is true, then a tuple is returned whose first
        element is the array of components and whose second element is
        an array of rates at which the components are changing.

        """
        scalar = not getattr(tdb, 'shape', 0) and not getattr(tdb2, 'shape', 0)
        if scalar:
            tdb = array((tdb,))

        try:
            initial_epoch, interval_length, coefficients = self._data
        except AttributeError:
            self._data = self._load()
            initial_epoch, interval_length, coefficients = self._data

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
        if scalar:
            components = components[:,0]
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
        if scalar:
            rates = rates[:,0]
        return components, rates
