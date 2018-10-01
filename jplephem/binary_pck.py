"""Compute things from a NASA SPICE binary PCK kernel file.

ftp://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/pck.html

"""
from numpy import array, empty, empty_like, rollaxis
from .daf import DAF
from .names import target_names

T0 = 2451545.0
S_PER_DAY = 86400.0


def jd(seconds):
    """Convert a number of seconds since J2000 to a Julian Date."""
    return T0 + seconds / S_PER_DAY


class BinaryPCK(object):
    """A JPL binary PCK (extension ``.bcp``) kernel.

    You can load a binary PCK file by specifying its filename::

        kernel = BinaryPCK.open('moon_pa_de421_1900-2050.bpc')

    Run ``print(kernel)`` see which segments are inside and iterate
    across ``kernel.segments`` to access them each in turn.

    To see the text comments, call ``kernel.comments()``.

    """
    def __init__(self, daf):
        self.daf = daf
        self.segments = [Segment(self.daf, source, descriptor)
                         for source, descriptor in self.daf.summaries()]

    @classmethod
    def open(cls, path):
        """Open the file at `path` and return a binary PCK instance."""
        return cls(DAF(open(path, 'rb')))

    def close(self):
        """Close this file."""
        self.daf.file.close()
        for segment in self.segments:
            if hasattr(segment, '_data'):
                del segment._data  # TODO: explicitly close each memory map

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
    """A single segment of a binary PCK file.

    There are several items of information about each segment that are
    loaded from the underlying PCK file, and made available as object
    attributes:

    segment.source - official ephemeris name, like 'DE-0430LE-0430'
    segment.start_second - initial epoch, as seconds from J2000
    segment.end_second - final epoch, as seconds from J2000
    segment.start_jd - start_second, converted to a Julian Date
    segment.end_jd - end_second, converted to a Julian Date
    segment.body - integer body identifier
    segment.frame - integer frame identifier
    segment.data_type - integer data type identifier
    segment.start_i - index where segment starts
    segment.end_i - index where segment ends

    """
    def __init__(self, daf, source, descriptor):
        self.daf = daf
        self.source = source
        (self.start_second, self.end_second, self.body, self.frame,
         self.data_type, self.start_i, self.end_i) = descriptor
        self.start_jd = jd(self.start_second)
        self.end_jd = jd(self.end_second)
        self._data = None

    def __str__(self):
        return self.describe(verbose=False)

    def describe(self, verbose=True):
        """Return a textual description of the segment."""
        body = titlecase(target_names.get(self.body, 'Unknown body'))
        text = ('{0.start_jd:.2f}..{0.end_jd:.2f} frame={0.frame}'
                '  {1} ({0.body})'.format(self, body))
        if verbose:
            text += ('\n  data_type={0.data_type} source={1}'
                     .format(self, self.source.decode('ascii')))
        return text

    def _load(self):
        """Map the coefficients into memory using a NumPy array.

        """
        if self.data_type == 2:
            component_count = 3
        else:
            raise ValueError('only binary PCK data type 2 is supported')

        init, intlen, rsize, n = self.daf.read_array(self.end_i - 3, self.end_i)
        initial_epoch = jd(init)
        interval_length = intlen / S_PER_DAY
        coefficient_count = int(rsize - 2) // component_count
        coefficients = self.daf.map_array(self.start_i, self.end_i - 4)

        coefficients.shape = (int(n), int(rsize))
        coefficients = coefficients[:,2:]  # ignore MID and RADIUS elements
        coefficients.shape = (int(n), component_count, coefficient_count)
        coefficients = rollaxis(coefficients, 1)
        return initial_epoch, interval_length, coefficients

    def compute(self, tdb, tdb2, derivative=True):
        """Generate angles and derivatives for time `tdb` plus `tdb2`.

        If ``derivative`` is true, return a tuple containing both the
        angle and its derivative; otherwise simply return the angles.

        """
        scalar = not getattr(tdb, 'shape', 0) and not getattr(tdb2, 'shape', 0)
        if scalar:
            tdb = array((tdb,))

        data = self._data
        if data is None:
            self._data = data = self._load()

        initial_epoch, interval_length, coefficients = data
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

        if not derivative:
            return components

        # Chebyshev differentiation.

        dT = empty_like(T)
        dT[0] = 0.0
        dT[1] = 1.0
        if coefficient_count > 2:
            dT[2] = twot1 + twot1
            for i in range(3, coefficient_count):
                dT[i] = twot1 * dT[i-1] - dT[i-2] + T[i-1] + T[i-1]
        dT *= 2.0
        dT /= interval_length

        rates = (dT.T * coefficients).sum(axis=2)
        if scalar:
            rates = rates[:,0]

        return components, rates


def titlecase(name):
    """Title-case body `name` if it looks safe to do so."""
    return name if name.startswith(('1', 'C/', 'DSS-')) else name.title()
