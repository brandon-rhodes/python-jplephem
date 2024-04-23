"""Compute positions from a NASA SPICE SPK ephemeris kernel file.

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

"""
from numpy import array, interp, rollaxis
from .calendar import compute_calendar_date
from .daf import DAF
from .descriptorlib import reify
from .exceptions import OutOfRangeError
from .names import target_names

T0 = 2451545.0
S_PER_DAY = 86400.0

def jd(seconds):
    """Convert a number of seconds since J2000 to a Julian Date."""
    return T0 + seconds / S_PER_DAY

class SPK(object):
    """A JPL SPK ephemeris kernel for computing positions and velocities.

    You can load an SPK by specifying its filename::

        kernel = SPK.open('de431.bsp')

    Run ``print(kernel)`` to list the ephemeris segments.  You can also
    loop across all of the segments in the list ``kernel.segments`` or,
    as a convenience, you can select a particular segment by providing a
    center and target integer in square brackets.  So ``kernel[3,399]``
    will select the segment that computes the distance between the
    Earth-Moon barycenter (3) and the Earth itself (399).

    To extract the text comments from the SPK use ``kernel.comments()``.

    """
    def __init__(self, daf):
        self.daf = daf
        self.segments = [
            build_segment(self.daf, source, descriptor)
            for source, descriptor in self.daf.summaries()
        ]
        self.pairs = dict(((s.center, s.target), s) for s in self.segments)

    @classmethod
    def open(cls, path):
        """Open the file at `path` and return an SPK instance."""
        return cls(DAF(open(path, 'rb')))

    def close(self):
        """Close this SPK file."""
        self.daf.file.close()
        for segment in self.segments:
            if '_data' in segment.__dict__:
                del segment._data
        self.daf._array = None
        self.daf._map = None

    def __str__(self):
        daf = self.daf
        d = lambda b: b.decode('latin-1')
        lines = [
            'File type {0} and format {1} with {2} segments:'
            .format(d(daf.locidw), d(daf.locfmt), len(self.segments))
        ]
        lines.extend(str(segment) for segment in self.segments)
        return '\n'.join(lines)

    def __getitem__(self, key):
        """Given (center, target) integers, return the last matching segment."""
        return self.pairs[key]

    def comments(self):
        """Return the file comments, as a string."""
        return self.daf.comments()

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()

def build_segment(daf, source, descriptor):
    data_type = descriptor[5]
    cls = _segment_classes.get(data_type, BaseSegment)
    return cls(daf, source, descriptor)

class BaseSegment(object):
    """A single segment of an SPK file.

    There are several items of information about each segment that are
    loaded from the underlying SPK file, and made available as object
    attributes:

    segment.source - official ephemeris name, like 'DE-0430LE-0430'
    segment.start_second - initial epoch, as seconds from J2000
    segment.end_second - final epoch, as seconds from J2000
    segment.start_jd - start_second, converted to a Julian Date
    segment.end_jd - end_second, converted to a Julian Date
    segment.center - integer center identifier
    segment.target - integer target identifier
    segment.frame - integer frame identifier
    segment.data_type - integer data type identifier
    segment.start_i - index where segment starts
    segment.end_i - index where segment ends

    """
    _data = None

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
        """Return a textual description of the segment."""
        start = '%d-%02d-%02d' % compute_calendar_date(self.start_jd + 0.5)
        end = '%d-%02d-%02d' % compute_calendar_date(self.end_jd + 0.5)
        center = titlecase(target_names.get(self.center, 'Unknown center'))
        target = titlecase(target_names.get(self.target, 'Unknown target'))
        text = ('{1}..{2}  Type {0.data_type}'
                '  {3} ({0.center}) -> {4} ({0.target})'
                .format(self, start, end, center, target))
        if verbose:
            text += ('\n  frame={0.frame} source={1}'
                     .format(self, self.source.decode('ascii')))
        return text

    def compute(self, tdb, tdb2=0.0):
        """Compute the component values for the time `tdb` plus `tdb2`."""
        raise ValueError(
            'jplephem has not yet learned how to compute positions'
            ' from an ephemeris segment with data type {0}'
            .format(self.data_type)
        )

    def compute_and_differentiate(self, tdb, tdb2=0.0):
        """Compute components and differentials for time `tdb` plus `tdb2`."""
        raise ValueError(
            'jplephem has not yet learned how to compute positions and'
            ' velocities from an ephemeris segment with data type {0}'
            .format(self.data_type)
        )


class Segment(BaseSegment):
    # Type 2 or type 3 segment.

    def compute(self, tdb, tdb2=0.0):
        """Compute the component values for the time `tdb` plus `tdb2`."""
        for position in self.generate(tdb, tdb2):
            return position

    def compute_and_differentiate(self, tdb, tdb2=0.0):
        """Compute components and differentials for time `tdb` plus `tdb2`."""
        return tuple(self.generate(tdb, tdb2))

    @reify
    def _data(self):
        """Map the coefficients into memory using a NumPy array.

        """
        if self.data_type == 2:
            component_count = 3
        elif self.data_type == 3:
            component_count = 6
        else:
            raise ValueError('this class only supports SPK data types 2 and 3')

        init, intlen, rsize, n = self.daf.read_array(self.end_i - 3, self.end_i)
        coefficient_count = int(rsize - 2) // component_count
        coefficients = self.daf.map_array(self.start_i, self.end_i - 4)

        coefficients.shape = (int(n), int(rsize))
        coefficients = coefficients[:,2:]  # ignore MID and RADIUS elements
        coefficients.shape = (int(n), component_count, coefficient_count)
        coefficients = rollaxis(coefficients, 1)
        coefficients = rollaxis(coefficients, 2)
        coefficients = coefficients[::-1]
        return init, intlen, coefficients

    def load_array(self):
        init, intlen, coefficients = self._data
        initial_epoch = jd(init)
        interval_length = intlen / S_PER_DAY
        coefficients = coefficients[::-1]
        coefficients = rollaxis(coefficients, 2)
        coefficients = rollaxis(coefficients, 2)
        return initial_epoch, interval_length, coefficients

    def generate(self, tdb, tdb2):
        """Generate components and differentials for time `tdb` plus `tdb2`.

        Most uses will simply want to call the `compute()` method or the
        `compute_differentials()` method, for convenience.  But in those
        cases (see Skyfield) where you want to compute a position and
        examine it before deciding whether to proceed with the velocity,
        but without losing all of the work that it took to get to that
        point, this generator lets you get them as two separate steps.

        """
        scalar = not getattr(tdb, 'shape', 0) and not getattr(tdb2, 'shape', 0)
        if scalar:
            tdb = array((tdb,))

        init, intlen, coefficients = self._data
        coefficient_count, component_count, n = coefficients.shape

        # Keeping fractions strictly separate from whole numbers
        # maintains the highest possible precision.

        index1, offset1 = divmod((tdb - T0) * S_PER_DAY - init, intlen)
        index2, offset2 = divmod(tdb2 * S_PER_DAY, intlen)
        index3, offset = divmod(offset1 + offset2, intlen)
        index = (index1 + index2 + index3).astype(int)

        if (index < 0).any() or (index > n).any():
            raise OutOfRangeError(
                'segment only covers dates %d-%02d-%02d through %d-%02d-%02d'
                % (compute_calendar_date(self.start_jd + 0.5) +
                   compute_calendar_date(self.end_jd + 0.5)),
                out_of_range_times=(index < 0) | (index > n),
            )

        omegas = (index == n)
        index[omegas] -= 1
        offset[omegas] += intlen

        coefficients = coefficients[:,:,index]

        # Chebyshev polynomial.

        s = 2.0 * offset / intlen - 1.0
        s2 = 2.0 * s

        w0 = w1 = 0.0
        wlist = []

        for coefficient in coefficients[:-1]:
            w2 = w1
            w1 = w0
            w0 = coefficient + (s2 * w1 - w2)
            wlist.append(w1)

        components = coefficients[-1] + (s * w0 - w1)

        if scalar:
            components = components[:,0]

        yield components

        # Chebyshev differentiation.

        dw0 = dw1 = 0.0

        for coefficient, w1 in zip(coefficients[:-1], wlist):
            dw2 = dw1
            dw1 = dw0
            dw0 = 2.0 * w1 + dw1 * s2 - dw2

        rates = w0 + s * dw0 - dw1
        rates /= intlen
        rates *= 2.0
        rates *= S_PER_DAY

        if scalar:
            rates = rates[:,0]

        yield rates

class Type9Segment(BaseSegment):
    """Lagrange Interpolation - Unequal Time Steps"""

    def map_arrays(self):
        """Raw coefficients and epochs as memory-mapped NumPy arrays."""
        i = self.end_i
        polynomial_degree, number_of_states = self.daf.read_array(i - 1, i)
        if polynomial_degree != 1:
            raise ValueError('jplephem does not yet support Type 9 segments'
                             ' with a polynomial degree of {0}'
                             .format(polynomial_degree))
        number_of_states = int(number_of_states)
        i = self.start_i
        j = i + 6 * number_of_states - 1
        coefficients = self.daf.map_array(i, j)
        coefficients.shape = number_of_states, 6
        coefficients = coefficients.T
        epochs = self.daf.map_array(j + 1, j + number_of_states)
        return coefficients, epochs

    @reify
    def _data(self):
        """Cached arrays that are ready for interpolation."""
        coefficients, epochs = self.map_arrays()

        # Make iteration faster by pre-creating tuples of separate arrays.
        positions = tuple(coefficients[:3])
        and_velocities = tuple(coefficients)
        epochs = jd(epochs)
        return positions, and_velocities, epochs

    def compute(self, tdb, tdb2=0.0):
        """Interpolate [x y z] at time `tdb` plus `tdb2`.

        A standard JPL Type 9 ephemerides will return kilometers.

        """
        positions, and_velocities, epochs = self._data
        return array([interp(tdb, epochs, c) for c in positions])

    def compute_and_differentiate(self, tdb, tdb2=0.0):
        """Interpolate [x y z dx dy dz] at time `tdb` plus `tdb2`.

        A standard JPL Type 9 ephemerides will return kilometers and
        kilometers per second.

        """
        positions, and_velocities, epochs = self._data
        return array([interp(tdb, epochs, c) for c in and_velocities])

def titlecase(name):
    """Title-case target `name` if it looks safe to do so."""
    return name if name.startswith(('1', 'C/', 'DSS-')) else name.title()

_segment_classes = {
    2: Segment,
    3: Segment,
    9: Type9Segment,
}
