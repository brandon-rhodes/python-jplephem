"""Compute positions from a NASA SPICE SPK ephemeris kernel file.

http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html

"""
from numpy import array, interp, rollaxis, zeros, reshape
from .calendar import compute_calendar_date
from .daf import DAF
from .descriptorlib import reify
from .exceptions import OutOfRangeError
from .names import target_names

T0 = 2451545.0
S_PER_DAY = 86400.0

def _jd(seconds):
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
        f = open(path, 'rb')
        try:
            return cls(DAF(f))
        except Exception:
            f.close()
            raise

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
        self.start_jd = _jd(self.start_second)
        self.end_jd = _jd(self.end_second)

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
        initial_epoch = _jd(init)
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
        epochs = _jd(epochs)
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

class Type1Segment(BaseSegment):
    """Type 1: Modified Difference Arrays

    Original Repository: https://github.com/whiskie14142/spktype01
    Original Author: Shushi Uetsuki (whiskie14142)
    Modified: xiaozhongguo@gmail.com
    """

    def __init__(self, daf, source, descriptor):
        super().__init__(daf, source, descriptor)

        # initialize arrays for spke01
        self.G = zeros(15)
        self.REFPOS = zeros(3)
        self.REFVEL = zeros(3)
        self.KQ = array([0, 0, 0])
        self.FC = zeros(15)
        self.FC[0] = 1.0
        self.WC = zeros(13)
        self.W = zeros(17)
        
        # initialize for compute_type01
        self.mda_record_exist = False

    def _get_MDA_record(self, time_sec):
        """Return a Modified Difference Array(MDA) record for the time to 
        evaluate with its effective time boundaries (lower and upper).
        Inputs:
            time_sec - epoch for computation, seconds from J2000
        Returns: mda_record, lower_boundary, upper_boundary
            mda_record: A Modified Difference Array record
            lower_boundary: lower boundary of the record, seconds since J2000
            upper_boundary: upper boundary of the record, seconds since J2000
        """

        # Number of records in this segment
        entry_count = int(self.daf.map_array(self.end_i, self.end_i))
        
        # Number of entries in epoch directory 
        epoch_dir_count = entry_count // 100
        
        # serch target epoch in epoch directory to narrow serching aria
        if epoch_dir_count >= 1:
            epoch_dir = self.daf.map_array(self.end_i - epoch_dir_count,
                                            self.end_i - 1)
            found = False
            for i in range(1, epoch_dir_count + 1):
                if epoch_dir[i-1] > time_sec:
                    found = True
                    break
            if found:
                serch_last_index = i * 100
                serch_start_index = (i - 1) * 100 + 1
            else:
                serch_last_index = entry_count
                serch_start_index = epoch_dir_count * 100 + 1
        else:
            serch_last_index = entry_count
            serch_start_index = 1

        # epoch_table contains epochs for all records in this segment        
        epoch_table = self.daf.map_array(self.start_i + (entry_count * 71),
                                       self.start_i + (entry_count * 71) + entry_count - 1)

        # serch target epoch in epoch_table
        found = False
        for i in range(serch_start_index, serch_last_index + 1):
            if epoch_table[i-1] > time_sec:
                found = True
                break
        if not found:
            i = serch_last_index
        record_index = i
        upper_boundary = epoch_table[i-1]
        if i != 1:
            lower_boundary = epoch_table[i-2]
        else:
            lower_boundary = self.start_second
        
        mda_record = self.daf.map_array(self.start_i + ((record_index - 1) * 71),
                                        self.start_i + (record_index * 71) - 1)

        # mda_record : one record of MDA
        # lower_boundary : lower boundary of epoch in this MDA record
        # upper_boundary : upper boundary of epoch in this MDA record
        return mda_record, lower_boundary, upper_boundary

    def get_MDA_record(self, eval_sec):
        """Return a MDA record for defined epoch.
        Inputs:
            eval_sec - epoch for computation, seconds from J2000
        Returns:
            MDA record - a Numpy array of 71 floating point numbers
        Exception:
            ValueError will be raised when:
                eval_sec is outside of SPK data
        """
        
        # chech last segment can be used
        if eval_sec < self.start_second or eval_sec > self.end_second:
            raise ValueError('Invalid Time to evaluate')
        
        return self._get_MDA_record(eval_sec)
        
    def spke01(self, ET, RECORD):
        """Compute position and velocity from a Modified Difference Array record
        
        Inputs:
            ET: Epoch time to evaluate position and velocity (seconds since J2000)
            RECORD: A record of Modified Difference Array
        Returns: STATE
            STATE: A numpy array which contains position and velocity
        """
        
        # This method has been translated from SPKE01 of SPICE Toolkit and
        # modified by Shushi Uetsuki.
        #
        # SPICE Toolkit for FORTRAN : http://naif.jpl.nasa.gov/naif/toolkit_FORTRAN.html
        # SPK Required Reading : http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/C/req/spk.html
        #
        # Original FORTRAN code uses 'SAVE' directive, and it means all variable
        # should be saved for next call.  So i decided to make almost all 
        # variables to be instance variable.  Some of them are initialized in 
        # __init__ method.
        
        STATE = zeros(6)

#     Variable  I/O  Description
#     --------  ---  --------------------------------------------------
#     ET         I   Target epoch.
#     RECORD     I   Data record.
#     STATE      O   State (position and velocity).
#
#$ Detailed_Input
#
#     ET          is a target epoch, at which a state vector is to
#                 be computed.
#
#     RECORD      is a data record which, when evaluated at epoch ET,
#                 will give the state (position and velocity) of some
#                 body, relative to some center, in some inertial
#                 reference frame.
#
#$ Detailed_Output
#
#     STATE       is the state. Units are km and km/sec.
#
#$ Parameters
#
#     None.
#
#$ Exceptions
#
#     None.
#
#$ Files
#
#     None.
#
#$ Particulars
#
#     The exact format and structure of type 1 (difference lines)
#     segments are described in the SPK Required Reading file.
#
#     Difference lines (DL's) are generated by JPL navigation
#     system programs P and PV. Each data record is equivalent
#     to the (slightly rearranged) 'P' portion of a NAVIO PV file
#     data record.
#
#     SPKE01 is a specialized version of Fred Krogh's subroutine DAINT.
#     Only the calling sequence has been changed.
#
#     Because the original version was undocumented, only Fred 
#     knows how this really works.
#
#$ Examples
#
#     None.
#
#$ Restrictions
#
#     Unknown.
#
#$ Literature_References
#
#     NAIF Document 168.0, "S- and P- Kernel (SPK) Specification and
#     User's Guide"
#
#$ Author_and_Institution
#
#     F.T. Krogh      (JPL)
#     I.M. Underwood  (JPL)
#
#$ Version
#
#-    SPICELIB Version 1.1.0, 14-FEB-1997 (WLT) 
#     
#        The goto's were removed and loop and if structures 
#        revealed.  We still don't know exactly what's going 
#        on, but at least the bones of this routine have been 
#        cleaned off and are ready for assembly. (WLT)
#
#-    SPICELIB Version 1.0.4, 30-OCT-1996 (WLT)
#
#        Removed redundant SAVE statements from the declaration
#        section.  Thanks to Steve Schlaifer for finding this
#        error.
#
#-    SPICELIB Version 1.0.3, 10-MAR-1992 (WLT)
#
#        Comment section for permuted index source lines was added
#        following the header.
#
#-    SPICELIB Version 1.0.2, 23-AUG-1991 (HAN)
#
#        SPK01 was removed from the Required_Reading section of the
#        header. The information in the SPK01 Required Reading file
#        is now part of the SPK Required Reading file.
#
#-    SPICELIB Version 1.0.1, 22-MAR-1990 (HAN)
#
#        Literature references added to the header.
#
#-    SPICELIB Version 1.0.0, 31-JAN-1990 (IMU) (FTK)
#
#-&
# 
#$ Index_Entries
#
#     evaluate type_1 spk segment
#
#-&


        
#
#     Unpack the contents of the MDA array.
#
#        Name    Dimension  Description
#        ------  ---------  -------------------------------
#        TL              1  Final epoch of record
#        G              15  Stepsize function vector
#        REFPOS          3  Reference position vector
#        REFVEL          3  Reference velocity vector
#        DT         15,NTE  Modified divided difference arrays
#        KQMAX1          1  Maximum integration order plus 1
#        KQ            NTE  Integration order array
#
#     For our purposes, NTE is always 3.
#
        self.TL = RECORD[0]
        self.G = RECORD[1:1+15]
#     
#     Collect the reference position and velocity.
#     
        self.REFPOS[0] = RECORD[16]
        self.REFVEL[0] = RECORD[17]
        
        self.REFPOS[1] = RECORD[18]
        self.REFVEL[1] = RECORD[19]
        
        self.REFPOS[2] = RECORD[20]
        self.REFVEL[2] = RECORD[21]
        
        self.DT = reshape(RECORD[22:22+45], (15, 3), order='F')
        
        self.KQMAX1 = int(RECORD[67])
        self.KQ[0] = int(RECORD[68])
        self.KQ[1] = int(RECORD[69])
        self.KQ[2] = int(RECORD[70])
#     
#     Next we set up for the computation of the various differences
#     
        self.DELTA = ET - self.TL
        self.TP = self.DELTA
        self.MQ2 = self.KQMAX1 - 2
        self.KS = self.KQMAX1 - 1
#
#     This is clearly collecting some kind of coefficients.  
#     The problem is that we have no idea what they are...
#     
#     The G coefficients are supposed to be some kind of step size 
#     vector. 
#     
#     TP starts out as the delta t between the request time 
#     and the time for which we last had a state in the MDL file. 
#     We then change it from DELTA  by the components of the stepsize 
#     vector G.  
#
        for J in range(1, self.MQ2 + 1):
            self.FC[J] = self.TP / self.G[J-1]
            self.WC[J-1] = self.DELTA / self.G[J-1]
            self.TP = self.DELTA + self.G[J-1]
#
#     Collect KQMAX1 reciprocals. 
#   
        for J in range(1, self.KQMAX1 + 1):
            self.W[J-1] = 1.0 / float(J)
#
#     Compute the W(K) terms needed for the position interpolation
#     (Note,  it is assumed throughout this routine that KS, which 
#     starts out as KQMAX1-1 (the ``maximum integration'') 
#     is at least 2.
#
        self.JX = 0
        self.KS1 = self.KS - 1
        
        while self.KS >= 2:
            
            self.JX = self.JX + 1
            
            for J in range(1, self.JX + 1):
                self.W[J+self.KS-1] = self.FC[J] * self.W[J+self.KS1-1] - self.WC[J-1] * self.W[J+self.KS-1]
            
            self.KS = self.KS1
            self.KS1 = self.KS1 - 1
#
#     Perform position interpolation: (Note that KS = 1 right now.
#     We don't know much more than that.)
#
        for I in range(1, 3 + 1):
            
            self.KQQ = self.KQ[I-1]
            self.SUM = 0.0
            
            for J in range(self.KQQ, 0, -1):
                self.SUM = self.SUM + self.DT[J-1, I-1] * self.W[J+self.KS-1]
            
            STATE[I-1] = self.REFPOS[I-1] + self.DELTA * (self.REFVEL[I-1] + self.DELTA * self.SUM)
#
#     Again we need to compute the W(K) coefficients that are 
#     going to be used in the velocity interpolation. 
#     (Note, at this point, KS = 1, KS1 = 0.)
#      
        for J in range(1, self.JX + 1):
            self.W[J+self.KS-1] = self.FC[J] * self.W[J+self.KS1-1] - self.WC[J-1] * self.W[J+self.KS-1]
        
        self.KS = self.KS - 1
        
#
#     Perform velocity interpolation:
#
        for I in range(1, 3 + 1):
            self.KQQ = self.KQ[I-1]
            self.SUM = 0.0
            
            for J in range(self.KQQ, 0, -1):
                self.SUM = self.SUM + self.DT[J-1, I-1] * self.W[J+self.KS-1]
            
            STATE[I+3-1] = self.REFVEL[I-1] + self.DELTA * self.SUM
        
#
#     That's all folks.  We don't know why we did anything, but 
#     at least we can tell structurally what we did.
#      
        
        return STATE

    def compute_and_differentiate(self, tdb, tdb2=0.0):
        """Compute position and velocity of target from SPK data (data type 1).
        Inputs:
            tdb, tdb2 - Julian date of epoch for computation.  (tdb + tdb2) will 
                be used for computation.  If you want precise definition of 
                epoch, tdb should be an integer or a half integer, and tdb2
                should be a relatively small floating point number.
        Returns:
            Position (X, Y, Z) and velocity (XD, YD, ZD) of the target at 
            epoch.  Position and velocity are provided as Numpy arrays 
            respectively.
        """
        scalar = not getattr(tdb, 'shape', 0) and not getattr(tdb2, 'shape', 0)

        eval_sec = (tdb - T0)
        eval_sec = (eval_sec + tdb2) * S_PER_DAY

        if scalar:
            eval_sec = array([eval_sec])

        results = zeros((6, eval_sec.size))

        for i, eval_seci in enumerate(eval_sec):
            if self.mda_record_exist:
                if eval_seci >= self.mda_lb and eval_seci < self.mda_ub:
                    result = self.spke01(eval_seci, self.mda_record)
                    results[:, i] = result
                    continue

            self.mda_record, self.mda_lb, self.mda_ub = self.get_MDA_record(eval_seci)
            self.mda_record_exists = True

            result = self.spke01(eval_seci, self.mda_record)
            results[:, i] = result
        if scalar:
            return results[0:3, 0], results[3:, 0]
        return results[0:3], results[3:]

    def compute(self, tdb, tdb2=0.0):
        """Compute position and velocity of target from SPK data (data type 1)."""
        position, _ = self.compute_and_differentiate(tdb, tdb2)
        return position

def titlecase(name):
    """Title-case target `name` if it looks safe to do so."""
    return name if name.startswith(('1', 'C/', 'DSS-')) else name.title()

_segment_classes = {
    1: Type1Segment,
    2: Segment,
    3: Segment,
    9: Type9Segment,
}
