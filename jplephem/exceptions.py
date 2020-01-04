"""A set of special exceptions that can be thrown by the jplephem library"""

class OutOfRangeError(ValueError):
    """One or more time values given were out of range for the ephemeris.

    This exception is thrown if any input times are out of the range of
    times supported by an ephemeris.  It has an extra attribute:

    - `out_of_range_times`: if the input `tdb` of times is an array,
      this provides an array of booleans of the same length where `True`
      means the corresponding date is out of range.

    """
    def __init__(self, message, out_of_range_times):
        self.args = message,
        self.out_of_range_times = out_of_range_times
