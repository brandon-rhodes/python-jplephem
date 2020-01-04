"""A set of special exceptions that can be thrown by the jplephem library"""

class OutOfRangeTimestampError(ValueError):
    """
    This exception is thrown if any input times are out of the range of
    times jplephem can compute ephemeris for.
    It has for properties:

    - `message` is a string explaining what happened,
    - `min_timestamp` and `max_timestamp` are floats giving the minimum and
      maximum supported times,
    - `out_of_range_times` is an array of booleans where `True` means the
      corresponding date in the input array is out of range and `False` means
      it is correct.
    """
    def __init__(self, message, min_timestamp, max_timestamp, out_of_range_times):
        self.message = message
        self.min_timestamp = min_timestamp
        self.max_timestamp = max_timestamp
        self.out_of_range_times = out_of_range_times
