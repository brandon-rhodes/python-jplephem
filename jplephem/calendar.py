"""Routines for dealing with Julian dates."""

def compute_calendar_date(jd_integer, julian_before=None):
    """Convert Julian day ``jd_integer`` to ``(year, month, day)``.

    Uses the proleptic Gregorian calendar unless ``julian_before`` is
    set to a specific Julian day, in which case the Julian calendar is
    used for dates older than that.

    """
    use_gregorian = (julian_before is None) or (jd_integer >= julian_before)

    # See the Explanatory Supplement to the Astronomical Almanac 15.11.
    f = jd_integer + 1401
    f += use_gregorian * ((4 * jd_integer + 274277) // 146097 * 3 // 4 - 38)
    e = 4 * f + 3
    g = e % 1461 // 4
    h = 5 * g + 2
    day = h % 153 // 5 + 1
    month = (h // 153 + 2) % 12 + 1
    year = e // 1461 - 4716 + (12 + 2 - month) // 12
    return year, month, day

def compute_julian_date(year, month=1, day=1.0):
    """Given a proleptic Gregorian date, return a Julian date float."""
    return compute_julian_day(year, month, day) - 0.5

def compute_julian_day(year, month=1, day=1):
    """Given a proleptic Gregorian date, return a Julian day int."""
    janfeb = month < 3
    return (+ 1461 * (year + 4800 - janfeb) // 4
            + 367 * (month - 2 + janfeb * 12) // 12
            - 3 * ((year + 4900 - janfeb) // 100) // 4
            - 32075
            + day)
