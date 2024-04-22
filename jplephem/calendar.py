"""Routines for dealing with Julian dates."""

def compute_julian_day(year, month=1, day=1):
    """Given a proleptic Gregorian calendar date, return a Julian day int."""
    janfeb = month < 3
    return (day
            + 1461 * (year + 4800 - janfeb) // 4
            + 367 * (month - 2 + janfeb * 12) // 12
            - 3 * ((year + 4900 - janfeb) // 100) // 4
            - 32075)
