"""JPL Planetary and Lunar Ephemeris DE421 for the jplephem package.

This is a recent short-period ephemeris published by the Jet Propulsion
Laboratory.  It requires only 27 MB of storage and is specially accurate
with respect to the position of Earth's Moon.

:Name: DE421 (February 2008)
:Years: 1900 through 2050
:Planets: Yes
:Sun/Moon: Yes
:Nutations: Yes
:Librations: Yes
:Report: `Folkner, Williams, Boggs (2009) [PDF] <http://tmo.jpl.nasa.gov/progress_report/42-178/178C.pdf>`_
:Size: 27 MB

The JPL called this ephemeris is a "significant advance" over
predecessors like `DE405 <http://pypi.python.org/pypi/de405>`_ / `DE406
<http://pypi.python.org/pypi/de406>`_ and cited accuracies that are in
many cases ten times greater, such as giving the position of Venus
within 200m and the positions of Earth and Mars within 300m over the
last decade.  Note that even greater accuracies are achieved, for
Mercury and Venus in particular (but not for the Moon), by `DE423
<http://pypi.python.org/pypi/de423>`_ which also has the advantage of
covering a 400-year period instead of only 150 years.

Greater accuracy can also be expected from the long-term ephemeris
`DE422 <http://pypi.python.org/pypi/de422>`_ since it incorporates more
spacecraft observations than DE421.  It also covers a period of 6000
years, making it useful to astronomy historians.  But as it requires a
half-gigabyte of disk space, some users may prefer DE421.

To compute using this ephemeris in Python, see the `jplephem
<http://pypi.python.org/pypi/jplephem>`_ package.

"""
