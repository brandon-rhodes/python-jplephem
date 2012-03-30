"""JPL Planetary and Lunar Ephemeris DE405 for the `jplephem` package.

This ephemeris, published in May 1997, has been the basis for the
Astronomical Almanac of the United States Naval Observatory since 2003,
and powers the Jet Propulsion Laboratory's online HORIZONS service.
While the more recent DE422 ephemeris (September 2009) provides even
greater accuracy, it requires ten times the disk space and memory; this
ephemeris requires only 54 MB of storage.

Positions computed using this ephemeris should be accurate to within
roughly 0.001 arcseconds for the inner planets, and 0.1 arcseconds for
the outer planets.

Coordinates and velocities are provided for the major planets, the Sun,
and Earth's Moon for the years 1600 through 2200.  To use this ephemeris
in Python, see the `jplephem <http://pypi.python.org/pypi/jplephem>`_
package on the Python Package Index.

"""
