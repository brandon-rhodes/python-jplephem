"""JPL Planetary and Lunar Ephemeris DE406 for the `jplephem` package.

This is the long-period ephemeris, published in May 1997, which is used
by the Jet Propulsion Laboratory's online HORIZONS service for dates far
in the past or future.  To limit its size to 190 MB, this ephemeris
allows greater error than `DE405 <http://pypi.python.org/pypi/de405>`_
amounting to up to 1m for the position of the Moon and 25m for all other
solar system bodies.

Coordinates and velocities are provided for the major planets, the Sun,
and Earth's Moon for the years -3000 through 3000.  To use this ephemeris
in Python, see the `jplephem <http://pypi.python.org/pypi/jplephem>`_
package on the Python Package Index.

"""

