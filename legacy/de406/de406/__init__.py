"""JPL Planetary and Lunar Ephemeris DE406 for the jplephem package.

This is the long-period ephemeris which is used by the Jet Propulsion
Laboratory's online HORIZONS service for dates far in the past or
future.

:Name: DE406 (May 1997)
:Years: -3000 through 3000
:Planets: Yes
:Sun/Moon: Yes
:Nutations: No
:Librations: No
:Report: `Standish (1998) [PDF] <http://iau-comm4.jpl.nasa.gov/de405iom/de405iom.pdf>`_
:Size: 190 MB

This ephemeris is like `DE405 <http://pypi.python.org/pypi/de405>`_ but
covers a longer time period.  To keep its files from being too large,
its size was reduced by permitting up to 1 meter of interpolation error
for the position of the Moon, and 25 meters for all other solar system
bodies.  Since DE405 itself is often only accurate to within a few
kilometers for planetary positions, the difference was not important for
many users.  You can get much higher accuracy from the more recent
long-term ephemeris `DE422 <http://pypi.python.org/pypi/de422>`_ but at
the cost of three times the RAM and disk space.

To compute using this ephemeris in Python, see the `jplephem
<http://pypi.python.org/pypi/jplephem>`_ package.

"""
