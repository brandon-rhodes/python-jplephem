"""JPL Planetary and Lunar Ephemeris DE405 for the jplephem package.

This ephemeris has been the basis for the Astronomical Almanac of the
United States Naval Observatory since 2003 and powers the Jet Propulsion
Laboratory's online HORIZONS service. Coordinates and velocities are
provided for the major planets, the Sun, and Earth's Moon.

:Name: DE405 (May 1997)
:Years: 1600 through 2200
:Planets: Yes
:Sun/Moon: Yes
:Nutations: Yes
:Librations: Yes
:Report: `Standish (1998) [PDF] <ftp://ssd.jpl.nasa.gov/pub/eph/planets/ioms/de405.iom.pdf>`_
:Size: 54 MB

The accuracy of this ephemeris is good; the position of the Earth or
Mars, for example, should be accurate to within 2 km.  But the more
recent `DE421 <http://pypi.python.org/pypi/de421>`_ ephemeris provides
greater accuracy, especially with respect to the Moon, and you should
use it instead if you are planning a space mission.  For missions to
Mercury or Venus, `DE423 <http://pypi.python.org/pypi/de421>`_ will be
an even better choice.

See `DE406 <http://pypi.python.org/pypi/de406>`_ if you are interested
in a similar ephemeris for dates far in the past or future, or `DE422
<http://pypi.python.org/pypi/de422>`_ if you want high accuracy over a
long time period (and have enough disk space).

To compute using this ephemeris in Python, see the `jplephem
<http://pypi.python.org/pypi/jplephem>`_ package.

"""
