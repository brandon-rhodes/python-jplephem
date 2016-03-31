"""JPL Planetary and Lunar Ephemeris DE423 for the jplephem package.

This is the most recent short-period ephemeris published by the Jet
Propulsion Laboratory.  It achieves especially high accuracy for Mercury
and Venus since it was prepared for the MESSENGER mission, and should be
able to predict the coordinates of those planets within a few
milliarcseconds.

:Name: DE423 (February 2010)
:Years: 1800 through 2200
:Planets: Yes
:Sun/Moon: Yes
:Nutations: Yes
:Librations: Yes
:Report: `Folkner (2010) [PDF] <ftp://ssd.jpl.nasa.gov/pub/eph/planets/ioms/de423.iom.pdf>`_
:Size: 36 MB

While lunar missions will want to use the slightly older ephemeris
`DE421 <http://pypi.python.org/pypi/de421>`_ because of the extra
accuracy that it provides for the Moon, general purpose users will
probably find this ephemeris more useful, especially as it covers a full
400 years instead of only 150 years.

Similar accuracy can be expected from the long-term ephemeris `DE422
<http://pypi.python.org/pypi/de422>`_ which also covers a period of 6000
years, making it useful to astronomy historians.  But as it requires a
half-gigabyte of disk space, some users may prefer DE423.

To compute using this ephemeris in Python, see the `jplephem
<http://pypi.python.org/pypi/jplephem>`_ package.

"""
