"""Package supporting JPL planetary ephemeris computations.

This package lets you consult a Jet Propulsion Laboratory ephemeris for
the position and velocity of one of the planets, or the magnitude and
rate-of-change of the Earth's nutation or the Moon's libration.  To
determine the position of Mars using the DE423 ephemeris, for example,
you would start by installing two packages::

    pip install jplephem
    pip install de423

Then you can compute positions using a script like this::

    import de423
    from jplephem import Ephemeris

    e = Ephemeris(de423)
    jed = 2444391.5  # 1980.06.01
    print e.compute('mars', jed)

The result should be a tuple providing the object's position in the
Solar System given in kilometers, as well as its velocity in kilometers
per second::

    (x, y, z, xrate, yrate, zrate)

The ephemerides currently available as Python packages (the following
links explain the differences between them) are:

* `DE405 <http://pypi.python.org/pypi/de421>`_ (May 1997)
* `DE406 <http://pypi.python.org/pypi/de421>`_ (May 1997)
* `DE421 <http://pypi.python.org/pypi/de421>`_ (February 2008)
* `DE422 <http://pypi.python.org/pypi/de421>`_ (September 2009)
* `DE423 <http://pypi.python.org/pypi/de421>`_ (February 2010)

"""
from .ephem import Ephemeris
