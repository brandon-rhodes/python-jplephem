"""Package supporting JPL planetary ephemeris computations.

This package lets you consult a Jet Propulsion Laboratory ephemeris for
the position and velocity of one of the planets, or the magnitude and
rate-of-change of the Earth's nutation or the Moon's libration.  To
determine the position of Mars using the DE421 ephemeris, for example,
you would start by installing two packages::

    pip install jplephem
    pip install de421

Then you can compute positions using a script like this::

    import de421
    from jplephem import Ephemeris

    eph = Ephemeris(de421)
    jd = 2444391.5  # 1980.06.01
    print eph.compute('mars', jd)

The result will be a 6-element NumPy array providing the object's
position in the Solar System, given in kilometers along the axes of the
ICRF (a more precise reference frame than J2000 but oriented in the same
direction), as well as its velocity along those axes in kilometers per
day::

    (x, y, z, xrate, yrate, zrate)

The string that you provide to ``e.compute()``, like ``mars`` in the
example above, actually names the data file that you want loaded and
used from the ephemeris package.  To see the full list of data files
that an ephemeris provides, you can simply list the files in its
directory.  Most of the ephemerides provide thirteen data sets::

    earthmoon   mercury    pluto   venus
    jupiter     moon       saturn
    librations  neptune    sun
    mars        nutations  uranus

The ephemerides currently available as Python packages (the following
links explain the differences between them) are:

* `DE405 <http://pypi.python.org/pypi/de405>`_ (May 1997)
* `DE406 <http://pypi.python.org/pypi/de406>`_ (May 1997)
* `DE421 <http://pypi.python.org/pypi/de421>`_ (February 2008)
* `DE422 <http://pypi.python.org/pypi/de422>`_ (September 2009)
* `DE423 <http://pypi.python.org/pypi/de423>`_ (February 2010)

"""
from .ephem import Ephemeris, DateError
