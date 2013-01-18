"""Use a JPL planetary ephemeris to predict planet positions.

This package uses a Jet Propulsion Laboratory ephemeris to predict the
position and velocity of a planet, or the magnitude and rate-of-change
of the Earth's nutation or the Moon's libration.  Its only dependency is
``NumPy``.  To take the smallest and most convenient ephemeris as an
example, you can install this package alongside ephemeris DE421 with
these commands::

    pip install jplephem
    pip install de421

Loading DE421 and computing a position require one line of Python each,
given a barycentric dynamical time expressed as a Julian date::

    import de421
    from jplephem import Ephemeris

    eph = Ephemeris(de421)
    x, y, z = eph.position('mars', 2444391.5)  # 1980.06.01

The result of calling ``position()`` is a 3-element NumPy array giving
the planet's position in the solar system in kilometers along the three
axes of the ICRF (a more precise reference frame than J2000 but oriented
in the same direction).  If you also want to know the planet's velocity,
call ``compute()`` instead::

    x, y, z, dx, dy, dz = eph.compute('mars', 2444391.5)

Velocities are returned as kilometers per day.

Both of these methods will also accept a NumPy array, which is the most
efficient way of computing a series of positions or velocities.  For
example, the position of Mars at each midnight over an entire year can
be computed with::

    import numpy as np
    t0 = 2444391.5
    t = np.arange(t0, t0 + 366.0, 1.0)
    x, y, z = eph.position('mars', 2444391.5)

You will find that ``x``, ``y``, and ``z`` in this case are each a NumPy
array of the same length as your input ``t``.

The string that you provide to ``e.compute()``, like ``'mars'`` in the
example above, actually names the data file that you want loaded from
the ephemeris package.  To see the list of data files that an ephemeris
provides, call its ``names()`` method.  Most of the ephemerides provide
thirteen data sets::

    earthmoon   mercury    pluto   venus
    jupiter     moon       saturn
    librations  neptune    sun
    mars        nutations  uranus

Each ephemeris covers a specific range of dates, beyond which it cannot
provide reliable predictions of each planet's position.  These limits
are available as attributes of the ephemeris::

    t0, t1 = eph.jalpha, eph.jomega

The ephemerides currently available as Python packages (the following
links explain the differences between them) are:

* `DE405 <http://pypi.python.org/pypi/de405>`_ (May 1997)
* `DE406 <http://pypi.python.org/pypi/de406>`_ (May 1997)
* `DE421 <http://pypi.python.org/pypi/de421>`_ (February 2008)
* `DE422 <http://pypi.python.org/pypi/de422>`_ (September 2009)
* `DE423 <http://pypi.python.org/pypi/de423>`_ (February 2010)

"""
from .ephem import Ephemeris, DateError
