# -*- encoding: utf-8 -*-

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
call ``position_and_velocity()`` instead::

    position, velocity = eph.position_and_velocity('mars', 2444391.5)
    x, y, z = position            # a NumPy array
    xdot, ydot, zdot = velocity   # another array

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
provides, consult its ``names`` attribute.  Most of the JPL ephemerides
provide thirteen data sets::

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

High-Precision Dates
--------------------

Since all modern Julian dates are numbers larger than 2.4 million, a
standard 64-bit Python or NumPy float necessarily leaves only a limited
number of bits available for the fractional part.  An unpublished paper
(“AA Technical Note 2011-02”) by the United States Naval Observatory's
Astronomical Applications Department suggests that the `precision
possible with a 64-bit floating point Julian date is around 20.1µs
<http://ad.usno.navy.mil/edboard/110308.txt>`_.

If you need to supply times and receive back planetary positions with
greater precision than 20.1µs, then you have two options.

First, you can supply times using the special `float96` NumPy type,
which is also aliased to the name `longfloat`.  If you provide either a
`float96` or a `float96` array as your `TB` parameter to any `jplephem`
routine, you should get back a high-precision result.

Second, you can split your date or dates into two pieces, and supply
them as a pair of arguments two `tdb` and `tdb2`; one popular approach
for how to split your date is to use the `tdb` float for the integer
Julian date, and `tdb2` for the fraction that specifies the time of day.
Nearly all `jplephem` routines accept this optional `tdb2` argument if
you wish to provide it, thanks to the work of Marten van Kerkwijk!

Waiting To Compute Velocity
---------------------------

When a high-level astronomy library computes the distance between an
observer and a solar system body, it typically measures the light travel
delay to body and then uses a loop to take several steps backwards in
time, because a human with a telescope will not see the planet as it is
right now.  Instead, the planet will appear at the position it occupied
when the light that is now reaching Earth left the planet's surface or
clouds.

To make such loops less computationally expensive — loops that only need
to compute the planet position repeatedly, where computing the velocity
can wait until the loop's conclusion — `jplephem` provides a way to
split the `position_and_velocity()` call into two pieces.  This lets you
examine the position *before* deciding whether to proceed with also
with the expense of computing the velocity.

The key is the special `compute_bundle()` method, with returns a tuple
containing the coefficients and intermediate results that are needed by
*both* the position and the velocity computations.  There is nothing
wasted in calling `compute_bundle()` whether you are going to ask for
the position, the velocity, or both as your next computing step!

So your loop can look something like this::

    bundle = eph.compute_bundle('mars', tdb)
    position = eph.position_from_bundle(bundle)

    while True:
        # ...determine whether you are happy...
        if you_are_happy:
            break
        # otherwise, adjust `tdb` and re-compute:
        bundle = eph.compute_bundle('mars', tdb)
        position = eph.position_from_bundle(bundle)

    # Now we re-use the values in `bundle`, for free!
    velocity = eph.velocity_from_bundle(bundle)

This is especially important when the number of dates in `tdb` is large,
since even NumPy vector operations over hundreds of thousands of values
is going to take a noticeable amount of time, and every mass operation
that is avoided will help shepherd your program toward completion.

"""
from .ephem import Ephemeris, DateError

__all__ = ['Ephemeris', 'DateError']
