# -*- encoding: utf-8 -*-

"""Use a JPL ephemeris to predict planet positions.

This package can load and use a Jet Propulsion Laboratory (JPL)
ephemeris for predicting the position and velocity of a planet or other
Solar System body.  It only needs `NumPy <http://www.numpy.org/>`_,
which ``pip`` will automatically attempt to install alongside
``pyephem`` when you run::

    $ pip install jplephem

If you see NumPy compilation errors, then try downloading and installing
NumPy directly from `its web site <http://www.numpy.org/>`_ or simply
use a distribution of Python with science tools already installed, like
`Anaconda <http://continuum.io/downloads>`_.

Note that ``jplephem`` offers only the logic necessary to produce plain
three-dimensional vectors.  Most programmers interested in astronomy
will want to look at `Skyfield <http://rhodesmill.org/skyfield/>`_
instead, which uses ``jplephem`` but converts the numbers into more
traditional measurements like right ascension and declination.

Most users will use ``jplephem`` with the Satellite Planet Kernel (SPK)
files that the NAIF facility at NASA JPL offers for use with their own
SPICE toolkit.  They have collected their most useful kernels beneath
the directory:

http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/

To learn more about SPK files, the official `SPK Required Reading
<http://naif.jpl.nasa.gov/pub/naif/toolkit_docs/FORTRAN/req/spk.html>`_
document is available from the NAIF facility’s web site under the NASA
JPL domain.


Command Line Tool
-----------------

If you have downloaded a ``.bsp`` file and want to learn what ephemeris
segments are stored inside of it, you can have ``jplephem`` print them
out by invoking the module directly from the command line::

    python -m jplephem de430.bsp

This will print out a summary identical to the one shown in the
following section, but without requiring that you type and run any
Python code.


Getting Started With DE430
--------------------------

The recent DE430 ephemeris is a useful starting point.  It weighs in at
115 MB, but provides predictions across the generous range of years
1550–2650:

http://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/de430.bsp

After the kernel has downloaded, you can use ``jplephem`` to load this
SPK file and learn about the segments it offers:

>>> from jplephem.spk import SPK
>>> kernel = SPK.open('de430.bsp')
>>> print(kernel)
File type DAF/SPK and format LTL-IEEE with 14 segments:
2287184.50..2688976.50  Solar System Barycenter (0) -> Mercury Barycenter (1)
2287184.50..2688976.50  Solar System Barycenter (0) -> Venus Barycenter (2)
2287184.50..2688976.50  Solar System Barycenter (0) -> Earth Barycenter (3)
2287184.50..2688976.50  Solar System Barycenter (0) -> Mars Barycenter (4)
2287184.50..2688976.50  Solar System Barycenter (0) -> Jupiter Barycenter (5)
2287184.50..2688976.50  Solar System Barycenter (0) -> Saturn Barycenter (6)
2287184.50..2688976.50  Solar System Barycenter (0) -> Uranus Barycenter (7)
2287184.50..2688976.50  Solar System Barycenter (0) -> Neptune Barycenter (8)
2287184.50..2688976.50  Solar System Barycenter (0) -> Pluto Barycenter (9)
2287184.50..2688976.50  Solar System Barycenter (0) -> Sun (10)
2287184.50..2688976.50  Earth Barycenter (3) -> Moon (301)
2287184.50..2688976.50  Earth Barycenter (3) -> Earth (399)
2287184.50..2688976.50  Mercury Barycenter (1) -> Mercury (199)
2287184.50..2688976.50  Venus Barycenter (2) -> Venus (299)

Each segment of the file lets you predict the position of an object with
respect to some other reference point.  If you want the coordinates of
Mars at 2457061.5 (2015 February 8) with respect to the center of the
solar system, this ephemeris only requires you to take a single step:

>>> position = kernel[0,4].compute(2457061.5)
>>> print(position)
[  2.05700211e+08   4.25141646e+07   1.39379183e+07]

But learning the position of Mars with respect to the Earth takes three
steps, from Mars to the Solar System barycenter to the Earth-Moon
barycenter and finally to Earth itself:

>>> position = kernel[0,4].compute(2457061.5)
>>> position -= kernel[0,3].compute(2457061.5)
>>> position -= kernel[3,399].compute(2457061.5)
>>> print(position)
[  3.16065185e+08  -4.67929557e+07  -2.47554111e+07]

You can see that the output of this ephemeris is in kilometers.  If you
use another ephemeris, check its documentation to be sure of the units
that it employs.

If you supply the date as a NumPy array, then each component that is
returned will itself be a vector as long as your date:

>>> import numpy as np
>>> jd = np.array([2457061.5, 2457062.5, 2457063.5, 2457064.5])
>>> position = kernel[0,4].compute(jd)
>>> print(position)
[[  2.05700211e+08   2.05325363e+08   2.04928663e+08   2.04510189e+08]
 [  4.25141646e+07   4.45315179e+07   4.65441136e+07   4.85517457e+07]
 [  1.39379183e+07   1.48733243e+07   1.58071381e+07   1.67392630e+07]]

Some ephemerides include velocity inline by returning a 6-vector instead
of a 3-vector.  For an ephemeris that does not, you can ask for the
Chebyshev polynomial to be differentiated to produce a velocity, which
is delivered as a second return value:

>>> position, velocity = kernel[0,4].compute_and_differentiate(2457061.5)
>>> print(position)
[  2.05700211e+08   4.25141646e+07   1.39379183e+07]
>>> print(velocity)
[ -363896.06287889  2019662.99596519   936169.77271558]


Details of the API
------------------

Here are a few details for people ready to go beyond the high-level API
provided above and read through the code to learn more.

* Instead of reading an entire ephemeris into memory, ``jplephem``
  memory-maps the underlying file so that the operating system can
  efficiently page into RAM only the data that your code is using.

* Once the metadata has been parsed from the binary SPK file, the
  polynomial coefficients themselves are loaded by building a NumPy
  array object that has access to the raw binary file contents.
  Happily, NumPy already knows how to interpret a packed array of
  double-precision floats.  You can learn about the underlying DAF
  “Double Precision Array File” format, in case you ever need to open
  other such array files in Python, through the ``DAF`` class in the
  module ``jplephem.daf``.

* An SPK file is made of segments.  When you first create an ``SPK``
  kernel object ``k``, it examines the file and creates a list of
  ``Segment`` objects that it keeps in a list under an attribute named
  ``k.segments`` which you are free to examine in your own code by
  looping over it.

* There is more information about each segment beyond the one-line
  summary that you get when you print out the SPK file, which you can
  see by asking the segment to print itself verbosely:

  >>> segment = kernel[3,399]
  >>> print(segment.describe())
  2287184.50..2688976.50  Earth Barycenter (3) -> Earth (399)
    frame=1 data_type=2 source=DE-0430LE-0430

* Each ``Segment`` loaded from the kernel has a number of attributes
  that are loaded from the SPK file:

  >>> help(segment)
  Help on Segment in module jplephem.spk object:
  ...
   |  segment.source - official ephemeris name, like 'DE-0430LE-0430'
   |  segment.start_second - initial epoch, as seconds from J2000
   |  segment.end_second - final epoch, as seconds from J2000
   |  segment.start_jd - start_second, converted to a Julian Date
   |  segment.end_jd - end_second, converted to a Julian Date
   |  segment.center - integer center identifier
   |  segment.target - integer target identifier
   |  segment.frame - integer frame identifier
   |  segment.data_type - integer data type identifier
   |  segment.start_i - index where segment starts
   |  segment.end_i - index where segment ends
  ...

* The square-bracket lookup mechanism ``kernel[3,399]`` is a
  non-standard convenience that returns only the last matching segment
  in the file.  While the SPK standard does say that the last segment
  takes precedence, it also says that earlier segments for a particular
  center-target pair should be fallen back upon for dates that the last
  segment does not cover.  So, if you ever tackle a complicated kernel,
  you will need to implement fallback rules that send some dates to the
  final segment for a given center and target, but that send other dates
  to earlier segments that are qualified to cover them.

* If you are accounting for light travel time and require repeated
  computation of the position, but then need the velocity at the end,
  and want to avoid repeating the expensive position calculation, then
  try out the ``segment.generate()`` method - it will let you ask for
  the position, and then only proceed to the velocity once you are sure
  that the light-time error is now small enough.


High-Precision Dates
--------------------

Since all modern Julian dates are numbers larger than 2.4 million, a
standard 64-bit Python or NumPy float necessarily leaves only a limited
number of bits available for the fractional part.  *Technical Note
2011-02* from the United States Naval Observatory's Astronomical
Applications Department suggests that the `precision possible with a
64-bit floating point Julian date is around 20.1 µs
<http://jplephem.s3.amazonaws.com/JD_precision_test.pdf>`_.

If you need to supply times and receive back planetary positions with
greater precision than 20.1 µs, then you have two options.

First, you can supply times using the special ``float96`` NumPy type,
which is also aliased to the name ``longfloat``.  If you provide either
a ``float96`` scalar or a ``float96`` array as your ``tdb`` parameter to
any ``jplephem`` routine, you should get back a high-precision result.

Second, you can split your date or dates into two pieces, and supply
them as a pair of arguments two ``tdb`` and ``tdb2``.  One popular
approach for how to split your date is to use the ``tdb`` float for the
integer Julian date, and ``tdb2`` for the fraction that specifies the
time of day.  Nearly all ``jplephem`` routines accept this optional
``tdb2`` argument if you wish to provide it, thanks to the work of
Marten van Kerkwijk!


Legacy Ephemeris Packages
-------------------------

Back before I learned about SPICE and SPK files, I had run across the
text-file formatted JPL ephemerides at:

ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/

I laboriously assembled the data in these text files into native NumPy
array files and wrapped them each in a Python package so that users
could install an ephemeris with a simple command::

    pip install de421

If you want to use one of these pip-installable ephemerides, you will be
using a slightly older API, and will lose the benefit of the efficient
memory-mapping that the newer SPK code performs.  With the old API, here
is how you would load DE421 and compute a position, given a barycentric
dynamical time expressed as a Julian date::

    import de421
    from jplephem import Ephemeris

    eph = Ephemeris(de421)
    x, y, z = eph.position('mars', 2444391.5)  # 1980.06.01

For more information about the legacy API, consult the ``jplephem``
entry on PyPI for the final release of the 1.x series:

https://pypi.python.org/pypi/jplephem/1.2

The ephemerides that were made available as Python packages (the
following links explain the differences between them) are:

* `DE405 <http://pypi.python.org/pypi/de405>`_ (May 1997)
  — 54 MB covering years 1600 through 2200
* `DE406 <http://pypi.python.org/pypi/de406>`_ (May 1997)
  — 190 MB covering years -3000 through 3000
* `DE421 <http://pypi.python.org/pypi/de421>`_ (February 2008)
  — 27 MB covering years 1900 through 2050
* `DE422 <http://pypi.python.org/pypi/de422>`_ (September 2009)
  — 531 MB covering years -3000 through 3000
* `DE423 <http://pypi.python.org/pypi/de423>`_ (February 2010)
  — 36 MB covering years 1800 through 2200


Reporting issues
----------------

You can report any issues, bugs, or problems at the GitHub repository:

https://github.com/brandon-rhodes/python-jplephem/


Changelog
---------

**2016 December 19 — Version 2.6**

* Fixed the ability to invoke the module from the command line with
  ``python -m jplephem``, and added a test to keep it fixed.

**2015 November 9 — Version 2.5**

* Move ``fileno()`` call out of the ``DAF`` constructor to support
  fetching at least summary information from ``StringIO`` objects.

**2015 November 1 — Version 2.4**

* Add Windows compatibility by switching ``mmap()`` from using
  ``PAGESIZE`` to ``ALLOCATIONGRANULARITY``.

* Avoid a new NumPy deprecation warning by being careful to use only
  integers in the NumPy ``shape`` tuple.

* Add names "TDB" and "TT" to the names database for DE430.

**2015 August 16 — Version 2.3**

* Added auto-detection and support for old NAIF/DAF kernels like
  ``de405.bsp`` to the main ``DAF`` class itself, instead of requiring
  the awkward use of an entirely different alternative class.

**2015 August 5 — Version 2.2**

* You can now invoke ``jplephem`` from the command line.

* Fixes an exception that was raised for SPK segments with a coefficient
  count of only 2, like the DE421 and DE430 segments that provide the
  offset of Mercury from the Mercury barycenter.

* Supports old NAIF/DAF kernels like ``de405.bsp``.

* The ``SPK()`` constructor is now simpler, taking a ``DAF`` object
  instead of an open file.  This is considered an internal API change —
  the public API is the constructor ``SPK.open()``.

**2015 February 24 — Version 2.1**

* Switched from mapping an entire SPK file into memory at once to
  memory-mapping each segment separately on demand.

**2015 February 8 — Version 2.0**

* Added support for SPICE SPK kernel files downloaded directly from
  NASA, and designated old Python-packaged ephemerides as “legacy.”

**2013 November 26 — Version 1.2**

* Helge Eichhorn fixed the default for the ``position_and_velocity()``
  argument ``tdb2`` so it defaults to zero days instead of 2.0 days.
  Tests were added to prevent any future regression.

**2013 July 10 — Version 1.1**

* Deprecates the old ``compute()`` method in favor of separate
  ``position()`` and ``position_and_velocity()`` methods.

* Supports computing position and velocity in two separate phases by
  saving a “bundle” of coefficients returned by ``compute_bundle()``.

* From Marten van Kerkwijk: a second ``tdb2`` time argument, for users
  who want to build higher precision dates out of two 64-bit floats.

**2013 January 18 — Version 1.0**

* Initial release


References
----------

The Jet Propulsion Laboratory's “Solar System Dynamics” page introduces
the various options for doing solar system position computations:
http://ssd.jpl.nasa.gov/?ephemerides

The plain ASCII format element sets from which the ``jplephem`` Python
ephemeris packages are built, along with documentation, can be found at:
ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/

Equivalent FORTRAN code for using the ephemerides be found at the same
FTP site: ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/

"""
from .ephem import Ephemeris, DateError

__all__ = ['Ephemeris', 'DateError']
