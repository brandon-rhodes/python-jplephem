# -*- encoding: utf-8 -*-

"""Use a JPL ephemeris to predict planet positions.

Cite as: `Astrophysics Source Code Library, record ascl:1112.014
<https://ascl.net/1112.014>`_

This package can load and use a Jet Propulsion Laboratory (JPL)
ephemeris for predicting the position and velocity of a planet or other
Solar System body.  It currently supports binary SPK files (extension
``.bsp``) like `those distributed by the Jet Propulsion Laboratory
<https://ssd.jpl.nasa.gov/ftp/eph/planets/bsp/>`_ that are:

* **Type 2** — positions stored as Chebyshev polynomials, with velocity
  derived by computing their derivative.

* **Type 3** — positions and velocities both stored explicitly as
  Chebyshev polynomials.

* **Type 9** — a series of discrete positions and velocities, with
  separate timestamps that do not need to be equally spaced.  Currently
  there is only support for linear interpolation: for Type 9 ephemerides
  of polynomial degree 1, not of any higher degrees.

Note that even if an ephemeris isn’t one of the above types, you can
still use ``jplephem`` to read its text comment and list the segments
inside, using the subcommands ``comment`` and ``daf`` described below.

Installation
------------

The only third-party package that ``jplephem`` depends on is `NumPy
<http://www.numpy.org/>`_, which ``pip`` will automatically attempt to
install alongside ``pyephem`` when you run::

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

If you have downloaded a ``.bsp`` file, you can run ``jplephem`` from
the command line to display the data inside of it::

    python -m jplephem comment de421.bsp
    python -m jplephem daf de421.bsp
    python -m jplephem spk de421.bsp

You can also take a large ephemeris and produce a smaller excerpt by
limiting the range of dates that it covers::

    python -m jplephem excerpt 2018/1/1 2018/4/1 de421.bsp excerpt421.bsp

You will get an error if your starting year is negative, because Unix
commands expect a list of options when they see a dash.  The fix is to
provide a special argument ``--`` which says “I’m done passing options,
even if the next argument stars with a dash”::

    python -m jplephem excerpt -- -800/1/1 800/1/1 de422.bsp excerpt422.bsp

You can also filter by the integer codes for the targets you need.
Unrecognized targets will not raise an error, to let you apply a master
list of targets to a whole series of SPK files that might or might not
each have all of the targets::

    python -m jplephem excerpt --targets 1,2,3 2018/1/1 2018/4/1 de421.bsp excerpt421.bsp

If the input ephemeris is a URL, then ``jplephem`` will try to save
bandwidth by fetching only the blocks of the remote file that are
necessary to cover the dates you have specified.  For example, the
Jupiter satellite ephemeris ``jup310.bsp`` is famously large, weighing
in a nearly a gigabyte.  But if all you need are Jupiter's satellites
for a few months, you can download considerably less data::

    $ python -m jplephem excerpt 2018/1/1 2018/4/1 \\
        https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/satellites/jup310.bsp \\
        excerpt.bsp
    $ ls -lh excerpt.bsp
    -rw-r----- 1 brandon brandon 1.2M Feb 11 13:36 excerpt.bsp

In this case only about one-thousandth of the ephemeris's data needed to
be downloaded.

Getting Started With DE421
--------------------------

The DE421 ephemeris is a useful starting point.  It weighs in at 17 MB,
but provides predictions over the years 1900–2050:

https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp

After the kernel has downloaded, you can use ``jplephem`` to load this
SPK file and learn about the segments it offers:

>>> from jplephem.spk import SPK
>>> kernel = SPK.open('de421.bsp')
>>> print(kernel)
File type DAF/SPK and format LTL-IEEE with 15 segments:
1899-07-29..2053-10-09  Type 2  Solar System Barycenter (0) -> Mercury Barycenter (1)
1899-07-29..2053-10-09  Type 2  Solar System Barycenter (0) -> Venus Barycenter (2)
1899-07-29..2053-10-09  Type 2  Solar System Barycenter (0) -> Earth Barycenter (3)
1899-07-29..2053-10-09  Type 2  Solar System Barycenter (0) -> Mars Barycenter (4)
1899-07-29..2053-10-09  Type 2  Solar System Barycenter (0) -> Jupiter Barycenter (5)
1899-07-29..2053-10-09  Type 2  Solar System Barycenter (0) -> Saturn Barycenter (6)
1899-07-29..2053-10-09  Type 2  Solar System Barycenter (0) -> Uranus Barycenter (7)
1899-07-29..2053-10-09  Type 2  Solar System Barycenter (0) -> Neptune Barycenter (8)
1899-07-29..2053-10-09  Type 2  Solar System Barycenter (0) -> Pluto Barycenter (9)
1899-07-29..2053-10-09  Type 2  Solar System Barycenter (0) -> Sun (10)
1899-07-29..2053-10-09  Type 2  Earth Barycenter (3) -> Moon (301)
1899-07-29..2053-10-09  Type 2  Earth Barycenter (3) -> Earth (399)
1899-07-29..2053-10-09  Type 2  Mercury Barycenter (1) -> Mercury (199)
1899-07-29..2053-10-09  Type 2  Venus Barycenter (2) -> Venus (299)
1899-07-29..2053-10-09  Type 2  Mars Barycenter (4) -> Mars (499)

Since the next few examples involve vector output, let’s tell NumPy to
make vector output attractive.

>>> import numpy as np
>>> np.set_printoptions(precision=3)

Each segment of the file lets you predict the position of one body with
respect to another for a given Julian date.  A small routine is provided
to convert Gregorian calendar dates to Julian dates:

>>> from jplephem.calendar import compute_julian_date
>>> compute_julian_date(2015, 2, 8)
2457061.5

Here is how to compute the coordinates of Mars (target 4) relative to
the Solar System barycenter (target 0) at midnight 2015 February 8 TDB
(Barycentric Dynamical Time), using the Julian date we just computed:

>>> position = kernel[0,4].compute(2457061.5)
>>> print(position)
[2.057e+08 4.251e+07 1.394e+07]

By contrast, it takes three steps to learn the position of Mars with
respect to the Earth: from Mars to the Solar System barycenter, to the
Earth-Moon barycenter (3), and finally to Earth itself (399).

>>> position = kernel[0,4].compute(2457061.5)
>>> position -= kernel[0,3].compute(2457061.5)
>>> position -= kernel[3,399].compute(2457061.5)
>>> print(position)
[ 3.161e+08 -4.679e+07 -2.476e+07]

You can see that the output of this ephemeris DE421 is in kilometers.
If you use another ephemeris, check its documentation to be sure of the
units that it employs.

If you supply the date as a NumPy array, then each component that is
returned will itself be a vector as long as your date:

>>> jd = np.array([2457061.5, 2457062.5, 2457063.5, 2457064.5])
>>> position = kernel[0,4].compute(jd)
>>> print(position)
[[2.057e+08 2.053e+08 2.049e+08 2.045e+08]
 [4.251e+07 4.453e+07 4.654e+07 4.855e+07]
 [1.394e+07 1.487e+07 1.581e+07 1.674e+07]]

Some ephemerides include velocity inline by returning a 6-vector instead
of a 3-vector.  For an ephemeris that does not, you can ask for the
Chebyshev polynomial to be differentiated to produce a velocity, which
is delivered as a second return value:

>>> position, velocity = kernel[0,4].compute_and_differentiate(2457061.5)
>>> print(position)
[2.057e+08 4.251e+07 1.394e+07]
>>> print(velocity)
[-363896.059 2019662.996  936169.773]

The velocity will by default be distance traveled per day, in whatever
units for distance the ephemeris happens to use.  To get a velocity per
second, simply divide by the number of seconds in a day:

>>> velocity_per_second = velocity / 86400.0
>>> print(velocity_per_second)
[-4.212 23.376 10.835]

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
  1899-07-29..2053-10-09  Type 2  Earth Barycenter (3) -> Earth (399)
    frame=1 source=DE-0421LE-0421

* Each ``Segment`` loaded from the kernel has a number of attributes
  that are loaded from the SPK file:

  >>> from jplephem.spk import BaseSegment
  >>> help(BaseSegment)
  Help on class BaseSegment in module jplephem.spk:
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

* If you want to access the raw coefficients, use the segment
  ``load_array()`` method.  It returns two floats and a NumPy array:

  >>> initial_epoch, interval_length, coefficients = segment.load_array()
  >>> print(coefficients.shape)
  (3, 14080, 13)

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

Support for Binary PCKs
-----------------------

You can also load and produce rotation matrices from a binary PCK file.
Its segments are available through the ``segments`` attributes of the
returned object.

>>> from jplephem.pck import PCK
>>> p = PCK.open('moon_pa_de421_1900-2050.bpc')
>>> p.segments[0].body
31006
>>> p.segments[0].frame
1
>>> p.segments[0].data_type
2

Given a solary system barycenter Julian date, the segment will return
the three angles necessary to build a rotation matrix: right ascension
of the pole, declination of the pole, and cumulative rotation of the
body’s axis.  Typically these will all be in radians.

>>> tdb = 2454540.34103
>>> print(p.segments[0].compute(tdb, 0.0, False))
[3.928e-02 3.878e-01 3.253e+03]

You can ask for velocity as well.

>>> r, v = p.segments[0].compute(tdb, 0.0, True)
>>> print(r)
[3.928e-02 3.878e-01 3.253e+03]
>>> print(v)
[6.707e-09 4.838e-10 2.655e-06]

Closing an ephemeris
--------------------

To release all open files and memory maps associated with an ephemeris,
call its ``close()`` method.

>>> kernel.close()
>>> p.close()

Reporting issues
----------------

You can report any issues, bugs, or problems at the GitHub repository:

https://github.com/brandon-rhodes/python-jplephem/

Changelog
---------

**2024 April 24 — Version 2.22**

* When printed, segments now print their start and end dates using the
  Gregorian calendar instead of printing raw Julian dates.

* A small ``compute_julian_date`` routine is now provided for converting
  calendar dates into Julian dates.

* Fixed the text of the ``ValueError`` that is raised when the PCK
  segment ``compute()`` method is given an out-of-range date; it was
  reporting incorrectly large numbers for the Julian date range, because
  a PCK counts time using seconds before or after J2000, not years.

**2023 December 1 — Version 2.21**

* Tweaked an import to avoid a fatal exception under Python 2, in case
  anyone is still using it.

**2023 November 13 — Version 2.20**

* Each segment is now protected by a lock, in case two threads
  simultaneously trigger the code that performs the initial load of the
  segment’s data; the symptom was a rare exception ``ValueError: cannot
  reshape array``.

**2023 September 6 — Version 2.19**

* Fixed a bug in the ``excerpt`` command that was causing it to truncate
  its output when the input ephemeris had more than about two dozen
  segments.  The command’s output should now include all matching
  segments from even a very large ephemeris.

* Fixed the ``excerpt`` command so the calendar dates specified on the
  command line produce Julian dates ending with the fraction ``.5``,
  which makes excerpt endpoints more exact.

**2022 September 28 — Version 2.18**

* Added support for big-endian processors, and created a GitHub Actions
  CI build that includes both a big- and a little-endian architecture.

**2021 December 31 — Version 2.17**

* Fixed an ``AttributeError`` in the ``excerpt`` command.

**2021 July 3 — Version 2.16**

* Fixed a ``ValueError`` raised in the ``excerpt`` command when an
  ephemeris segment needs to be entirely skipped because it has no
  overlap with the user-specified range of dates.

* Added a ``__version__`` constant to the package’s top level.

**2020 September 2 — Version 2.15**

* The ``excerpt`` subcommand now accepts a ``--targets`` option to save
  space by copying only matching segments into the output SPK file.

* The Julian day fraction ``tdb2`` is handled even more carefully than
  before, providing a smoother delta between successive positions when
  the difference between successive times is down around 0.1 µs.

**2020 March 26 — Version 2.14**

* Fall back to plain file I/O on platforms that support ``fileno()`` but
  that don’t support ``mmap()``, like the `Pyodide platform
  <https://github.com/iodide-project/pyodide>`_.

**2020 February 22 — Version 2.13**

* The exception raised when a segment is given a Julian date outside the
  segment’s date range is now an instance of the ``ValueError`` subclass
  ``OutOfRangeError`` that reminds the caller of the range of dates
  supported by the SPK segment, and carries an array attribute
  indicating which input dates were at fault.

**2019 December 13 — Version 2.12**

* Replaced use of NumPy ``flip()`` with a reverse slice ``[::-1]`` after
  discovering the function was a recent addition that some user installs
  of NumPy do not support.

**2019 December 13 — Version 2.11**

* Reverse the order in which Chebyshev polynomials are computed to
  slightly increase speed, to simplify the code, and in one case
  (comparing PCK output to NASA) to gain a partial digit of extra
  precision.

**2019 December 11 — Version 2.10**

* Document and release support for ``.bcp`` binary PCK kernel files
  through the new ``jplephem.pck`` module.

**2019 January 3 — Version 2.9**

* Added the ``load_array()`` method to the segment class.

**2018 July 22 — Version 2.8**

* Switched to a making a single memory map of the entire file, to avoid
  running out of file descriptors when users load an ephemeris with
  hundreds of segments.

**2018 February 11 — Version 2.7**

* Expanded the command line tool, most notably with the ability to fetch
  over HTTP only those sections of a large ephemeris that cover a
  specific range of dates, producing a smaller ``.bsp`` file.

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

Equivalent FORTRAN code for using the ephemerides be found at the same
FTP site: ftp://ssd.jpl.nasa.gov/pub/eph/planets/fortran/

"""
from .ephem import Ephemeris, DateError
__version__ = '2.22'

__all__ = ['Ephemeris', 'DateError', '__version__']
