"""JPL Planetary and Lunar Ephemeris DE422 for the jplephem package.

This is the most recent long-period ephemeris published by the Jet
Propulsion Laboratory.  While requiring more than half a gigabyte of
space, it achieves quite high accuracy.

:Name: DE422 (September 2009)
:Years: -3000 through 3000
:Planets: Yes
:Sun/Moon: Yes
:Nutations: Yes
:Librations: Yes
:Report: `Jones, Fomalont, Dhawan, Romney, Folkner, Lanyi, Border, Jacobson (2010) [PDF] <http://arxiv.org/pdf/1012.0264>`_
:Size: 531 MB

This ephemeris incorporates ranging data from recent spacecraft
missions, including the Cassini mission to Saturn, which results in an
accuracy for Saturn of a few milliarcseconds for positions over the past
decade.  If this ephemeris is too large for your application, take a
look at `DE421 <http://pypi.python.org/pypi/de421>`_ and `DE423
<http://pypi.python.org/pypi/de423>`_ as alternatives.

To compute using this ephemeris in Python, see the `jplephem
<http://pypi.python.org/pypi/jplephem>`_ package.

"""
