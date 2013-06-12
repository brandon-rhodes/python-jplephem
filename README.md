
Welcome to the repository for the `jplephem` Python library!

The package is a Python implementation of the math that standard JPL
ephemerides use to predict raw (x,y,z) planetary positions.  If you are
simply interested in using it, please head over to its documentation and
download link on the Python Package Index:

https://pypi.python.org/pypi/jplephem

This repository is where the development version of `jplephem` is
maintained.  You will find its source code beneath the `jplephem`
directory, alongside its `setup.py` file.

This repository also contains the tools for packaging JPL ephemerides so
they are easy to install and import as Python packages.  For example, to
build your own copies of the ephemerides that are listed down at the
bottom of the `jplephem` documentation page, you could run the following
commands (skip the first two of these commands if you have already
cloned this repository and are sitting in its directory):

    git clone https://github.com/brandon-rhodes/python-jplephem.git
    cd jplephem
    bin/fetch.sh
    bin/build.py

The result should be fully-compiled ephemerides ready to be installed!

If, instead of building the standard set of publicly supported
ephemerides, you are interested in trying some of the more obscure
offerings from the JPL, then select one of the directory names listed at
the following URL:

ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/

Some of the JPL ephemeris directories include a README file with more
details, and the properties of the main ephemerides are documented in
the top-level README at:

ftp://ssd.jpl.nasa.gov/pub/eph/planets/README.txt

Once you have selected an ephemeris and know its JPL directory name, go
edit your copy of `bin/fetch.sh` and replace the standard ephemeris
directory names with the alternative ephemeris that you would like to
try out, and then re-run both `bin/fetch.sh` and `bin/build.py` as shown
above.  The result should be a new directory complete with:

* A `setup.py` file.
* An `__init__.py` file one directory lower.
* Several `.npy` files with actual ephemeris data.

You can then install the ephemeris locally with:

    python setup.py install

Or, you can ask for it to be packaged as a `.tar.gz` file with:

    python setup.py sdist

Either way, your new ephemeris should work like the existing ephemerides
that are already listed on the Python Package Index for use with
`jplephem` — enjoy!
