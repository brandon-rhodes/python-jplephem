
Welcome to the repository for the `jplephem` Python library!

The package is a Python implementation of the math that standard JPL
ephemerides use to predict raw (x,y,z) planetary positions.  It is one
of the foundations of the Skyfield astronomy library for Python:

http://rhodesmill.org/skyfield/

But you can also use `jplephem` standalone to generate raw vectors.  If
that is your use case, then simply head over to its documentation and
download link on the Python Package Index:

https://pypi.python.org/pypi/jplephem

If you want to install it with `conda`, there is a recipe at:

https://github.com/conda-forge/jplephem-feedstock

This repository is where `jplephem` is maintained.  You will find its
source code beneath the `jplephem` directory that sits alongside the
`setup.py` file.  You can run its tests with:

    wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de405.bsp
    wget https://naif.jpl.nasa.gov/pub/naif/generic_kernels/spk/planets/a_old_versions/de421.bsp
    pip install de421
    python -m unittest discover jplephem

Enjoy!
