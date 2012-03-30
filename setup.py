import os
import sys
from distutils.core import setup

versions = {'jplephem': '0.1', 'de405': '1997', 'de406': '1997',
            'de421': '2008', 'de422': '2009', 'de423': '2010'}

# If this setup.py is being run from inside of an "sdist" that has been
# distributed on the Python Package Index, then it should be ready to
# simply install the package that it was distributed with.

if os.path.exists('PKG-INFO'):
    with open('PKG-INFO') as f:
        words = f.read().split()
    name = words[words.index('Name:') + 1]

# Otherwise, a developer must be running this "setup.py" from inside of
# a checkout of the "jplephem" repository.  In that case, they need to
# specify the package they mean using the command line.

else:
    if len(sys.argv) < 2 or sys.argv[1] not in versions:
        print 'usage: setup.py %s ...' % '|'.join(versions.keys())
        sys.exit(2)
    name = sys.argv.pop(1)

# The ephemeris packages are built very simply.

module = __import__(name)
description, long_description = module.__doc__.split('\n', 1)

setup(name = name,
      version = versions[name],
      description = description,
      long_description = long_description,
      license = 'MIT',
      author = 'Brandon Rhodes',
      author_email = 'brandon@rhodesmill.org',
      classifiers = [
        'Development Status :: 5 - Production/Stable' if name.startswith('d')
        else 'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
      packages = [name],
      package_data = {name: ['*.npy',]},
      )
