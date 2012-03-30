import sys
from distutils.core import setup

versions = {'jplephem': '0.1', 'de405': '1997', 'de406': '1997',
            'de421': '2008', 'de422': '2009', 'de423': '2010'}
if len(sys.argv) < 2 or sys.argv[1] not in versions:
    print 'usage: setup.py %s ...' % '|'.join(versions.keys())
    sys.exit(2)

name = sys.argv.pop(1)

if name.startswith('d'):
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
            'Development Status :: 5 - Production/Stable',
            'Intended Audience :: Science/Research',
            'License :: OSI Approved :: MIT License',
            'Topic :: Scientific/Engineering :: Astronomy',
            ],
          packages = [name],
          package_data = {name: ['*.npy',]},
          )
