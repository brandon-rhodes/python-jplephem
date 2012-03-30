import sys
from distutils.core import setup

packages = ['jplephem', 'de405', 'de406', 'de422', 'de423']
if len(sys.argv) < 2 or sys.argv[1] not in packages:
    print 'usage: setup.py %s ...' % '|'.join(packages)
    sys.exit(2)

name = sys.argv[1]
del sys.argv[1]

if name.startswith('d'):
    module = __import__(name)
    description, long_description = name.__doc__.split('\n', 1)
    setup(name = name,
          version = '1.0',
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
