from distutils.core import setup

# Fake the presence of numpy so import can succeed.
import sys
sys.modules['numpy'] = sys.modules['sys']

import jplephem
description, long_description = jplephem.__doc__.split('\n', 1)

setup(name = 'jplephem',
      version = jplephem.__version__,
      description = description,
      long_description = long_description,
      license = 'MIT',
      author = 'Brandon Rhodes',
      author_email = 'brandon@rhodesmill.org',
      url = 'https://github.com/brandon-rhodes/python-jplephem/',
      classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
      packages = ['jplephem'],
      install_requires = ['numpy'],
      )
