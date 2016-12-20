from distutils.core import setup

# Fake the presence of numpy so import can succeed.
import sys
sys.modules['numpy'] = sys.modules['sys']

import jplephem
description, long_description = jplephem.__doc__.split('\n', 1)

setup(name = 'jplephem',
      version = '2.6',
      description = description,
      long_description = long_description,
      license = 'MIT',
      author = 'Brandon Rhodes',
      author_email = 'brandon@rhodesmill.org',
      classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.6',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
      packages = ['jplephem'],
      requires = ['numpy'],
      )
