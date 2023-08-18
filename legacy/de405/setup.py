from setuptools import setup

name = 'de405'
module = __import__(name)
description, long_description = module.__doc__.split('\n', 1)

setup(name = name,
      version = '1997.1',
      description = description,
      long_description = long_description,
      license = 'MIT',
      author = 'Brandon Rhodes',
      author_email = 'brandon@rhodesmill.org',
      url = 'https://github.com/brandon-rhodes/python-jplephem/',
      classifiers = [
        'Development Status :: 5 - Production/Stable',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 3',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
      packages = [name],
      package_data = {name: ['*.npy',]},
      )
