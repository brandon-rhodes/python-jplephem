from distutils.core import setup

versions = {'jplephem': '0.1', 'de405': '1997', 'de406': '1997',
            'de421': '2008', 'de422': '2009', 'de423': '2010'}

# Fake the presence of numpy so __import__() always succeeds.
import sys
sys.modules['numpy'] = sys.modules['sys']

name = 'jplephem'
module = __import__(name)
description, long_description = module.__doc__.split('\n', 1)

setup(name = name,
      version = versions[name],
      description = description,
      long_description = long_description,
      license = 'MIT',
      author = 'Brandon Rhodes',
      author_email = 'brandon@rhodesmill.org',
      url = ('http://jplephem.s3.amazonaws.com/packages.html'
             if name.startswith('d') else None),
      classifiers = [
        'Development Status :: 5 - Production/Stable' if name.startswith('d')
        else 'Development Status :: 4 - Beta',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 2',
        'Programming Language :: Python :: 2.7',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.2',
        'Topic :: Scientific/Engineering :: Astronomy',
        ],
      packages = [name],
      package_data = {name: ['*.npy',]},
      install_requires = ['numpy'],
      )
