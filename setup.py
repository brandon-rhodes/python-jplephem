from setuptools import setup

import os, sys
sys.path.append(os.path.abspath(os.path.dirname(__file__)))

import jplephem
description, long_description = jplephem.__doc__.split('\n', 1)

setup(
    description = description,
    long_description = long_description,
)
