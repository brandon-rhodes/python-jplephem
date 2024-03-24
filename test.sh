#!/bin/bash

set -e

PYTHON=${PYTHON:-python}

cd "$(dirname "${BASH_SOURCE[0]}")"
projectdir=$PWD
cd ci
PYTHONPATH=$projectdir $PYTHON -m unittest test "$@"
