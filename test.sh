#!/bin/bash

set -e

PYTHON=${PYTHON:-python}

cd "$(dirname "${BASH_SOURCE[0]}")"
cd ci
PYTHONPATH=$(dirname $PWD) $PYTHON -m unittest test.py "$@"
