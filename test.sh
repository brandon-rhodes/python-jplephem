#!/bin/bash

set -e

cd "$(dirname "${BASH_SOURCE[0]}")"
cd ci
unit2 jplephem.test "$@"
