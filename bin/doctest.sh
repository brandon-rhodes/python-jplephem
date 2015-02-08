#!/bin/bash

cd $(dirname "$0")
cd ..
python -m doctest -o ELLIPSIS /dev/stdin < jplephem/__init__.py
