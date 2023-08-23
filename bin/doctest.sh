#!/bin/bash

cd $(dirname "$0")
cd ..
cd ci
export PYTHONPATH=..
python3 -m doctest -o ELLIPSIS /dev/stdin < ../jplephem/__init__.py
