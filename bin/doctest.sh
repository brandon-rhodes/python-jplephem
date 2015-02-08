#!/bin/bash

cd $(dirname "$0")
cd ..
python -m doctest /dev/stdin < jplephem/__init__.py
