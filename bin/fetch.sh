#!/bin/bash

cd $(dirname "$0")/..
pwd
for n in 405 406 421 422 423
do
    wget -rc ftp://ssd.jpl.nasa.gov/pub/eph/planets/ascii/de$n/
done
