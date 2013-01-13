#!/bin/bash

# Download and unzip the ephemeris package files back into their source
# directories, to avoid having to do a full rebuild when testing (which
# would otherwise require a fetch.sh followed by a build.sh).

set -e

for url in \
    http://jplephem.s3.amazonaws.com/de405-1997.tar.gz \
    http://jplephem.s3.amazonaws.com/de406-1997.tar.gz \
    http://jplephem.s3.amazonaws.com/de421-2008.tar.gz \
    http://jplephem.s3.amazonaws.com/de422-2009.tar.gz \
    http://jplephem.s3.amazonaws.com/de423-2010.tar.gz
do
    filename=$(basename $url)
    name=$(echo $filename | sed 's/-.*//')

    wget -c $url
    mkdir tmp
    cd tmp
    tar xvfz ../$filename
    mv */*/* ../$name/$name
    cd ..
    rm -rf tmp
done

