#!/bin/bash

for file in \
    de405-1997.tar.gz \
    de406-1997.tar.gz \
    de421-2008.tar.gz \
    de422-2009.tar.gz \
    de423-2010.tar.gz \

do
    s3cmd put -P dist/$file s3://jplephem/$f
done

s3cmd put -P packages.html s3://jplephem/packages.html
