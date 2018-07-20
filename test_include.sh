#!/bin/bash

# Tests that all hpp-files can be compiled by themselves, to check that they
# do not have missing includes.
for f in *.hpp
do
    echo "Testing ${f}"
    if ! g++ "${f}" -o /dev/null -std=c++11 > /dev/null 2> /dev/null
    then
        echo "Failure in ${f}"
        exit 1
    fi
done

echo "Success"
exit 0
