#!/bin/bash

# Tests that two compilation units that include the library can be linked
# together, meaning that there are no non-inline functions or other symbols
# such that linking multiple instances of them is an error.

report_error() {
    echo "Test script internal error on line $1"
}

trap 'report_error $LINENO' ERR

set -e

# Create temporary directory
TMPDIR=`mktemp -d`

# Add library to its own directory
mkdir "${TMPDIR}/tinygeom2d"
for f in tinygeom2d/*.hpp
do
    ln -s "$(pwd)/${f}" "${TMPDIR}/${f}"
done

# Create two source files that both include all headers, and one of them
# contains main
for f in tinygeom2d/*.hpp
do
    echo "#include \"${f}\"" >> "${TMPDIR}/a.cpp"
    echo "#include \"${f}\"" >> "${TMPDIR}/b.cpp"
done
echo "int main(int, char**) { }" >> "${TMPDIR}/a.cpp"

# Try to compile the source files into a binary
cd "${TMPDIR}"
if g++ a.cpp b.cpp -o /dev/null -std=c++11 -w
then
    echo "Success"
    exit 0
else
    echo "Failure"
    exit 1
fi
