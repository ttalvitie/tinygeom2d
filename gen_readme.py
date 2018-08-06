#!/usr/bin/env python3

print("""\
# tinygeom2d
Tiny 2D geometry library

## API reference""")

headers = [
    ("geometry.hpp", "Geometry primitives"),
    ("domain.hpp", "Polygonal domains"),
    ("intersection.hpp", "Intersection detection"),
    ("visibility.hpp", "Visibility computations"),
]

for (header, title) in headers:
    print("### {}: {}".format(header, title))
