#!/usr/bin/env python3

import re
import sys
import subprocess
import tempfile

subprocess.check_call(["make"], cwd="examples", stdout=subprocess.DEVNULL)

print("""\
# tinygeom2d
tinygeom2d is a tiny C++ library for 2D geometry in polygonal domains. Features:

* Shortest paths
* Visibility polygons
* Visibility graphs
* Header-only library
* Exact geometry primitives for points with 63-bit integer coordinates
* Simplicity and robustness achieved through exact computations and symbolic perturbations

## Why integers?
Even though floating point arithmetic is more flexible than integer arithmetic, it requires greater care in the implementation of geometric algorithms to account for the imprecision of the computations. For this reason, the library uses integer arithmetic in most places. To use the library with floating point data, the input has to be transformed into integer data first (see function `normalizationFactor` in `geometry.hpp`). This should not greatly reduce the precision of the data, because the precision of coordinates in tinygeom2d is greater than the 53-bit precision of double-precision floating point numbers (points in polygonal domains are typically distributed uniformly, unlike floating point numbers that are denser close to zero). Alternative libraries [VisiLibity1](https://karlobermeyer.github.io/VisiLibity1/) and [CGAL](https://www.cgal.org/) can use floating point data directly.

Exact arithmetic also makes it possible to use symbolic perturbations that remove the need for handling degenerate cases (such as three points being collinear) by symbolically moving each point an infinitesimally small distance from its original location. For more information about the symbolic perturbation used by the library, see the comments in `geometry.hpp`.

## Usage
tinygeom2d is a header-only library, so you can start using it simply by copying the tinygeom2d directory to your project directory (or anywhere in the include path). C++11 support is required from the compiler. The library is licensed with the permissive MIT license, so you can use it however you like provided that you retain the copyright notice and license terms in all copies. See the examples and API reference below for instructions on how to use the library.

# Examples""")

examples = [
    ("domain", "Domains and points"),
    ("visibility", "Visibility"),
    ("shortestpath", "Shortest paths"),
]

for (example, title) in examples:
    print("## {}".format(title))
    with open("examples/{}.cpp".format(example)) as fp:
        code = ""
        noreadme = False
        for line in fp:
            if "BEGIN NOREADME" in line:
                noreadme = True
            elif "END NOREADME" in line:
                noreadme = False
            else:
                if not noreadme:
                    code += line
    code = code.strip()
    
    output = subprocess.check_output(["./" + example], cwd="examples")
    output = output.decode("UTF-8")
    output = output.strip()
    
    # Run the code with NOREADME parts removed and compare output to verify
    # that it still works
    with tempfile.TemporaryDirectory() as tmpdir:
        with open(tmpdir + "/code.cpp", "w") as fp:
            fp.write(code)
        subprocess.check_call(["g++", tmpdir + "/code.cpp", "-o", tmpdir + "/code", "-std=c++11", "-I."])
        output_cmp = subprocess.check_output([tmpdir + "/code"])
        output_cmp = output_cmp.decode("UTF-8")
        output_cmp = output_cmp.strip()
        assert output_cmp == output
    
    print("```c++")
    print(code)
    print("```")
    print("Output:")
    print("```")
    print(output)
    print("```")
    print("![](examples/{}.svg)".format(example))
    print()

print("# API reference")

headers = [
    ("geometry.hpp", "Geometry primitives"),
    ("intersection.hpp", "Intersection detection"),
    ("domain.hpp", "Polygonal domains"),
    ("visibility.hpp", "Visibility computations"),
    ("shortestpath.hpp", "Shortest path computations"),
    ("int64.hpp", "Portable 64-bit integer multiplication comparison"),
]

for (header, title) in headers:
    print("## {}: {}".format(header, title))
    with open("tinygeom2d/" + header) as fp:
        code = fp.read()
    
    # Ugly regex/manual string handling hack that removes unnecessary things from the code
    code = re.sub(r'^\s*#.*$', r'', code, flags=re.MULTILINE)
    
    code = code.strip()
    pos = -1
    while True:
        pos = code.find("{", pos + 1)
        if pos == -1:
            break
        pred = code[code.rfind("\n", 0, pos)+1:pos]
        
        if not ("struct" in pred or "class" in pred or "namespace" in pred):
            level = 0
            i = pos + 1
            while True:
                if code[i] == "{":
                    level += 1
                elif code[i] == "}":
                    level -= 1
                    if level == -1:
                        break
                i += 1
            code = code[:pos] + ";" + code[i + 1:]
    
    code = re.sub(r'\s*;', r';', code)
    code = re.sub(r'\)\s*:[^;]*;', r');', code)
    code = re.sub(r'^(\s*)inline ', r'\1', code, flags=re.MULTILINE)
    code = re.sub(r'^\s*$', r'', code, flags=re.MULTILINE)
    
    pos = -1
    while True:
        pos = code.find("private:", pos + 1)
        if pos == -1:
            break
        
        i = pos
        level = 0
        while True:
            if code[i] == "{":
                level += 1
            elif code[i] == "}":
                level -= 1
                if level == -1:
                    break
            i += 1
        a = pos
        if code[a-2:a] == "\n\n":
            a -= 1
        code = code[:a] + code[i:]
    
    while True:
        match = re.search("namespace\s+[a-zA-Z0-9_]*_detail", code)
        if match == None:
            break
        pos = match.start()
        level = 0
        i = pos
        while True:
            if code[i] == "{":
                level += 1
            elif code[i] == "}":
                level -= 1
                if level == 0:
                    break
            i += 1
        code = code[:pos] + code[i+1:]
    
    code = re.sub(r'^\s*$\n^int cmpMul.*$', r'', code, flags=re.MULTILINE)
    code = re.sub(r'(^\s*//.*$\n)+\s*$', r'', code, flags=re.MULTILINE)
    code = re.sub(r'(^\s*$\n)+', r'\n', code, flags=re.MULTILINE)
    
    print("```c++")
    print(code)
    print("```")
    print()

print("""\
# Credits
The library was written by Topi Talvitie in 2018, based on earlier code of the following two visualizations written as a research assistant in the computational geometry research group in University of Helsinki supervised by Valentin Polishchuk:

* J. Hershberger, V. Polishchuk, B. Speckmann, T. Talvitie: Geometric kth Shortest Paths: the Applet, SoCG 2014 videos [[paper and visualization](http://www.computational-geometry.org/SoCG-videos/socg14video/ksp/index.html)]
* T. Talvitie: Visualizing Quickest Visibility Maps, SoCG 2015 videos [[paper](http://drops.dagstuhl.de/opus/volltexte/2015/5090/)] [[visualization](https://www.cs.helsinki.fi/group/compgeom/qvm/)]

Support for kth shortest paths and quickest visibility paths may be added to the library in the future.""")
