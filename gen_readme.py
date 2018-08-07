#!/usr/bin/env python3

import re
import sys
import subprocess

subprocess.check_call(["make"], cwd="examples", stdout=subprocess.DEVNULL)

print("""\
# tinygeom2d
Tiny 2D geometry library
# Examples""")

examples = [
    ("domain", "Domains and points"),
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
    
    print("```c++")
    print(code)
    print("```")
    print("Output:")
    print("```")
    print(output)
    print("```")
    print("![](examples/{}.svg)".format(example))

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
