#!/usr/bin/env python3

import re
import sys

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
        match = re.search("namespace\s+[a-zA-Z_]*_detail", code)
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
    
    code = re.sub(r'(^\s*//.*$\n)+\s*$', r'', code, flags=re.MULTILINE)
    code = re.sub(r'(^\s*$\n)+', r'\n', code, flags=re.MULTILINE)
    
    print("```c++")
    print(code)
    print("```")
