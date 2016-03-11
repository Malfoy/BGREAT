#!/usr/bin/env python
import sys

with open(sys.argv[1]) as infile:
    for line in infile:
        string = line.rstrip()
        if ">" in string:
            header=string
        elif not "N" in string:
            print header
            print string
