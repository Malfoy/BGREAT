#!/usr/bin/env python
import sys

sequence = ""
with open(sys.argv[1]) as infile:  # do not load the whole file, one line at time
    for line in infile:
        string = line.rstrip()
        if ">" not in string:
            sequence += string
        else:
            if sequence != "":
                print sequence
            print string
            sequence = ""
    print sequence



# lines = [line.rstrip() for line in open(sys.argv[1], "r")]
# i = 0
# sequence = ""
# while i < len(lines):
#     string = lines[i]
#     if ">" not in string:
#         sequence += string
#     else:
#         if sequence != "":
#             print sequence
#         print string
#         sequence = ""
#     i += 1
# print sequence
