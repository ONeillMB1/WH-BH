#!/usr/bin/env python

import sys

def usage():
    if len(sys.argv) < 4:
        print("Usage: <program> <excludeBed> <inFile> <outfile>")
        sys.exit(0)

excludeBedFile, inFileName, outFileName = sys.argv[1:]

def make
