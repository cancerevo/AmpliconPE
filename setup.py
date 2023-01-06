#!/usr/bin/env python3
from setuptools import setup, Extension

# from distutils.core import setup, Extension
import os

# Command for troubleshooting the SSW library
# (setuptools doesn't provide much feedback when it fails)
#   gcc -Wall -pipe -O2 -fPIC -shared -rdynamic -o libssw.so ssw.c

ssw_dir = "Complete-Striped-Smith-Waterman-Library/src/"

if __name__ == "__main__":
    setup(
        ext_modules=[
            Extension(
                "AmpliconPE.libssw",
                include_dirs=[ssw_dir],
                sources=[os.path.join(ssw_dir, "ssw.c")],
            )
        ]
    )
