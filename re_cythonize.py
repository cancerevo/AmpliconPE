#!/usr/bin/env python3

from distutils.core import Extension
from Cython.Build import cythonize
import configparser

# import numpy

config = configparser.RawConfigParser()
config.read("setup.cfg")
name = config["metadata"].get("name")

extensions = [Extension(name, [f"{name}/{name}.pyx"], extra_compile_args=["-O2"])]
cythonize(extensions, compiler_directives={"language_level": "3"})
