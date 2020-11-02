from distutils.core import setup, Extension
from Cython.Build import cythonize


setup(ext_modules = cythonize(Extension(
           "core",                                # the extension name
           sources=["core.pyx", "hapcolor.cpp"], # the Cython source and
                                                  # additional C++ source files
           language="c++",                        # generate and compile C++ code
      )))