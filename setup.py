import sys
import os
import shutil
import subprocess

from distutils.core import setup, Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize

import numpy

# clean previous build
for root, dirs, files in os.walk("./src/python/", topdown=False):
    for name in dirs:
        if (name == "build"):
            shutil.rmtree(name)


include_dirs = [
    numpy.get_include(),
    "./include",
    "./include/c",
    ]

extra_link_args=[
    "-L./lib/c",
]

fftw_path = os.getenv("FFTW", None)
if fftw_path:
    extra_link_args.append("-L{}/lib".format(fftw_path))
    
if not os.path.exists("./lib/c/libssht.a"):
    print("Making...")
    subprocess.run(["make"])
    print(os.listdir("./lib/c/"))

          
if not os.path.exists("./lib/c/libssht.a"):
    sys.exit()


setup(
    classifiers=['Programming Language :: Python :: 2.7'],
    name = "pyssht",
    version = "2.0",
    cmdclass={'build_ext': build_ext},
    ext_modules=cythonize([Extension(
        "pyssht",
        package_dir=['src'],
        sources=["src/python/pyssht.pyx"],
        include_dirs=include_dirs,
        libraries=["ssht", "fftw3"],
        extra_link_args=extra_link_args,
        extra_compile_args=[]
    )])
)
