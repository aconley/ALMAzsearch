from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext

import sys
major, minor1, minor2, release, serial = sys.version_info

if (major < 3) and (minor1 < 7):
    raise SystemExit("ALMAzsearch requires at least python 2.7")

import numpy
ext_modules = [Extension("fnu", ["ALMAzsearch/fnu.pyx"],
                         include_dirs=[numpy.get_include()],
                         libraries=["m"])]

setup(
    name="ALMAzsearch",
    version="0.1.1",
    author="Alexander Conley",
    author_email="alexander.conley@colorado.edu",
    packages=["ALMAzsearch"],
    package_data={'ALMAzsearch':['resources/*tsys*']},
    license="GPL",
    description="ALMA z search simulator",
    classifiers=[
        "Development Status :: 3 - Alpha",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
    requires = ['numpy (>1.7.0)', 'scipy (>0.10.0)', 
                'astropy (>0.3.0)', 'cython (>0.11.0)'],
   cmdclass = {'build_ext': build_ext},
    ext_modules = ext_modules

)

