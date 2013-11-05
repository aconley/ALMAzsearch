from distutils.core import setup

import sys
major, minor1, minor2, release, serial = sys.version_info

if (major < 3) and (minor1 < 7):
    raise SystemExit("bethermin12_sim requires at least python 2.7")

setup(
    name="ALMAzsearch",
    version="0.1.0",
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
                'astropy (>0.3.0)']
)

