#!/usr/bin/env python

# Python 2.x
from distutils.command.build_py import build_py
from setuptools import setup
from besl import __version__ as version


setup(name='besl',
      version=version,
      description='besl - Brian Svobodda Library, a module of astronomy '
                  'related routines',
      author=['Brian Svoboda'],
      author_email=['svobodb@email.arizona.edu'],
      url='https://github.com/autocorr/besl',
      packages=['besl'],
      provides=['besl'],
      cmdclass={'build_py': build_py},
      keywords=['Scientific/Engineering', 'Astronomy'],
      classifiers=["Development Status :: 2 - Pre-alpha",
                   "Intended Audience :: Science/Research",
                   "Natural Language :: English",
                   "Operating System :: Unix",
                   "Programming Language :: Python",
                   "License :: OSI Approved :: "
                   "GNU General Public Licence v3 (GPLv3)"],
      requires=['numpy', 'scipy', 'pandas', 'matplotlib', 'astropy',
                'photutils', 'pyephem'])
