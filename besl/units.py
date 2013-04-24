"""
===================
Units and Constants
===================

Astrophysical units and constants.

References
----------
[1] : http://ssd.jpl.nasa.gov/?constants
[2] : http://wikipedia.org
[3] : http://ned.ipac.caltech.edu/level5/Glossary/lang_formulae.html
"""

from __future__ import division
import numpy as _np

class Number(object):
    """
    Numerical constants class.
    """
    def __init__(self):
        # angle
        self.arcsec = _np.deg2rad(1 / (60 * 60.))
        self.arcmin = _np.deg2rad(1 / 60.)
        self.degree = _np.deg2rad(1.)
        self.hour = _np.deg2rad(360 / 24.)
        # composite
        self.FWHM = _np.sqrt(8 * _np.log(2.))
        self.HWHM = self.FWHM / 2.

class CGS(object):
    """
    Centimeter/gram/second unit system class.
    """
    def __init__(self):
        # fundamental
        self.c = 29979245800. # cm s^-1, [1]
        self.k = 1.3806504e-16 # erg K^-1, [2]
        self.G = 6.67428e-8 # cm^3 g^-1 s^-2, [2]
        self.h = 6.62606885e-27 # erg s, [2]
        self.hbar = self.h / (2 * _np.pi) # erg s, [2]
        self.amu = 1.660538782e-24 # g, [2]
        self.bohr_magneton_emu = 9.27400915e-21 # erg G^-1 (EMU), [2]
        self.bohr_radius = 5.2917720859e-9 # cm, [2]
        self.emass = 9.10938215e-28 # g, [2]
        self.fine_struct = 7.297352570e-3 # [2]
        # astronomy
        self.jday = 86400. # s, [1]
        self.jyear = 365.25 * self.jday # s, [1]
        self.au = 14959787070000. # cm, [1]
        self.pc = 648000 / _np.pi * self.au # cm, [2]
        self.kpc = 1e3 * self.pc
        self.Mpc = 1e6 * self.pc
        self.Gpc = 1e9 * self.pc

class MKS(object):
    """
    Meter/kilogram/second unit system class.
    """
    def __init__(self):
        # astronomy
        self.jday = 86400. # s, [1]
        self.jyear = 365.25 * self.jday # s, [1]
        self.au = 149597870700. # m, [1]
        self.pc = 648000 / _np.pi * self.au # m, [2]
        self.kpc = 1e3 * self.pc
        self.Mpc = 1e6 * self.pc
        self.Gpc = 1e9 * self.pc
        pass

class Astro(object):
    """
    Astrophysical constants class.
    """
    def __init__(self):
        # TODO
        #self.msol =
        #self.mearth =
        #self.mmoon =
        #self.mjup =
        pass

class Molecule(object):
    """
    Molecular constants class.
    """
    def __init__(self):
        # TODO mass, dipole moments, frequencies
        # fundamental
        self.amu = 1.660538782e-24 # g, [2]
        # mass
        #self.h2_m =
        #self.co_m =
        #self.oh_m =
        #self.cs_m =
        #self.cn_m =
        #self.h2o_m =
        #self.hcn_m =
        #self.nh3_m =
        #self.ccs_m =
        #self.h2co_m =
        # dipole moments
        # transition frequencies
        pass

# Initialize classes
num = Number()
cgs = CGS()
mks = MKS()
astro = Astro()
molec = Molecule()
