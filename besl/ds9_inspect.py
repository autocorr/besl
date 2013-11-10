#!/usr/bin/env python
# encoding: utf-8

import os
import ds9
from astropy import wcs
from astropy.io import fits
import numpy as np
import pandas as pd


class Inspector(object):
    rind_contour_color = 'yellow'

    def __init__(self, cnum):
        self.cnum = cnum
        self.d = ds9.ds9()

    def _get_rind(self)
        pass

    def show_rind(self):
        self._get_find()
        self.d.set('frame 2')
        self.d.set('contour limits 0 1')
        self.d.set('contour nlevels 1')
        self.d.set('contour smooth 1')
        self.d.set('contour method block')
        self.d.set('contour')
        self.d.set('contour copy')
        self.d.set('frame 1')
        self.d.set('contour paste {0}'.format(self.rind_contour_color))
        self.d.set('frame delete 2')


class HiGalInspector(Inspector):
    """
    Visually inspect HiGal cutouts in DS9.
    """
    # Directories and file names
    self.root = '/mnt/eld_data/HiGal/'
    self.data = os.path.join(self.root, 'dat_files')
    self.data_file = os.path.join(self.data,
        'source_{cnum}_blue_svoboda_photall_err.dat')
    self.img = os.path.join(self.root, 'img_files')
    self.img_file = os.path.join(self.img,
        '{img_type}_source_{cnum}_blue_svoboda.fits')

    def __init__(self, cnum):
        super(HigalInspector, self).__init__(cnum)
        self.cnum = cnum


    def get_img(self, cnum, img_type):
        self.filen = self.img_file.format(cnum=cnum, img_type=img_type)
        self.img = fits.open(self.filen)
        self.d.set_pyfits(self.img)

    def view(self, cnum, img_type):
        self.show_rind(cnum)
        self.get_img(cnum, img_type)


