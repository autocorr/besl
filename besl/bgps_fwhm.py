"""
=========
BGPS FWHM
=========

Compute BGPS FWHM sizes and fluxes.

"""
import numpy as np
from scipy.ndimage import (binary_closing, watershed_ift)
import catalog
import image


# TODO
# Create maps of FWHM masks for evo matching


class FwhmSet(object):
    """
    Set of FWHM Solvers.

    Examples
    --------
    >>> fs = FwhmSet()
    >>> fs.solve_all()
    >>> fs.write()
    """

    def __init__(self, df=None, v=210, verbose=False):
        """
        Parameters
        ----------
        df : pd.DataFrame
        v : number, default 210
            BGPS version number to use for images
        verbose : bool, default False
            Print current clump
        """
        if df is None:
            df = catalog.read_cat('bgps_v210').set_index('v210cnum')
        self.df = df
        self.verbose = verbose
        self.fluxlib = image.BgpsLib(exten='map20', v=v)
        self.rindlib = image.BgpsLib(exten='labelmask', v=v)
        self.rmslib = image.BgpsLib(exten='noisemap20', v=v)
        print ':: Reading in Images'
        self.fluxlib.read_images()
        self.rindlib.read_images()
        self.rmslib.read_images()
        self.solvers = {}
        # DataFrame
        cols = ['npix', 'peak_flux', 'err_peak_flux', 'sangle', 'eqangle',
                'fwhm_npix', 'fwhm_sangle', 'fwhm_eqangle', 'fwhm_flux',
                'err_fwhm_flux']
        for col in cols:
            self.df[col] = np.nan

    def solve_all(self):
        """
        Solve for each clump and assign to the BGPS DataFrame.
        """
        print ':: Run Solver'
        for cnum in self.df.index:
            # Run solver
            field = self.df.loc[cnum, 'field']
            flux_hdu = self.fluxlib.get_hdu(field)
            rind_hdu = self.rindlib.get_hdu(field)
            rms_hdu = self.rmslib.get_hdu(field)
            solver = FwhmSolver(cnum, flux_hdu, rind_hdu, rms_hdu)
            solver.solve()
            # Assign parameters
            self.df.loc[cnum, 'npix'] = solver.npix
            self.df.loc[cnum, 'peak_flux'] = solver.peak
            self.df.loc[cnum, 'err_peak_flux'] = solver.err_peak
            self.df.loc[cnum, 'sangle'] = solver.sangle
            self.df.loc[cnum, 'eqangle'] = solver.eqangle
            self.df.loc[cnum, 'fwhm_npix'] = solver.fwhm_npix
            self.df.loc[cnum, 'fwhm_sangle'] = solver.fwhm_sangle
            self.df.loc[cnum, 'fwhm_eqangle'] = solver.fwhm_eqangle
            self.df.loc[cnum, 'fwhm_flux'] = solver.fwhm_flux
            self.df.loc[cnum, 'err_fwhm_flux'] = solver.err_fwhm_flux
            self.solvers[cnum] = solver
            if self.verbose & (cnum % 100 == 0):
                print '-- ', cnum
        fd = FwhmDeconvolve(self.df.copy())
        self.df = fd.df

    def write(self, outname='bgps_fwhm'):
        """
        Parameters
        ----------
        outname : str, default 'bgps_fwhm'
        """
        self.df.to_csv(outname + '.csv')


class FwhmSolver(object):
    """
    Compute the FWHM properties of a BGPS source.

    Examples
    --------
    >>> fs = FwhmSolver(5000, flux_hdu, rind_hdu)
    >>> fs = fs.solve()
    >>> fs.fwhm_flux
    0.9475
    """
    pixel_size = 7.2  # arcsec

    def __init__(self, cnum, flux_hdu, rind_hdu, rms_hdu):
        """
        Parameters
        ----------
        cnum : number
        flux : astropy.io.fits.Hdu
        rind : astropy.io.fits.Hdu
        """
        self.cnum = cnum
        self.ppbeam = flux_hdu[0].header['PPBEAM']
        self.flux = flux_hdu[0].data
        self.rind = rind_hdu[0].data
        self.rms = rms_hdu[0].data
        assert self.flux.shape == self.rind.shape == self.rms.shape

    def _clip_maps(self):
        """
        Clip the maps to the size of the labelmasks to save computation cost.
        """
        pix = np.argwhere(self.rind == self.cnum).T
        pix0 = np.unique(np.sort(pix[0]))
        pix1 = np.unique(np.sort(pix[1]))
        self.flux = self.flux.take(pix0, axis=0).take(pix1, axis=1)
        self.rind = self.rind.take(pix0, axis=0).take(pix1, axis=1)
        self.rms = self.rms.take(pix0, axis=0).take(pix1, axis=1)

    def _mask_on_rind(self):
        self.mflux = np.copy(self.flux)
        self.mflux[self.rind != self.cnum] = 0
        self.mrms = np.copy(self.rms)
        self.mrms[self.rind != self.cnum] = 0

    def _sample_peak(self):
        """
        Sample the peak flux value from the map in a 5-point average to
        calculated the half-max value. Uses uniform weight between the points.

        Returns
        -------
        fwhm : number
        """
        # Incase max on edge, pad edges (repeats edge value)
        mfpad = np.pad(self.mflux, pad_width=1, mode='edge')
        mrpad = np.pad(self.mrms, pad_width=1, mode='edge')
        self.pmax1, self.pmax0 = np.argwhere(mfpad == self.mflux.max())[0]
        # Sample five-point pattern
        pix0 = np.array([0, -1, 1, 0, 0]) + self.pmax0
        pix1 = np.array([0, 0, 0, -1, 1]) + self.pmax1
        self.peak = mfpad[pix1, pix0].mean()
        self.err_peak = np.sqrt(np.sum(mrpad[pix1, pix0]**2) / self.ppbeam +
                                (0.06 * mfpad[pix1, pix0].sum() / self.ppbeam)**2) / 5.
        self.fwhm = self.peak / 2.
        # Remove offset from padded array
        self.max0 = self.pmax0 - 1
        self.max1 = self.pmax1 - 1

    def _mask_on_fwhm(self):
        self.brind = np.copy(self.rind)
        # Mask the pixels below the threshold
        self.brind[self.rind != self.cnum] = 0
        self.brind[self.rind == self.cnum] = 1
        self.brind[self.mflux < self.fwhm] = 0
        # Closing to dilate and erode the masked pixels
        self.fwhm_mask = binary_closing(self.brind)
        # Watershed to select pixels contiguous with the peak pixel
        markers = np.zeros(self.fwhm_mask.size).reshape(self.fwhm_mask.shape)
        markers[self.max1, self.max0] = 2
        markers[~self.fwhm_mask] = 1
        self.fwhm_mask = watershed_ift(self.fwhm_mask.astype(np.uint8),
                                       markers.astype(np.int16))
        # Remask where watershed
        self.fwhm_mask[self.fwhm_mask != 2] = 0
        self.fwhm_mask[self.fwhm_mask == 2] = 1

    def _get_props(self):
        # Full
        self.npix = self.rind[self.rind == self.cnum].size  # pixels
        self.sangle = self.npix * self.pixel_size**2  # arcsec^2
        self.eqangle = np.sqrt(self.sangle / np.pi)  # arcsec
        # Peak
        self.peak /= self.ppbeam
        # FWHM
        fwhmv = self.mflux[self.fwhm_mask == 1]
        efwhmv = self.mrms[self.fwhm_mask == 1]
        self.fwhm_npix = fwhmv.size
        self.fwhm_sangle = self.fwhm_npix * self.pixel_size**2  # arcsec^2
        self.fwhm_eqangle = np.sqrt(self.fwhm_sangle / np.pi)  # arcsec
        self.fwhm_flux = fwhmv.sum() / self.ppbeam  # Jy
        self.err_fwhm_flux = np.sqrt(np.sum(efwhmv**2) / self.ppbeam +
                                     (0.06 * fwhmv.sum() / self.ppbeam)**2)  # Jy

    def solve(self):
        self._clip_maps()
        self._mask_on_rind()
        self._sample_peak()
        self._mask_on_fwhm()
        self._get_props()


class FwhmDeconvolve(object):
    theta_mb = 33 / 2.

    def __init__(self, df=None):
        if df is None:
            df = catalog.read_cat('bgps_v210_fwhm').set_index('v210cnum')
        self.df = df
        self._calc_values()

    def _calc_values(self):
        self.df['eqangled'] = self.df['eqangle'].apply(self.beam_subtract)
        self.df['fwhm_eqangled'] = self.df['fwhm_eqangle'].apply(self.beam_subtract)
        self.df['sangled'] = np.pi * self.df['eqangled']**2
        self.df['fwhm_sangled'] = np.pi * self.df['fwhm_eqangled']**2
        self.df['fwhm_sangled_ratio'] = self.df['fwhm_sangled'] / self.df['sangled']

    def beam_subtract(self, theta_eq):
        if theta_eq > self.theta_mb:
            return (theta_eq**2 - self.theta_mb**2)**(0.5)
        else:
            return np.nan


