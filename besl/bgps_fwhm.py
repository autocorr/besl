"""
=========
BGPS FWHM
=========

Compute BGPS FWHM sizes and fluxes.

"""
import numpy as np
from scipy.ndimage import (binary_closing, watershed_ift)
from besl.catalog import read_cat
from besl.image import BgpsLib


# TODO
# Create maps of FWHM masks for evo matching


class FwhmSet(object):
    """
    Set of FWHM Solvers.
    """
    bgps = read_cat('bgps_v210').set_index('v210cnum')

    def __init__(self, v=210, verbose=False):
        """
        Parameters
        ----------
        v : number, default 210
            BGPS version number to use for images
        verbose : bool, default False
            Print current clump
        """
        self.verbose = verbose
        self.fluxlib = BgpsLib(exten='map20', v=v)
        self.rindlib = BgpsLib(exten='labelmask', v=v)
        self.fluxlib.read_images()
        self.rindlib.read_images()
        self.solvers = {}
        # DataFrame
        self.bgps['npix'] = np.nan
        self.bgps['fwhm_flux'] = np.nan
        self.bgps['fwhm_angle'] = np.nan
        self.bgps['fwhm_sangle'] = np.nan
        self.bgps['fwhm_npix'] = np.nan

    def solve_all(self):
        """
        Solve for each clump and assign to the BGPS DataFrame.
        """
        for cnum in self.bgps.index:
            # Run solver
            field = self.bgps.loc[cnum, 'field']
            flux_hdu = self.fluxlib.get_hdu(field)
            rind_hdu = self.rindlib.get_hdu(field)
            solver = FwhmSolver(cnum, flux_hdu, rind_hdu)
            solver.solve()
            # Assign parameters
            self.bgps.loc[cnum, 'fwhm_flux'] = solver.fwhm_flux
            self.bgps.loc[cnum, 'fwhm_angle'] = solver.fwhm_angle
            self.bgps.loc[cnum, 'fwhm_sangle'] = solver.fwhm_sangle
            self.bgps.loc[cnum, 'fwhm_npix'] = solver.fwhm_npix
            self.bgps.loc[cnum, 'npix'] = solver.npix
            self.solvers[cnum] = solver
            if self.verbose & (cnum % 100 == 0):
                print '-- ', cnum

    def write(self, outname='bgps_v210_fwhm'):
        """
        Parameters
        ----------
        outname : str, default 'bgps_fwhm'
        """
        self.bgps.to_csv(outname + '.csv')


class FwhmSolver(object):
    """
    Compute the FWHM properties of a BGPS source.

    Examples
    --------
    >>> fs = FwhmSolver(5000, flux_hdu, rind_hdu)
    >>> fs = fs.solve()
    >>> fs.fwhm_flux
    0.904
    """
    def __init__(self, cnum, flux_hdu, rind_hdu):
        """
        Parameters
        ----------
        cnum : number
        flux : astropy.io.fits.Hdu
        rind : astropy.io.fits.Hdu
        """
        self.cnum = cnum
        self.flux_hdu = flux_hdu
        self.rind_hdu = rind_hdu
        self.ppbeam = flux_hdu[0].header['PPBEAM']
        self.flux = flux_hdu[0].data
        self.rind = rind_hdu[0].data
        assert flux.shape == rind.shape

    def _clip_maps(self):
        """
        Clip the maps to the size of the labelmasks to save computation cost.
        """
        pix = np.argwhere(self.rind == self.cnum).T
        pix0 = np.unique(np.sort(pix[0]))
        pix1 = np.unique(np.sort(pix[1]))
        self.flux = self.flux.take(pix0, axis=0).take(pix1, axis=1)
        self.rind = self.rind.take(pix0, axis=0).take(pix1, axis=1)

    def _mask_on_rind(self):
        self.mflux = np.copy(self.flux)
        self.mflux[self.rind != self.cnum] = 0

    def _sample_peak(self):
        """
        Sample the peak flux value from the map in a 5-point average to
        calculated the half-max value. Uses uniform weight between the points.

        Returns
        -------
        fwhm : number
        """
        # Incase max on edge, pad edges (repeats edge value)
        mpad = np.pad(self.mflux, pad_width=1, mode='edge')
        self.pmax1, self.pmax0 = np.argwhere(mpad == self.mflux.max())[0]
        # Sample five-point pattern
        pix0 = np.array([0, -1, 1, 0, 0]) + self.pmax0
        pix1 = np.array([0, 0, 0, -1, 1]) + self.pmax1
        self.peak = mpad[pix1, pix0].mean()
        self.fwhm = self.peak / 2.
        # Remove offset from padded array
        self.max0 = self.pmax0 - 1
        self.max1 = self.pmax1 - 1
        return self.fwhm

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
        fwhm_vals = self.mflux[self.fwhm_mask == 1]
        self.npix = self.rind[self.rind == self.cnum].size
        self.fwhm_npix = fwhm_vals.size
        self.fwhm_sangle = self.fwhm_npix * 7.2**2  # in arcsec^2
        self.fwhm_angle = np.sqrt(self.fwhm_sangle / np.pi)  # in arcsec
        self.fwhm_flux = fwhm_vals.sum() / self.ppbeam  # in Jy

    def solve(self):
        self._clip_maps()
        self._mask_on_rind()
        self._sample_peak()
        self._mask_on_fwhm()
        self._get_props()


