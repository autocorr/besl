"""
=========
BGPS FWHM
=========

Compute BGPS FWHM sizes and fluxes.

"""
import numpy as np
from scipy.ndimage import binary_closing
from besl.catalog import read_bgps
from besl.image import BgpsLib


## TODO Notes
# for cnum in bgps.index:
# * get rind
# * get image
#   mask rind on cnum
#   mask flux on (cnum masked rind)
#   get pos of peak pixel
#   sample and weighted average of peak pixel
#   copy mask where flux less than avg max
#   grow the flux mask by closure
#   sample the flux mask
#   out sum flux
#   out number of pixels for area
#   write to table
#   (add pixels back to full rind and save image per field for later use)


class FwhmSet(object):
    """
    Set of FWHM Solvers.
    """
    def __init__(self, v=201, verbose=False):
        self.bgps = read_bgps(v=v).set_index('cnum')
        self.fluxlib = BgpsLib(exten='map20', v=v)
        self.rindlib = BgpsLib(exten='labelmask', v=v)
        self.solvers = {}
        # DataFrame
        self.bgps['fwhm_flux'] = np.nan
        self.bgps['fwhm_angle'] = np.nan
        self.bgps['fwhm_sangle'] = np.nan

    def solve_all(self):
        for cnum in self.bgps.index:
            field = self.bgps.loc[cnum, 'field']
            flux = self.fluxlib.get_img(field)
            rind = self.rindlib.get_lib(field)
            solver = FwhmSolver(cnum, flux, rind)
            solver.solve()
            self.bgps.loc[cnum, 'fwhm_flux'] = solver.fwhm_flux
            self.bgps.loc[cnum, 'fwhm_angle'] = solver.fwhm_angle
            self.bgps.loc[cnum, 'fwhm_sangle'] = solver.fwhm_sangle
            self.solvers[cnum] = solver
            if verbose:
                print cnum

    def write(self, outname='bgps_fwhm'):
        self.bgps.to_csv(outname)


class FwhmSolver(object):
    """
    Compute the FWHM properties of a BGPS source.
    """
    def __init__(self, cnum, flux, rind):
        """
        Parameters
        ----------
        cnum : number
        flux : astropy.io.fits.Hdu
        rind : astropy.io.fits.Hdu
        """
        assert flux.shape == rind.shape
        self.cnum = cnum
        self.flux = flux
        self.rind = rind

    def _clip_maps(self):
        """
        Clip the maps to the size of the labelmasks to save computation cost.
        """
        pix = np.argwhere(self.rind == self.cnum).T[::-1]
        pix0 = np.unique(np.sort(pix[0]))
        pix1 = np.unique(np.sort(pix[1]))
        import ipdb; ipdb.set_trace()
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
        max1, max0 = np.argwhere(self.mflux == self.mflux.max())[0]
        pix0 = np.array([0, -1, 1, 0, 0]) + max0
        pix1 = np.array([0, 0, 0, -1, 1]) + max1
        self.peak = np.self.mflux[pix1, pix0].mean()
        self.fwhm = self.peak / 2.
        return self.fwhm

    def _mask_on_fwhm(self):
        self.brind = np.copy(self.rind)
        self.brind[self.rind != self.cnum] = 0
        self.brind[self.rind == self.cnum] = 1
        self.brind[self.mflux < self.fwhm] = 0
        self.fwhm_mask = binary_closing(self.brind)

    def _get_props(self):
        fwhm_vals = self.mflux[self.fwhm_mask]
        self.fwhm_flux = fwhm_vals.sum()
        self.fwhm_npix = fwhm_vals.size
        self.fwhm_sangle = self.fwhm_npix * 7.2**2  # in arcsec^2
        self.fwhm_angle = np.sqrt(self.sangle / np.pi)

    def solve(self):
        self._clip_maps()
        self._mask_on_rind()
        self._sample_peak()
        self._mask_on_fwhm()
        self._get_props


