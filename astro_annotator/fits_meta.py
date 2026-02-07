from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
import astropy.units as u
import numpy as np


class FitsMeta:
    def __init__(self, filename):
        self.filename = filename
        self.hdulist = fits.open(filename)
        self.header = self.hdulist[0].header
        self.data = self.hdulist[0].data

        self._init_data()
        self._init_wcs()
        self._compute_geometry()

    def _init_data(self):
        if self.data.ndim == 3:
            self.data_2d = self.data[0, :, :]
        else:
            self.data_2d = self.data

        self.ny, self.nx = self.data_2d.shape

    def _init_wcs(self):
        self.wcs = WCS(self.header, naxis=2)

    def _compute_geometry(self):
        corners_pix = np.array([
            [0, 0],
            [0, self.ny - 1],
            [self.nx - 1, 0],
            [self.nx - 1, self.ny - 1]
        ])

        corners_world = self.wcs.wcs_pix2world(corners_pix, 0)

        self.center = SkyCoord(
            ra=np.mean(corners_world[:, 0]) * u.deg,
            dec=np.mean(corners_world[:, 1]) * u.deg
        )

        top = SkyCoord(
            ra=[corners_world[0, 0], corners_world[2, 0]] * u.deg,
            dec=[corners_world[0, 1], corners_world[2, 1]] * u.deg
        )

        left = SkyCoord(
            ra=[corners_world[0, 0], corners_world[1, 0]] * u.deg,
            dec=[corners_world[0, 1], corners_world[1, 1]] * u.deg
        )

        #self.radius = 0.5 * min(
        #    top[0].separation(top[1]),
        #    left[0].separation(left[1])
        #)
        # full image diagonal (top-left to bottom-right)
        tl = SkyCoord(ra=corners_world[0, 0]*u.deg, dec=corners_world[0, 1]*u.deg)
        br = SkyCoord(ra=corners_world[3, 0]*u.deg, dec=corners_world[3, 1]*u.deg)
        diagonal = tl.separation(br)

        # radius = half diagonal
        self.radius = diagonal / 2

    def world_to_pixel(self, skycoord):
        return self.wcs.world_to_pixel(skycoord)
