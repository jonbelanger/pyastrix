from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy.coordinates import Angle
from astropy.wcs.utils import proj_plane_pixel_scales
import numpy as np



class FitsMeta:
    """
    Handles FITS file loading, WCS extraction, and geometry for annotation.
    Automatically detects celestial WCS axes for standard and PixInsight-style RGB FITS.
    """
    def __init__(self, filename):
        self.filename = filename
        self.hdulist = fits.open(filename)
        hdu_index, celestial_axes = self._find_celestial_hdu()
        self.header = self.hdulist[hdu_index].header
        self.data = self.hdulist[hdu_index].data
        self.no_valid_wcs = celestial_axes is None
        self.celestial_axes = celestial_axes
        self._init_data()
        self._init_wcs()
        self._compute_geometry()
        self._extract_obs_time()

    def _find_celestial_hdu(self):
        """
        Returns (hdu_index, celestial_axes) where celestial_axes is (1,2) for standard or (2,3) for PixInsight-style.
        If no valid WCS is found, returns (0, None).
        """
        for idx, hdu in enumerate(self.hdulist):
            header = hdu.header
            ctype1 = header.get('CTYPE1', '').upper()
            ctype2 = header.get('CTYPE2', '').upper()
            ctype3 = header.get('CTYPE3', '').upper()
            # Standard: celestial axes are 1/2
            if (ctype1.startswith('RA') or ctype1.startswith('GLON')) and (ctype2.startswith('DEC') or ctype2.startswith('GLAT')):
                return idx, (1, 2)
            # PixInsight-style: celestial axes are 2/3
            if (ctype2.startswith('RA') or ctype2.startswith('GLON')) and (ctype3.startswith('DEC') or ctype3.startswith('GLAT')):
                return idx, (2, 3)
        return 0, None

    def _init_data(self):
        if self.data.ndim == 3:
            self.data_2d = self.data[0, :, :]
        else:
            self.data_2d = self.data
        self.ny, self.nx = self.data_2d.shape

    def _init_wcs(self):
        # For 3D (RGB) images, use correct axes for celestial WCS
        if self.data.ndim == 3:
            if self.celestial_axes == (2, 3):
                self.wcs = WCS(self.header, naxis=2, key=' ')
            else:
                self.wcs = WCS(self.header, naxis=2)
        else:
            self.wcs = WCS(self.header)


    def _extract_obs_time(self):
        """Extract observation time from FITS header."""
        try:
            date_obs = self.header.get('DATE-OBS', None)
            exptime = self.header.get('EXPTIME', 0.0)
            if date_obs:
                start_time = Time(date_obs)
                self.obs_time = start_time + exptime / 2.0 * u.s
            else:
                self.obs_time = Time.now()
        except Exception:
            self.obs_time = Time.now()

    def world_to_pixel(self, skycoord):
        return self.wcs.world_to_pixel(skycoord)

    def _init_data(self):
        if self.data.ndim == 3:
            self.data_2d = self.data[0, :, :]
        else:
            self.data_2d = self.data

        self.ny, self.nx = self.data_2d.shape

    def _init_wcs(self):
        # For 3D (RGB) images, use correct axes for celestial WCS
        if self.data.ndim == 3:
            # Default: use last two axes, but override if PixInsight-style
            naxis = 2
            if hasattr(self, 'celestial_axes') and self.celestial_axes == (2, 3):
                #print("[DEBUG] Constructing WCS with naxis=2 (axes 2 and 3, PixInsight-style)")
                # Use only axes 2 and 3 for WCS
                self.wcs = WCS(self.header, naxis=2, key=' ')
                # Astropy will use the last two axes by default, which matches axes 2/3
            else:
                #print("[DEBUG] Constructing WCS with naxis=2 (last two axes, standard)")
                self.wcs = WCS(self.header, naxis=2)
        else:
            self.wcs = WCS(self.header)
        #print(f"[DEBUG] WCS summary: ctype={self.wcs.wcs.ctype}, naxis={self.wcs.wcs.naxis}")

    def _compute_geometry(self):
        corners_pix = np.array([
            [0, 0],
            [0, self.ny - 1],
            [self.nx - 1, 0],
            [self.nx - 1, self.ny - 1]
        ])

        # Pad pixel coordinates to match WCS axes (e.g., for 3D WCS)
        n_axes = self.wcs.wcs.naxis
        if corners_pix.shape[1] < n_axes:
            # Pad with reference pixel for extra axes
            pad_val = [self.header.get(f'CRPIX{i+1}', 0.0) for i in range(2, n_axes)]
            corners_pix = np.hstack([
                corners_pix,
                np.tile(pad_val, (corners_pix.shape[0], 1))
            ])

        corners_world = self.wcs.wcs_pix2world(corners_pix, 0)

        # Normalize angles to valid ranges
        ras = Angle(corners_world[:, 0] * u.deg).wrap_at(180 * u.deg)
        decs = Angle(corners_world[:, 1] * u.deg)
        
        # Clamp declinations to valid range [-90, 90]
        decs_clamped = np.clip(decs.deg, -90, 90)
        
        # Create normalized corner coordinates
        corners_norm = np.column_stack([ras.deg, decs_clamped])
        
        # Create SkyCoord objects from corners
        corners_coord = SkyCoord(
            ra=ras,
            dec=decs_clamped * u.deg
        )
        
        # Compute center as the mean of the SkyCoord objects
        center_ra = corners_coord.ra.mean()
        center_dec = corners_coord.dec.mean()
        
        self.center = SkyCoord(ra=center_ra, dec=center_dec)

        # Use normalized coordinates for top, left, and diagonal
        top = SkyCoord(
            ra=[corners_norm[0, 0], corners_norm[2, 0]] * u.deg,
            dec=[corners_norm[0, 1], corners_norm[2, 1]] * u.deg
        )

        left = SkyCoord(
            ra=[corners_norm[0, 0], corners_norm[1, 0]] * u.deg,
            dec=[corners_norm[0, 1], corners_norm[1, 1]] * u.deg
        )

        try:
            from astropy.wcs.utils import proj_plane_pixel_scales
            pixel_scales = proj_plane_pixel_scales(self.wcs)  # deg/pixel
            print(f"[INFO] Pixel scales: {pixel_scales} deg/pixel, Image size: {self.nx}x{self.ny}")
            mean_pixscale = np.sqrt(pixel_scales[0] * pixel_scales[1])  # deg/pixel
            diag_pix = np.sqrt(self.nx**2 + self.ny**2)
            diag_deg = diag_pix * mean_pixscale  # deg
            self.radius = (diag_deg / 2) * u.deg
            print(f"[INFO] Computed FOV: {diag_deg:.2f} deg, pixel scale: {mean_pixscale*3600:.3f} arcsec/pixel, radius: {self.radius}")
        except Exception as e:
            import traceback
            print(f"[WARNING] Could not compute pixel scale from WCS. Using fallback radius. Exception: {e}")
            traceback.print_exc()
            tl = SkyCoord(ra=corners_norm[0, 0]*u.deg, dec=corners_norm[0, 1]*u.deg)
            br = SkyCoord(ra=corners_norm[3, 0]*u.deg, dec=corners_norm[3, 1]*u.deg)
            diagonal = tl.separation(br)
            self.radius = diagonal / 2

    def _extract_obs_time(self):
        """Extract observation time from FITS header."""
        try:
            date_obs = self.header.get('DATE-OBS', None)
            exptime = self.header.get('EXPTIME', 0.0)
            if date_obs:
                start_time = Time(date_obs)
                self.obs_time = start_time + exptime / 2.0 * u.s
            else:
                self.obs_time = Time.now()
        except Exception:
            self.obs_time = Time.now()

    def world_to_pixel(self, skycoord):
        return self.wcs.world_to_pixel(skycoord)
