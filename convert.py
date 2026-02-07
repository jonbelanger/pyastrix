from PIL import Image
import numpy as np
from astropy.io import fits

# Load PNG
img = Image.open("messier42_friedman_zs_61a_20260112.png").convert("RGB")
rgb = np.array(img)

# Convert to FITS-friendly format: (channels, y, x)
data = np.moveaxis(rgb, -1, 0)
    
with fits.open("messier42_friedman_zs_61a_20260112.wcs", ignore_missing_simple=True) as hdul:
    header = hdul[0].header.copy()

# Optional but polite
header['NAXIS']  = 3
header['NAXIS1'] = data.shape[2]
header['NAXIS2'] = data.shape[1]
header['NAXIS3'] = data.shape[0]
header['CTYPE3'] = 'RGB'

# Write FITS
fits.writeto(
    "messier42_friedman_zs_61a_20260112.fits",
    data,
    header,
    overwrite=True,
    output_verify='silentfix'
)