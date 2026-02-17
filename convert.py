
import sys
from astropy.wcs import WCS
from PIL import Image
import numpy as np
from astropy.io import fits

def main():
    if len(sys.argv) != 4:
        print("Usage: python convert.py <input_image> <input_wcs_fits> <output_fits>")
        sys.exit(1)

    input_image = sys.argv[1]
    input_wcs_fits = sys.argv[2]
    output_fits = sys.argv[3]

    # Load image
    img = Image.open(input_image).convert("RGB")
    rgb = np.array(img)
    data = np.moveaxis(rgb, -1, 0)

    with fits.open(input_wcs_fits, ignore_missing_simple=True) as hdul:
        header = hdul[0].header.copy()

    # Optional but polite
    header['NAXIS']  = 3
    header['NAXIS1'] = data.shape[2]
    header['NAXIS2'] = data.shape[1]
    header['NAXIS3'] = data.shape[0]
    header['CTYPE3'] = 'RGB'

    # Write FITS
    fits.writeto(
        output_fits,
        data,
        header,
        overwrite=True,
        output_verify='silentfix'
    )

if __name__ == "__main__":
    main()