from astroquery.gaia import Gaia
from astropy.io import fits
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy.visualization import ImageNormalize, PercentileInterval, AsinhStretch
import astropy.units as u
import numpy as np
import matplotlib.pyplot as plt

# Load your FITS
filename = "messier42_friedman_zs_61a_20260112.fits"
hdulist = fits.open(filename)
header = hdulist[0].header
data = hdulist[0].data
print(data.shape)

# Pick the first plane if 3D
if data.ndim == 3:
    # Use only the last 2 axes for WCS
    wcs = WCS(hdulist[0].header, naxis=2)
    data_2d = data[0, :, :]  # first slice
else:
    wcs = WCS(hdulist[0].header)
    data_2d = data

data_2d = data[0, :, :]  # shape now (2022, 2865)

# pixel dimensions
ny, nx = data_2d.shape

# compute WCS coordinates for corners
corners_pix = np.array([
    [0, 0],
    [0, ny-1],
    [nx-1, 0],
    [nx-1, ny-1]
])
corners_world = wcs.wcs_pix2world(corners_pix, 0)  # 0-based

# compute center as mean of corners
center_ra = np.mean(corners_world[:, 0])
center_dec = np.mean(corners_world[:, 1])
center = SkyCoord(ra=center_ra*u.deg, dec=center_dec*u.deg)

# compute angular separation along each edge
top_edge = SkyCoord(ra=[corners_world[0,0], corners_world[2,0]]*u.deg,
                    dec=[corners_world[0,1], corners_world[2,1]]*u.deg)
left_edge = SkyCoord(ra=[corners_world[0,0], corners_world[1,0]]*u.deg,
                     dec=[corners_world[0,1], corners_world[1,1]]*u.deg)

# angular size along x (cols) and y (rows)
x_size = top_edge[0].separation(top_edge[1])
y_size = left_edge[0].separation(left_edge[1])

# radius = half of shortest side
radius = 0.5 * min(x_size, y_size)
print(f"Circle radius to fit inside shortest dimension: {radius.to(u.deg)}")

query = f"""
SELECT
    source_id,
    ra,
    dec,
    parallax,
    parallax_error,
    pmra,
    pmdec,
    phot_g_mean_mag,
    bp_rp,
    teff_gspphot,
    3261.56 / parallax AS distance_ly,
    (phot_g_mean_mag + 5*LOG10(parallax) - 10) AS abs_mag
FROM gaiadr3.gaia_source
WHERE CONTAINS(
    POINT('ICRS', ra, dec),
    CIRCLE('ICRS', {center_ra}, {center_dec}, {radius.to(u.deg).value})
    ) = 1
    AND parallax / parallax_error > 5
    AND phot_g_mean_mag <= 10
ORDER BY phot_g_mean_mag ASC
"""

job = Gaia.launch_job(query)
results = job.get_results()

print(results)

# --- convert Gaia RA/Dec to pixel coordinates ---
# Convert to plain float arrays
ra_vals = np.array(results['ra'], dtype=float)
dec_vals = np.array(results['dec'], dtype=float)

gaia_coords = SkyCoord(ra=ra_vals*u.deg, dec=dec_vals*u.deg)
x_pix, y_pix = wcs.world_to_pixel(gaia_coords)

# Suppose data.shape = (3, ny, nx) for a 3-layer FITS cube
if data.ndim == 3:
    ny, nx = data.shape[1], data.shape[2]
    rgb = np.zeros((ny, nx, 3), dtype=np.float32)
    
    for i in range(3):
        # normalize each plane independently
        norm = ImageNormalize(data[i,:,:],
                              interval=PercentileInterval(99.5),  # clip top 0.5%
                              stretch=AsinhStretch())             # asinh stretch
        rgb[:,:,i] = norm(data[i,:,:])
else:
    # single-layer grayscale -> RGB
    norm = ImageNormalize(data,
                          interval=PercentileInterval(99.5),
                          stretch=AsinhStretch())
    rgb = np.stack([norm(data)]*3, axis=-1)


# Now rgb is in [0,1], ready for plotting with annotations
plt.figure(figsize=(12,12))
plt.imshow(rgb, origin='lower')

# --- START OF ADJUSTTEXT BLOCK ---
from adjustText import adjust_text
texts = []

for i in range(len(results)):
    x = x_pix[i]
    y = y_pix[i]

    # Draw the hollow circle marker
    plt.plot(x, y,
         marker='o',
         markersize=5,
         markerfacecolor='none',
         markeredgecolor='red',
         markeredgewidth=.3)

    # Create the label
    distance_ly = float(results['distance_ly'][i])

    label = f"G={results['phot_g_mean_mag'][i]:.2f}, {distance_ly:.1f} ly"

    # Add text but let adjustText move it
    texts.append(plt.text(x, y + 10, label, color='red', fontsize=8,
                          fontfamily='sans-serif', ha='center', va='bottom'))

# Adjust all labels to avoid overlap
adjust_text(texts, only_move={'points':'y', 'text':'y'},
            arrowprops=dict(arrowstyle='-', color='red', lw=0.5))
# --- END OF ADJUSTTEXT BLOCK ---

plt.axis('off')
plt.savefig("messier42_annotated.png", dpi=200, bbox_inches='tight', pad_inches=0)
plt.close()