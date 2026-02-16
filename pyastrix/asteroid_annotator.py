# asteroid_annotator.py
from astroquery.imcce import Skybot
from astropy.coordinates import SkyCoord
from astropy.table import Table, vstack
import astropy.units as u
import numpy as np

from .catalog_base import CatalogQuery

class AsteroidAnnotator(CatalogQuery):
    """
    Query SkyBot for all asteroids visible in a FITS image (cone search).
    Returns standardized columns:
        - source_id: string, asteroid name or number
        - ra: degrees
        - dec: degrees
        - phot_g_mean_mag: apparent magnitude (V-band)
    """
    def __init__(self):
        self.name = "asteroids"

    def query(self, center, radius, obs_time=None):
        """
        Query all asteroids at a given position and time using SkyBot.
        Args:
            center: SkyCoord object with RA/Dec
            radius: search radius in degrees (astropy Quantity or float)
            obs_time: observation time (astropy.time.Time object)
        Returns:
            Astropy Table with asteroid positions and magnitudes
        """
        if obs_time is None:
            from astropy.time import Time
            obs_time = Time.now()

        # Ensure radius is a float in degrees
        if hasattr(radius, 'to'):
            radius_deg = float(radius.to(u.deg).value)
        else:
            radius_deg = float(radius)

        print(f"\n=== Asteroid Query (SkyBot) ===")
        print(f"Center: {center.ra.deg:.2f}, {center.dec.deg:.2f}")
        print(f"Search radius: {radius_deg:.2f} deg")
        print(f"Observation epoch: {obs_time.iso}")

        try:
            result = Skybot.cone_search(center, radius_deg * u.deg, obs_time)
            if result is None or len(result) == 0:
                print("No asteroids found in FOV.")
                return Table(names=["source_id", "ra", "dec", "phot_g_mean_mag"], dtype=[str, float, float, float])

            # Debug: print column names and first row
            print(f"[DEBUG] SkyBot result columns: {result.colnames}")
            if len(result) > 0:
                print(f"[DEBUG] First row: {result[0]}")

            # Standardize output, now with velocity columns
            out = Table(names=["source_id", "ra", "dec", "phot_g_mean_mag", "velocity_ra", "velocity_dec"], dtype=[str, float, float, float, float, float])
            for row in result:
                ra = row['RA'].to_value(u.deg) if hasattr(row['RA'], 'to_value') else float(row['RA'])
                dec = row['DEC'].to_value(u.deg) if hasattr(row['DEC'], 'to_value') else float(row['DEC'])
                mag = row['V'].value if hasattr(row['V'], 'value') else float(row['V']) if 'V' in row.colnames else None
                # Use RA_rate and DEC_rate (arcsec/hour), convert to arcsec/minute
                if 'RA_rate' in row.colnames and hasattr(row['RA_rate'], 'to_value'):
                    vel_ra = row['RA_rate'].to_value(u.arcsec/u.hour) / 60.0
                elif 'RA_rate' in row.colnames:
                    vel_ra = float(row['RA_rate']) / 60.0
                else:
                    vel_ra = 0.0
                if 'DEC_rate' in row.colnames and hasattr(row['DEC_rate'], 'to_value'):
                    vel_dec = row['DEC_rate'].to_value(u.arcsec/u.hour) / 60.0
                elif 'DEC_rate' in row.colnames:
                    vel_dec = float(row['DEC_rate']) / 60.0
                else:
                    vel_dec = 0.0
                out.add_row([
                    str(row['Number']) if row['Number'] else str(row['Name']),
                    ra,
                    dec,
                    mag,
                    vel_ra,
                    vel_dec
                ])
            print(f"Found {len(out)} asteroids in FOV.")
            return out
        except Exception as e:
            print(f"Error querying SkyBot: {e}")
            import traceback
            traceback.print_exc()
            return Table(names=["source_id", "ra", "dec", "phot_g_mean_mag"], dtype=[str, float, float, float])
