from astroquery.gaia import Gaia
from .catalog_base import CatalogQuery

class GaiaCatalog(CatalogQuery):
    def __init__(self, mag_limit=10, parallax_snr=5):
        self.mag_limit = mag_limit
        self.parallax_snr = parallax_snr

    def query(self, center, radius):
        query = f"""
        SELECT
            source_id,
            ra,
            dec,
            parallax,
            phot_g_mean_mag,
            bp_rp,
            3261.56 / parallax AS distance_ly
        FROM gaiadr3.gaia_source
        WHERE CONTAINS(
            POINT('ICRS', ra, dec),
            CIRCLE('ICRS', {center.ra.deg}, {center.dec.deg}, {radius.to('deg').value})
        ) = 1
        AND parallax / parallax_error > {self.parallax_snr}
        AND phot_g_mean_mag <= {self.mag_limit}
        ORDER BY phot_g_mean_mag ASC
        """

        job = Gaia.launch_job(query)
        return job.get_results()
