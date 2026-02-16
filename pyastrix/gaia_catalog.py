from astroquery.gaia import Gaia
from .catalog_base import CatalogQuery
import numpy as np

class GaiaCatalog(CatalogQuery):
    def __init__(self, mag_limit=10, parallax_snr=5):
        self.mag_limit = mag_limit
        self.parallax_snr = parallax_snr
        self.name = "gaia"

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
        table = job.get_results()

        # populate a single display label per row so renderer can stay catalog-agnostic
        labels = []
        for row in table:
            parts = []
            sid = str(row.get("source_id", ""))
            if sid and any(c.isalpha() for c in sid):
                parts.append(sid)

            mag = row.get("phot_g_mean_mag", None)
            if mag is not None and not np.isnan(mag):
                parts.append(f"G={mag:.2f}")

            distance = row.get("distance_ly", None) or row.get("distance", None)
            if distance is not None and not np.isnan(distance):
                parts.append(f"{float(distance):.0f} ly")

            labels.append(", ".join(parts) if parts else "")

        table["display_label"] = labels
        return table
