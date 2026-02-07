# simbad_catalog.py
from astroquery.simbad import Simbad
from astropy.coordinates import SkyCoord
from astropy.table import Table
import astropy.units as u
import numpy as np
from .catalog_base import CatalogQuery

import re

messier_pattern = re.compile(r"^M\d{1,3}$")  # M followed by 1–3 digits

class SimbadCatalog(CatalogQuery):
    """
    SIMBAD catalog adapter returning only "interesting" objects.

    Returns standardized columns:
        - source_id: string, either HD/TYC star or Messier/NGC/PGC/Abell token
        - ra: degrees
        - dec: degrees
    """

    def __init__(self, max_results=50000):
        self.max_results = max_results
        self.simbad = Simbad()
        self.simbad.TIMEOUT = 60
        self.simbad.reset_votable_fields()
        self.simbad.add_votable_fields(
            "main_id",
            "ra",
            "dec",
            "otype",
            "ids"
        )

        # Object type groups
        self.star_types = set([
            "*", "V*", "TT*", "Y*O", "Or*", "Ae*", "Be*", "RG*", "WD*", "PM*", "WD*", "RG*", "AB*", "LP*"
        ])
        self.ism_types = set([
            "ISM", "SFR", "HII", "Cld", "RNe", "MoC", "HVC", "cor", "SNR", "PN", "SN*"
        ])
        self.cluster_types = set([
            "Cl*", "OpC", "GlC", "As*"
        ])
        self.galaxy_types = set([
            "G", "EmG", "AGN", "SyG", "QSO", "SBG", "H2G", "Sy1", "Sy2", "bCG", "LSB"
        ])
        # Tokens for labeling non-stellar objects
        self.deep_sky_tokens = ["M", "NGC", "PGC", "Abell", "LEDA"]

    def query(self, center, radius):

        result = self.simbad.query_region(center, radius=radius)
        if result is None or len(result) == 0:
            return Table(names=["source_id", "ra", "dec"], dtype=[str, float, float])

        if len(result) > self.max_results:
            result = result[:self.max_results]

        table = self._standardize(result)
        table = self._filter_interesting(table)
        return table

    def _standardize(self, table):
        out = Table()
        colnames = [c.lower() for c in table.colnames]

        # RA/Dec
        coords = SkyCoord(table["ra"], table["dec"], unit=(u.hourangle, u.deg))
        out["ra"] = coords.ra.deg
        out["dec"] = coords.dec.deg

        # raw source_id
        out["source_id_raw"] = table["main_id"].astype(str)

        # Keep otype and IDS for filtering
        out["otype"] = table["otype"]
        out["IDS"] = table["ids"]

        return out

    def _filter_interesting(self, table):
        mask = np.zeros(len(table), dtype=bool)
        final_ids = []

        deep_sky_priority = ["M", "NGC", "Abell", "PGC", "LEDA"]

        for i, row in enumerate(table):
            otype = str(row["otype"]).strip()

            # check stellar types
            is_star = any(otype.startswith(t) for t in self.star_types)
            # check ISM, clusters, galaxies
            is_other = otype in self.ism_types or otype in self.cluster_types or otype in self.galaxy_types

            ids_field = str(row.get("IDS", ""))
            ids_tokens = [id_.replace(" ", "") for id_ in ids_field.split("|") if id_]

            source_id = None

            # --- Stars: HD preferred over TYC ---
            if is_star:
                for prefix in ["HD"]:
                    for tok in ids_tokens:
                        if tok.startswith(prefix):
                            source_id = tok
                            break
                    if source_id:
                        break
                if source_id is None:
                    mask[i] = False
                    final_ids.append("")
                    continue
                mask[i] = True

            # --- Non-stars: deep sky objects with priority ---
            elif is_other:
                for prefix in deep_sky_priority:
                    for tok in ids_tokens:
                        if prefix == "M":
                            if messier_pattern.match(tok):
                                source_id = tok
                                break
                        else:
                            if tok.startswith(prefix):
                                source_id = tok
                                break
                    if source_id:
                        break
                if source_id is None:
                    mask[i] = False
                    final_ids.append("")
                    continue
                mask[i] = True

            else:
                mask[i] = False
                final_ids.append("")
                continue

            final_ids.append(source_id)

        filtered = table[mask]
        filtered["source_id"] = [final_ids[i] for i, m in enumerate(mask) if m]

        # drop internal columns
        for c in ["otype", "IDS", "source_id_raw"]:
            if c in filtered.colnames:
                filtered.remove_column(c)

        return filtered
