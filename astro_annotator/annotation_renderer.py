# annotation_renderer.py
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
import matplotlib.patheffects as path_effects

class AnnotationRenderer:
    """Render catalog annotations on a FITS/RGB image using a provided Axes.

    Changes made:
      - Groups nearby detections (pixel-based clustering) and merges labels
      - Supports multiple label parts per object (e.g. name, magnitude, distance)
      - Renders each label part in a source-specific color (hard-coded mapping)
    """

    def __init__(self, wcs):
        self.wcs = wcs

        # colorization removed: labels will be rendered as single combined strings

    def _build_entries(self, results, nx, ny):
        """Build normalized entries with coords, pixel positions and a single label."""
        entries = []
        for row in results:
            try:
                coord = SkyCoord(row["ra"], row["dec"], unit="deg")
            except Exception:
                continue

            x_pix, y_pix = self.wcs.world_to_pixel(coord)
            if not (0 <= x_pix < nx and 0 <= y_pix < ny):
                continue

            src = row.get("source", None)
            disp_label = row.get("display_label", None) or row.get("display_name", None)
            if disp_label:
                label = str(disp_label)
            else:
                parts = self._extract_parts(row, src)
                if not parts:
                    continue
                label = ", ".join(parts)

            entries.append({
                "coord": coord,
                "x": float(x_pix),
                "y": float(y_pix),
                "label": label,
                "source": src,
            })

        return entries

    def _cluster_entries(self, entries, cluster_arcsec, cluster_pixels):
        """Cluster entries into groups using angular OR pixel thresholds.

        Returns list of groups with centroid coords and member lists.
        """
        import numpy as _np

        n = len(entries)
        if n == 0:
            return []
        if n == 1:
            e = entries[0]
            return [{
                "coord_centroid": e["coord"],
                "x": e["x"],
                "y": e["y"],
                "members": [e],
            }]

        coords = [e["coord"] for e in entries]
        xs = _np.array([e["x"] for e in entries], dtype=float)
        ys = _np.array([e["y"] for e in entries], dtype=float)

        parent = list(range(n))

        def find(i):
            while parent[i] != i:
                parent[i] = parent[parent[i]]
                i = parent[i]
            return i

        def union(i, j):
            ri = find(i)
            rj = find(j)
            if ri == rj:
                return
            parent[rj] = ri

        for i in range(n):
            for j in range(i + 1, n):
                sep = coords[i].separation(coords[j]).arcsecond
                pix_dist = _np.hypot(xs[i] - xs[j], ys[i] - ys[j])
                if sep <= cluster_arcsec or pix_dist <= cluster_pixels:
                    union(i, j)

        comp = {}
        for i in range(n):
            r = find(i)
            comp.setdefault(r, []).append(entries[i])

        groups = []
        for members in comp.values():
            ras = [m["coord"].ra.deg for m in members]
            decs = [m["coord"].dec.deg for m in members]
            xs_m = [m["x"] for m in members]
            ys_m = [m["y"] for m in members]
            groups.append({
                "coord_centroid": SkyCoord(ra=_np.mean(ras) * u.deg, dec=_np.mean(decs) * u.deg),
                "x": float(_np.mean(xs_m)),
                "y": float(_np.mean(ys_m)),
                "members": members,
            })

        return groups

    def _compose_group_label(self, group):
        """Compose a single label for a clustered group.

        SIMBAD-sourced labels are ordered first. Duplicate labels are removed.
        """
        simbad_members = []
        other_members = []
        seen = set()
        for m in group["members"]:
            lbl = m.get("label", "")
            if not lbl or lbl in seen:
                continue
            seen.add(lbl)
            if m.get("source") and "simbad" in str(m.get("source")).lower():
                simbad_members.append(lbl)
            else:
                other_members.append(lbl)

        member_strings = simbad_members + other_members
        if not member_strings:
            return None
        return "; ".join(member_strings)

    def _extract_parts(self, row, source):
        """Return list of (text, color) for a row depending on source fields."""
        parts = []
        src_key = (source or "").lower()

        # Gaia-style fields
        if "gaia" in src_key:
            mag = row.get("phot_g_mean_mag", None)
            if mag is not None and not np.isnan(mag):
                parts.append(f"G={mag:.2f}")

            distance = row.get("distance_ly", None) or row.get("distance", None)
            if distance is not None and not np.isnan(distance):
                parts.append(f"{float(distance):.0f} ly")

            # prefer a textual source id if it contains letters
            sid = str(row.get("source_id", ""))
            if sid and any(c.isalpha() for c in sid):
                parts.insert(0, sid)

        # SIMBAD-like catalogs often provide a 'main_id' or 'name'
        elif "simbad" in src_key:
            name = row.get("source_id", None)
            if name:
                parts.append(str(name))

        else:
            # generic fallback: include any textual id, magnitude, distance
            sid = row.get("source_id", None) or row.get("id", None)
            if sid:
                sid_s = str(sid)
                if any(c.isalpha() for c in sid_s):
                    parts.append(sid_s)

            mag = row.get("phot_g_mean_mag", None)
            if mag is not None and not np.isnan(mag):
                parts.append(f"G={mag:.2f}")

        return parts

    def render(self, ax, results, rgb, cluster_arcsec=.1, cluster_pixels=1.0, y_offset_pixels=10, spacing_pixels=4):
        """Render combined annotations.

        - `results` should be an Astropy Table or iterable with fields
            including 'ra', 'dec', and an optional 'source' column tagging origin.
        - Labels from different rows that fall within `cluster_arcsec` (arcseconds)
            will be merged. Using angular separation avoids accidental merges when
            RA/Dec differences are tiny but project to different pixel separations.
        """
        rgb = np.asarray(rgb, dtype=np.float32)
        if rgb.ndim == 2:
            rgb = np.stack([rgb] * 3, axis=-1)
        elif rgb.shape[2] != 3:
            raise ValueError("RGB image must have 3 channels")

        ax.imshow(rgb, origin="lower")

        ny, nx = rgb.shape[:2]

        # build entries (factored into helper)
        entries = self._build_entries(results, nx, ny)

        # cluster entries (factored into helper)
        groups = self._cluster_entries(entries, cluster_arcsec, cluster_pixels)
                
        texts = []
        fig = ax.figure
        # ensure renderer is available for text extent calculations
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

        for g in groups:
            cx, cy = g["x"], g["y"]

            # draw a single marker at the centroid
            ax.plot(
                cx,
                cy,
                marker="o",
                markersize=5,
                markerfacecolor="none",
                markeredgecolor="green",
                markeredgewidth=0.6,
            )

            # Build a single combined label for the group (helper)
            full_label = self._compose_group_label(g)
            if not full_label:
                continue

            # place single text at centroid (offset slightly upward)
            data_pos = (cx, cy)
            # convert to display coords to add pixel offset then back to data coords
            disp = ax.transData.transform(data_pos)
            disp = (disp[0], disp[1] + y_offset_pixels)
            data_pos = ax.transData.inverted().transform(disp)

            txt = ax.text(
                data_pos[0],
                data_pos[1],
                full_label,
                color="#ffffff",
                fontsize=10,
                fontfamily="sans-serif",
                fontweight="bold",
                ha="center",
                va="bottom",
                transform=ax.transData,
            )
            txt.set_path_effects([
                path_effects.Stroke(linewidth=2, foreground="black"),
                path_effects.Normal(),
            ])
            texts.append(txt)

        ax.axis("off")
        return texts

