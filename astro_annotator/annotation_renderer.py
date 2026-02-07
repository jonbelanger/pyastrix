# annotation_renderer.py
import numpy as np
from astropy.coordinates import SkyCoord
import matplotlib


class AnnotationRenderer:
    """Render catalog annotations on a FITS/RGB image using a provided Axes.

    Changes made:
      - Groups nearby detections (pixel-based clustering) and merges labels
      - Supports multiple label parts per object (e.g. name, magnitude, distance)
      - Renders each label part in a source-specific color (hard-coded mapping)
    """

    def __init__(self, wcs):
        self.wcs = wcs

        # simple hard-coded color map by catalog/source keyword
        self.color_map = {
            "gaia": "#4D9600",
            "simbad": "#1f77b4",
            "skybot": "#ff7f0e",
        }

    def _extract_parts(self, row, source):
        """Return list of (text, color) for a row depending on source fields."""
        parts = []
        src_key = (source or "").lower()

        color = self.color_map.get(src_key, "#ffffff")

        # Gaia-style fields
        if "gaia" in src_key:
            mag = row.get("phot_g_mean_mag", None)
            if mag is not None and not np.isnan(mag):
                parts.append((f"G={mag:.2f}", color))

            distance = row.get("distance_ly", None) or row.get("distance", None)
            if distance is not None and not np.isnan(distance):
                parts.append((f"{float(distance):.0f} ly", color))

            # prefer a textual source id if it contains letters
            sid = str(row.get("source_id", ""))
            if sid and any(c.isalpha() for c in sid):
                parts.insert(0, (sid, color))

        # SIMBAD-like catalogs often provide a 'main_id' or 'name'
        elif "simbad" in src_key:
            name = row.get("main_id", None) or row.get("name", None) or row.get("id", None)
            if name:
                parts.append((str(name), color))

        else:
            # generic fallback: include any textual id, magnitude, distance
            sid = row.get("source_id", None) or row.get("id", None)
            if sid:
                sid_s = str(sid)
                if any(c.isalpha() for c in sid_s):
                    parts.append((sid_s, color))

            mag = row.get("phot_g_mean_mag", None)
            if mag is not None and not np.isnan(mag):
                parts.append((f"G={mag:.2f}", color))

        return parts

    def render(self, ax, results, rgb, cluster_pixels=8, y_offset_pixels=10, spacing_pixels=4):
        """Render combined annotations.

        - `results` should be an Astropy Table or iterable with fields
          including 'ra', 'dec', and an optional 'source' column tagging origin.
        - Labels from different rows that fall within `cluster_pixels` will be merged.
        """
        rgb = np.asarray(rgb, dtype=np.float32)
        if rgb.ndim == 2:
            rgb = np.stack([rgb] * 3, axis=-1)
        elif rgb.shape[2] != 3:
            raise ValueError("RGB image must have 3 channels")

        ax.imshow(rgb, origin="lower")

        ny, nx = rgb.shape[:2]

        # build entries with pixel coords and parts
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
            parts = self._extract_parts(row, src)
            if not parts:
                continue

            entries.append({"x": float(x_pix), "y": float(y_pix), "parts": parts})

        # simple single-link clustering in pixel space
        groups = []
        for e in entries:
            placed = False
            for g in groups:
                dx = e["x"] - g["x"]
                dy = e["y"] - g["y"]
                if np.hypot(dx, dy) <= cluster_pixels:
                    g["members"].append(e)
                    # update centroid
                    n = len(g["members"])
                    g["x"] = (g["x"] * (n - 1) + e["x"]) / n
                    g["y"] = (g["y"] * (n - 1) + e["y"]) / n
                    placed = True
                    break
            if not placed:
                groups.append({"x": e["x"], "y": e["y"], "members": [e]})

        texts = []
        fig = ax.figure
        # ensure renderer is available for text extent calculations
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

        import matplotlib.patheffects as path_effects

        for g in groups:
            cx, cy = g["x"], g["y"]

            # draw a single marker at the centroid
            ax.plot(
                cx,
                cy,
                marker="o",
                markersize=5,
                markerfacecolor="none",
                markeredgecolor="white",
                markeredgewidth=0.6,
            )

            # collect parts, preserving insertion order and avoiding duplicates
            seen = set()
            parts = []
            for m in g["members"]:
                for txt, col in m["parts"]:
                    key = (txt, col)
                    if key in seen:
                        continue
                    seen.add(key)
                    parts.append((txt, col))

            if not parts:
                continue

            # start position in display coords (slightly above point)
            start_disp = ax.transData.transform((cx, cy + 0))
            start_disp[1] += y_offset_pixels

            cur_x = start_disp[0]
            y_disp = start_disp[1]

            for text_str, color in parts:
                # convert display back to data coords for placing text
                data_pos = ax.transData.inverted().transform((cur_x, y_disp))
                txt = ax.text(
                    data_pos[0],
                    data_pos[1],
                    text_str,
                    color=color,
                    fontsize=10,
                    fontfamily="sans-serif",
                    fontweight="bold",
                    ha="left",
                    va="bottom",
                    transform=ax.transData,
                )
                txt.set_path_effects([
                    path_effects.Stroke(linewidth=2, foreground="black"),
                    path_effects.Normal(),
                ])

                texts.append(txt)

                # measure width and advance
                bb = txt.get_window_extent(renderer=renderer)
                cur_x += bb.width + spacing_pixels

        ax.axis("off")
        return texts

