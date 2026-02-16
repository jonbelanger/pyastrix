# annotation_renderer.py
import numpy as np
from astropy.coordinates import SkyCoord
import astropy.units as u
import matplotlib
import matplotlib.patheffects as path_effects
from matplotlib.patches import Ellipse

class AnnotationRenderer:
    """Render catalog annotations on a FITS/RGB image using a provided Axes.

    Changes made:
      - Groups nearby detections (pixel-based clustering) and merges labels
      - Supports multiple label parts per object (e.g. name, magnitude, distance)
      - Renders each label part in a source-specific color (hard-coded mapping)
    """

    def __init__(self, wcs):
        self.wcs = wcs
        if self.wcs is None:
            print("[ERROR] No valid WCS present. AnnotationRenderer will not function.")
            return
        print("[DEBUG] WCS info:")
        print(f"  WCS object: {self.wcs}")
        try:
            print(f"  WCS ctype: {self.wcs.wcs.ctype}")
            print(f"  WCS naxis: {self.wcs.wcs.naxis}")
        except Exception as e:
            print(f"  [WARNING] Could not read WCS ctype/naxis: {e}")

        # colorization removed: labels will be rendered as single combined strings

    def _build_entries(self, results, nx, ny):
        """Build normalized entries with coords, pixel positions and a single label."""
        entries = []
        
        # Check which size fields exist in the table
        has_majaxis = "galdim_majaxis" in results.colnames
        has_minaxis = "galdim_minaxis" in results.colnames
        has_angle = "galdim_angle" in results.colnames
        has_vel_ra = "velocity_ra" in results.colnames
        has_vel_dec = "velocity_dec" in results.colnames
        
        for row in results:
            try:
                ra = float(row["ra"])
                dec = float(row["dec"])
                coord = SkyCoord(ra, dec, unit="deg", frame="icrs")
            except Exception:
                continue

            if self.wcs is None:
                raise RuntimeError("No valid WCS present. Cannot convert world to pixel coordinates.")
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

            entry = {
                "coord": coord,
                "x": float(x_pix),
                "y": float(y_pix),
                "label": label,
                "source": src,
            }
            
            # Copy size fields if they exist (for extended objects)
            if has_majaxis:
                try:
                    entry["galdim_majaxis"] = row["galdim_majaxis"]
                except Exception:
                    pass
            if has_minaxis:
                try:
                    entry["galdim_minaxis"] = row["galdim_minaxis"]
                except Exception:
                    pass
            if has_angle:
                try:
                    entry["galdim_angle"] = row["galdim_angle"]
                except Exception:
                    pass
            
            # Copy velocity fields if they exist (for asteroids)
            if has_vel_ra:
                try:
                    entry["velocity_ra"] = row["velocity_ra"]
                except Exception:
                    pass
            if has_vel_dec:
                try:
                    entry["velocity_dec"] = row["velocity_dec"]
                except Exception:
                    pass
            
            entries.append(entry)

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

    def _draw_marker(self, ax, cx, cy, group):
        """Draw either a circle for point sources or an ellipse for extended objects."""
        # Check if any member has size data (SIMBAD extended object)
        majaxis = None
        minaxis = None
        angle = None
        
        for member in group["members"]:
            # Check for SIMBAD galaxy dimension fields (in arcminutes)
            mx = member.get("galdim_majaxis", None)
            my = member.get("galdim_minaxis", None)
            a = member.get("galdim_angle", None)
            
            # Try to extract valid numeric values
            try:
                mx_val = None
                my_val = None
                
                # Handle masked values - check if it's a masked constant
                if mx is not None and str(mx) != "--":
                    if hasattr(mx, 'mask'):
                        if not mx.mask:  # if not masked
                            mx_val = float(mx)
                    else:
                        mx_val = float(mx)
                
                if my is not None and str(my) != "--":
                    if hasattr(my, 'mask'):
                        if not my.mask:  # if not masked
                            my_val = float(my)
                    else:
                        my_val = float(my)
                
                # Check if both are present, positive, and not NaN
                if mx_val is not None and my_val is not None:
                    if not np.isnan(mx_val) and not np.isnan(my_val):
                        if mx_val > 0 and my_val > 0:
                            majaxis = mx_val
                            minaxis = my_val
                            angle = float(a) if a is not None else 0
                            break  # Found valid size data
            except (TypeError, ValueError):
                pass
        
        # If we have size data, draw an ellipse; otherwise draw a circle
        if majaxis is not None and minaxis is not None:
            try:
                # Get pixel scale from WCS (returns Quantity objects)
                pixel_scales = self.wcs.proj_plane_pixel_scales()
                #print(f"DEBUG pixel_scales raw: {pixel_scales}, type: {type(pixel_scales[0])}")
                
                # Try to get in arcsec/pixel
                pix_scale = pixel_scales[0].to(u.arcsec)
                pix_scale_arcsec = float(pix_scale.value)
                #print(f"DEBUG pix_scale_arcsec (from .to(u.arcsec)): {pix_scale_arcsec}")
                
                # Size data is in arcminutes, convert to arcsec then to pixels
                majaxis_arcsec = majaxis * 60.0
                minaxis_arcsec = minaxis * 60.0
                
                width_pix = majaxis_arcsec / pix_scale_arcsec
                height_pix = minaxis_arcsec / pix_scale_arcsec
                
                #print(f"Drawing ellipse: majaxis={majaxis}', minaxis={minaxis}', pix_scale={pix_scale_arcsec:.3f} arcsec/pix, width={width_pix:.1f} pix, height={height_pix:.1f} pix")
                
                ellipse = Ellipse((cx, cy), width=width_pix, height=height_pix, 
                                 angle=angle, facecolor="none", 
                                 edgecolor="green", linewidth=0.6)
                ax.add_patch(ellipse)
            except Exception as e:
                # Fallback to circle if WCS scale extraction fails
                print(f"Ellipse drawing failed: {e}")
                import traceback
                traceback.print_exc()
                ax.plot(cx, cy, marker="o", markersize=5, markerfacecolor="none",
                       markeredgecolor="green", markeredgewidth=0.6)
        else:
            # Point source: draw small circle
            ax.plot(cx, cy, marker="o", markersize=5, markerfacecolor="none",
                   markeredgecolor="green", markeredgewidth=0.6)

    def _draw_velocity_vector(self, ax, cx, cy, group, scale_pixels=50):
        """Draw velocity vector for asteroids."""
        vel_ra = None
        vel_dec = None
        
        for member in group["members"]:
            vra = member.get("velocity_ra", None)
            vdec = member.get("velocity_dec", None)
            
            if vra is not None and vdec is not None:
                try:
                    vel_ra = float(vra)
                    vel_dec = float(vdec)
                    if not (np.isnan(vel_ra) and np.isnan(vel_dec)):
                        break
                except (TypeError, ValueError):
                    pass
        
        if vel_ra is not None and vel_dec is not None:
            print(f"[DEBUG] Drawing velocity vector at ({cx:.1f}, {cy:.1f}) with vel_ra={vel_ra}, vel_dec={vel_dec}")
            # Convert velocity (arcsec/minute) to pixels
            try:
                pixel_scales = self.wcs.proj_plane_pixel_scales()
                pix_scale = pixel_scales[0].to(u.arcsec)
                pix_scale_arcsec = float(pix_scale.value)

                # Convert from arcsec/min to pixels/min
                vel_x_pix = vel_ra / pix_scale_arcsec
                vel_y_pix = vel_dec / pix_scale_arcsec

                # Compute velocity magnitude in pixels/min
                vel_mag_pix = np.sqrt(vel_x_pix**2 + vel_y_pix**2)
                # Set a larger scaling factor for visibility (doubled)
                scaling = 40.0
                # Print only the effective arcsec/pixel after scaling
                # Since velocity is in arcsec/min, print effective arcsec/min per pixel
                effective_arcsec_per_pixel = pix_scale_arcsec / scaling
                print(f"[DEBUG] Effective scale: 1 pixel = {effective_arcsec_per_pixel:.4f} arcsec/min (after scaling)")
                # Optionally, set a minimum visible length
                min_length = 8.0
                # Scale the vector to make the arrow length proportional to velocity
                if vel_mag_pix > 0:
                    scale = max(scaling * vel_mag_pix, min_length) / vel_mag_pix
                else:
                    scale = 0.0
                vel_x_pix *= scale
                vel_y_pix *= scale
                print(f"[DEBUG] Arrow dx={vel_x_pix:.2f} px, dy={vel_y_pix:.2f} px (scaled)")
                # Draw arrow from current position to future position
                ax.arrow(cx, cy, vel_x_pix, vel_y_pix, 
                        head_width=2, head_length=1.5, fc='cyan', ec='cyan', alpha=0.7, linewidth=1)
            except Exception as e:
                print(f"[DEBUG] Exception drawing velocity vector: {e}")
                pass

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

    def _get_brightness_key(self, group):
        """Extract brightness metric from a group for sorting."""
        # Get the brightest magnitude from group members
        min_mag = float('inf')
        for member in group["members"]:
            mag = member.get("phot_g_mean_mag", None) or member.get("magnitude", None)
            if mag is not None and not np.isnan(mag):
                min_mag = min(min_mag, float(mag))
        return min_mag if min_mag != float('inf') else float('inf')

    def _texts_overlap(self, txt1, txt2, fig, margin_pixels=2):
        """Check if two text objects overlap with margin."""
        renderer = fig.canvas.get_renderer()
        
        try:
            bbox1 = txt1.get_window_extent(renderer=renderer)
            bbox2 = txt2.get_window_extent(renderer=renderer)
            
            # Expand bboxes by margin
            bbox1_expanded = bbox1.expanded(1 + margin_pixels/bbox1.width, 1 + margin_pixels/bbox1.height)
            bbox2_expanded = bbox2.expanded(1 + margin_pixels/bbox2.width, 1 + margin_pixels/bbox2.height)
            
            return bbox1_expanded.overlaps(bbox2_expanded)
        except Exception:
            return False

    def _filter_overlapping_labels(self, texts, groups, fig):
        """Remove overlapping labels intelligently, keeping brightest stars.
        
        Uses a heuristic: filter overlapping labels if the magnitude difference is small.
        This preserves the 4 brightest stars in crowded regions while removing nearby dim stars.
        """
        if not texts:
            return texts, groups
        
        brightness_list = []
        for i, (txt, group) in enumerate(zip(texts, groups)):
            brightness = self._get_brightness_key(group)
            brightness_list.append((brightness, txt, i))
        
        brightness_list.sort(key=lambda x: x[0])
        
        kept_indices = []
        kept_texts = []
        kept_brightness = []
        
        for brightness, txt, idx in brightness_list:
            should_filter = False
            
            for kept_txt, kept_bright in zip(kept_texts, kept_brightness):
                if self._texts_overlap(txt, kept_txt, fig, margin_pixels=1):
                    mag_diff = abs(brightness - kept_bright)
                    
                    if mag_diff < 3.0:
                        should_filter = True
                        break
            
            if not should_filter:
                kept_texts.append(txt)
                kept_brightness.append(brightness)
                kept_indices.append(idx)
        
        kept_groups = [groups[i] for i in kept_indices]
        
        return kept_texts, kept_groups

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
        groups_for_texts = []
        fig = ax.figure
        
        # ensure renderer is available for text extent calculations
        fig.canvas.draw()
        renderer = fig.canvas.get_renderer()

        for g in groups:
            cx, cy = g["x"], g["y"]

            # draw a marker (circle for point sources, ellipse for extended objects)
            self._draw_marker(ax, cx, cy, g)
            
            # draw velocity vector if this is an asteroid
            self._draw_velocity_vector(ax, cx, cy, g)

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
                fontsize=7.5,
                fontfamily="sans-serif",
                fontweight="bold",
                ha="center",
                va="bottom",
                transform=ax.transData,
                alpha=0.75,
            )

            txt.set_path_effects([
                path_effects.Stroke(linewidth=1.5, foreground="black"),
                path_effects.Normal(),
            ])

            txt.circle_center = (cx, cy)
            texts.append(txt)
            groups_for_texts.append(g)

        texts, groups_for_texts = self._filter_overlapping_labels(texts, groups_for_texts, fig)
        
        for txt_obj in ax.texts[:]:
            if txt_obj not in texts:
                txt_obj.remove()

        ax.axis("off")
        return texts

