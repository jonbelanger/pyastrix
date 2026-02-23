import time
import math
import shutil
import requests
import subprocess
import argparse
from pathlib import Path
from datetime import datetime
import sys
import json
import os
import win32com.client
import pythoncom
from watchdog.observers import Observer
from watchdog.events import FileSystemEventHandler
from astropy.io import fits
from astropy.table import vstack
import yaml

# Ensure local package import works when running from the repo root
MODULE_DIR = Path(__file__).resolve().parent
if str(MODULE_DIR) not in sys.path:
    sys.path.insert(0, str(MODULE_DIR))

try:
    from pyastrix import (
        FitsMeta,
        GaiaCatalog,
        ImageNormalizer,
        AnnotationRenderer,
        ImageWriter,
        SimbadCatalog,
        AsteroidAnnotator,
    )
except Exception as e:
    print("Warning: could not import pyastrix annotation modules:", e)

# =====================
# CONFIG PLACEHOLDERS (populated from automation.yaml at runtime)
# =====================
# These are intentionally left unset here; defaults should be supplied via
# the YAML config loaded by `load_config()` and applied in `main()`.
SHARPCAP_BASE = None
OUTPUT_DIR = None
ASTAP_CMD = None
MOVE_THRESHOLD_ARCMIN = None

# When True, annotation commands will be skipped. This may be overridden by
# the CLI `--skip-annotation` flag.
SKIP_ANNOTATION = False

def load_config(path: str):
    """Load YAML config from path and return a dict. Missing file returns {}.

    Expected keys: SHARPCAP_BASE, OUTPUT_DIR, ASTAP_CMD, MOVE_THRESHOLD_ARCMIN
    """
    p = Path(path)
    if not p.exists():
        return {}
    try:
        with open(p, 'r', encoding='utf-8') as fh:
            cfg = yaml.safe_load(fh) or {}
            return cfg
    except Exception as e:
        print(f"Warning: failed to load config {path}: {e}")
        return {}

# =====================
# MOUNT MONITOR
# =====================

class MountMonitor:
    def __init__(self):
        self.scope = None
        self.last_ra = None
        self.last_dec = None
        self.last_pier = None
        self.was_slewing = False
        self._connect()

    def _connect(self):
        """Connect or reconnect to ASCOM scope."""
        try:
            self.scope = win32com.client.Dispatch("ASCOM.DeviceHub.Telescope")
            self.scope.Connected = True
            print("ASCOM scope connected")
        except Exception as e:
            print(f"Failed to connect to ASCOM scope: {e}")
            self.scope = None

    def arcmin_delta(self, ra1, dec1, ra2, dec2):
        dra = (ra2 - ra1) * 15 * 60  # hours → arcmin
        ddec = (dec2 - dec1) * 60
        return math.sqrt(dra*dra + ddec*ddec)

    def update(self):
        # Initialize COM for this thread and create a local ASCOM telescope object.
        # This avoids using a COM object created on a different thread (watchdog thread).
        scope = None
        pythoncom.CoInitialize()
        try:
            try:
                scope = win32com.client.Dispatch("ASCOM.DeviceHub.Telescope")
                scope.Connected = True
                ra = scope.RightAscension
                dec = scope.Declination
                slewing = scope.Slewing
                pier = scope.SideOfPier
            except Exception as e:
                print(f"Could not read mount status from ASCOM: {e}")
                # Fall back to last-known coordinates if available; otherwise leave None
                # so the solver won't be passed bogus 0/0 coordinates.
                ra = self.last_ra if self.last_ra is not None else None
                dec = self.last_dec if self.last_dec is not None else None
                slewing = False
                pier = self.last_pier
        finally:
            try:
                if scope is not None:
                    scope.Connected = False
            except Exception:
                pass
            pythoncom.CoUninitialize()

        moved = False

        if self.last_ra is None or ra is None:
            # If we don't have a previous or current coordinate, force a solve
            moved = True
        elif pier != self.last_pier:
            print("Meridian flip detected")
            moved = True
        elif self.was_slewing and not slewing:
            print("Slew finished")
            moved = True
        else:
            # Only compute delta if both previous and current coordinates are present
            if self.last_ra is not None and self.last_dec is not None and ra is not None and dec is not None:
                delta = self.arcmin_delta(self.last_ra, self.last_dec, ra, dec)
                if delta > MOVE_THRESHOLD_ARCMIN:
                    moved = True

        self.last_ra = ra
        self.last_dec = dec
        self.last_pier = pier
        self.was_slewing = slewing

        return moved, ra, dec

# =====================
# WCS CACHE
# =====================

class WCSCache:
    def __init__(self):
        self.header = None

    def update_from_solved(self, solved_fits):
        hdr = fits.getheader(solved_fits)
        wcs_keys = {}
        for key in hdr:
            if key.startswith(("CR", "CD", "CTYPE", "CUNIT", "CDELT", "PC", "PV", "LONPOLE", "LATPOLE")):
                wcs_keys[key] = hdr[key]
        self.header = wcs_keys
        print("WCS cache updated")

    def apply_to(self, input_fits, output_fits):
        if self.header is None:
            print("No WCS cached yet")
            return None
        hdul = fits.open(input_fits)
        hdr = hdul[0].header
        for key, val in self.header.items():
            hdr[key] = val
        hdul.writeto(output_fits, overwrite=True)
        hdul.close()
        return output_fits

# =====================
# SOLVER
# =====================

class Solver:
    def __init__(self, plate_scale=None):
        self.last_wcs = None
        self.session = None
        self.blind = False
        # user-supplied plate scale in arcsec/pixel (None = don't pass)
        self.plate_scale = plate_scale

    def solve(self, fits_path, ra, dec):
        """Run ASTAP to solve the FITS file and return path to WCS-corrected file.

        This uses the `ASTAP_CMD` template defined in the file. The template must
        include a `{fits}` placeholder which will be replaced with the FITS path.
        """
        print(f"Solving with ASTAP: {fits_path.name}")
        OUTPUT_DIR.mkdir(exist_ok=True)
        out_path = OUTPUT_DIR / (fits_path.stem + "_wcs.fits")

        # Wait until the FITS file is stable (fully written) and accessible
        def _wait_for_stable_file(p, timeout=120, stable_secs=2):
            start = time.time()
            last_size = -1
            last_change = time.time()
            while time.time() - start < timeout:
                try:
                    size = os.path.getsize(p)
                except OSError:
                    size = -1
                if size != last_size:
                    last_size = size
                    last_change = time.time()
                else:
                    if time.time() - last_change >= stable_secs:
                        return True
                time.sleep(0.5)
            return False

        fits_abs = str(fits_path.resolve())
        if not _wait_for_stable_file(fits_abs, timeout=120, stable_secs=2):
            print(f"File {fits_abs} not stable/accessible after wait; aborting ASTAP solve.")
            return None

        # If user provided plate scale, write it to FITS header as SCALE keyword
        # ASTAP will read this from the header to avoid blind-scale searches
        if getattr(self, 'plate_scale', None) is not None:
            try:
                ps = float(self.plate_scale)
                hdul = fits.open(fits_abs)
                hdul[0].header['SCALE'] = (ps, "Plate scale in arcsec/pixel (user-supplied)")
                hdul.flush()
                hdul.close()
                print(f"Wrote SCALE keyword to FITS header: {ps:.3f} arcsec/pixel")
            except Exception as e:
                print(f"Warning: could not write plate scale to FITS header: {e}")
        
        # Build command and run ASTAP (minimal logging)
        base_cmd = [str(part).format(fits=fits_abs) for part in ASTAP_CMD]
        cmd = list(base_cmd)
        try:
            if not getattr(self, 'blind', False) and ra is not None and dec is not None:
                ra_str = f"{float(ra):.6f}"
                spd = float(dec) + 90.0
                spd_str = f"{spd:.6f}"
                cmd.extend(["-ra", ra_str, "-spd", spd_str])
        except Exception:
            pass

        solve_start = time.time()
        try:
            proc = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, timeout=300, cwd=str(fits_path.parent))
        except Exception:
            proc = None

        # If ASTAP failed (non-zero return), surface stderr for debugging
        if proc is not None and proc.returncode != 0:
            err = proc.stderr.decode(errors='ignore') if proc.stderr else ''
            print(f"ASTAP failed (rc={proc.returncode}): {err}")

        # Helper to inspect whether a FITS has WCS keys
        def _has_wcs(fpath):
            try:
                hdr = fits.getheader(fpath)
                for k in hdr:
                    if k.startswith(('CR', 'CD', 'CTYPE', 'CDELT')):
                        return True
            except Exception:
                return False
            return False

        # Prefer the input FITS if ASTAP modified it in-place (-update)
        if fits_path.exists() and fits_path.stat().st_mtime > solve_start and _has_wcs(fits_path):
            if fits_path.resolve() != out_path.resolve():
                shutil.copy(fits_path, out_path)
            self.last_wcs = out_path
            return out_path

        # Fallback: common output filenames
        candidates = [
            fits_path.with_name(fits_path.stem + "_wcs.fits"),
            fits_path.with_name(fits_path.stem + ".new.fits"),
            fits_path.with_suffix('.wcs'),
        ]
        for c in candidates:
            if isinstance(c, Path) and c.exists() and _has_wcs(c if c.suffix != '.wcs' else fits_path):
                # If candidate is ASCII .wcs, prefer copying the original FITS if it exists
                if c.suffix == '.wcs' and fits_path.exists() and _has_wcs(fits_path):
                    shutil.copy(fits_path, out_path)
                else:
                    shutil.copy(c, out_path)
                self.last_wcs = out_path
                return out_path

        return None
    

# =====================
# ANNOTATION
# =====================

class AnnotationManager:
    def __init__(self, skip=SKIP_ANNOTATION, **options):
        self.skip = skip
        # Annotation options with sensible defaults
        self.options = {
            'gaia': options.get('gaia', True),
            'simbad': options.get('simbad', True),
            'asteroids': options.get('asteroids', True),
            'gaia_mag_limit': options.get('gaia_mag_limit', 10.0),
            'gaia_parallax_snr': options.get('gaia_parallax_snr', 5.0),
            'cluster_arcsec': options.get('cluster_arcsec', 0.1),
            'cluster_pixels': options.get('cluster_pixels', 1),
            'output_dir': options.get('output_dir', None),
        }

    def run(self, fits_path):
        if self.skip:
            return
        try:
            annotate_native(
                fits_path,
                output_dir=self.options.get('output_dir'),
                gaia=self.options.get('gaia'),
                simbad=self.options.get('simbad'),
                asteroids=self.options.get('asteroids'),
                gaia_mag_limit=self.options.get('gaia_mag_limit'),
                gaia_parallax_snr=self.options.get('gaia_parallax_snr'),
                cluster_arcsec=self.options.get('cluster_arcsec'),
                cluster_pixels=self.options.get('cluster_pixels'),
            )
            return
        except Exception as e:
            print(f"Annotation failed (native): {e}")
            return

def annotate_native(fits_path, output_dir=None, gaia=True, simbad=True, asteroids=True,
                    gaia_mag_limit=10.0, gaia_parallax_snr=5.0, cluster_arcsec=0.1, cluster_pixels=1):
    """Minimal native annotation using pyastrix.

    Writes an annotated PNG to OUTPUT_DIR with suffix `_annotated.png`.
    """
    # Resolve output directory: prefer explicit arg, then global OUTPUT_DIR, else use FITS parent
    if output_dir is None:
        if 'OUTPUT_DIR' in globals() and OUTPUT_DIR is not None:
            output_dir = Path(OUTPUT_DIR)
        else:
            output_dir = Path(fits_path).parent
            print(f"Warning: OUTPUT_DIR not set; falling back to FITS directory: {output_dir}")
    else:
        output_dir = Path(output_dir)

    output_dir.mkdir(parents=True, exist_ok=True)
    out_image = output_dir / (fits_path.stem + "_annotated.png")

    print(f"Native annotation: loading {fits_path}")
    fm = FitsMeta(str(fits_path))

    normalizer = ImageNormalizer()
    rgb = normalizer.normalize(fm.data)

    if getattr(fm, 'no_valid_wcs', False):
        print("[ERROR] No valid WCS present. Skipping annotation rendering.")
        return

    renderer = AnnotationRenderer(fm.wcs)

    catalogs = []
    if gaia:
        catalogs.append(GaiaCatalog(mag_limit=gaia_mag_limit, parallax_snr=gaia_parallax_snr))
    if simbad:
        catalogs.append(SimbadCatalog())
    if asteroids:
        catalogs.append(AsteroidAnnotator())

    writer = ImageWriter()
    fig, ax = writer.begin(rgb)

    all_tables = []
    for catalog in catalogs:
        query_kwargs = {"center": fm.center, "radius": fm.radius}
        if isinstance(catalog, AsteroidAnnotator):
            query_kwargs["obs_time"] = fm.obs_time

        cat_name = getattr(catalog, "name", None) or catalog.__class__.__name__.lower()
        print(f"[DEBUG] Querying catalog: {cat_name}")
        try:
            results = catalog.query(**query_kwargs)
        except Exception as e:
            print(f"Catalog {cat_name} query failed: {e}")
            results = None

        print(f"[DEBUG] {cat_name} returned {len(results) if results is not None else 0} results.")

        if results is not None and "source" not in results.colnames:
            results["source"] = [cat_name] * len(results)

        if results is not None and "source_id" in results.colnames:
            try:
                results["source_id"] = [str(x) for x in results["source_id"]]
            except Exception:
                try:
                    results.remove_column("source_id")
                except Exception:
                    pass
                results["source_id"] = [str(x) for x in results.iterator()]

        if results is not None:
            all_tables.append(results)

    texts = []
    if all_tables:
        combined = vstack(all_tables)
        texts = renderer.render(ax, combined, rgb, cluster_arcsec=cluster_arcsec, cluster_pixels=cluster_pixels)

    writer.finish(fig, ax, str(out_image), texts)
    print(f"Annotated image written: {out_image}")

# =====================
# FILE WATCHER
# =====================

class FitsHandler(FileSystemEventHandler):
    def __init__(self, monitor, solver, wcs_cache, annotator):
        self.monitor = monitor
        self.solver = solver
        self.wcs_cache = wcs_cache
        self.annotator = annotator

    def on_created(self, event):
        if event.is_directory:
            return
        path = Path(event.src_path)
        # Ignore FITS created inside rawframes subdirectories (SharpCap temporary files)
        if any(part.lower() == 'rawframes' for part in path.parts):
            print(f"Skipping rawframes FITS: {path}")
            return
        if path.suffix.lower() != ".fits":
            return
        print(f"New FITS detected: {path}")
        time.sleep(3)  # wait for SharpCap to finish writing

        moved, ra, dec = self.monitor.update()

        if moved or self.solver.last_wcs is None:
            solved_path = self.solver.solve(path, ra, dec)
            if solved_path:
                self.wcs_cache.update_from_solved(solved_path)
                self.annotator.run(solved_path)
        else:
            # apply cached WCS
            output_path = OUTPUT_DIR / (path.stem + "_wcs.fits")
            self.wcs_cache.apply_to(path, output_path)
            self.annotator.run(output_path)
            print("Solve skipped — applied cached WCS")

# =====================
# MAIN
# =====================

class Coordinator:
    """Encapsulates service lifecycle: mount monitor, solver, WCS cache, annotator and watcher."""
    def __init__(self, watch_dir, plate_scale=None, skip_annotation=False, annotation_options=None):
        self.watch_dir = Path(watch_dir)
        self.monitor = MountMonitor()
        self.solver = Solver(plate_scale=plate_scale)
        self.wcs_cache = WCSCache()
        self.annotator = AnnotationManager(skip=skip_annotation, **(annotation_options or {}))
        self.observer = Observer()
        self.handler = FitsHandler(self.monitor, self.solver, self.wcs_cache, self.annotator)

    def start(self):
        self.watch_dir.mkdir(parents=True, exist_ok=True)
        self.observer.schedule(self.handler, str(self.watch_dir), recursive=True)
        self.observer.start()

    def stop(self):
        try:
            self.observer.stop()
        finally:
            try:
                self.observer.join()
            except Exception:
                pass

    def run_forever(self):
        try:
            while True:
                time.sleep(1)
        except KeyboardInterrupt:
            self.stop()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--skip-annotation", action="store_true", help="Skip running annotation commands")
    parser.add_argument("--plate-scale", type=float, help="Override plate scale (arcsec/pixel) to pass to ASTAP")
    parser.add_argument("--config", default="./automation.yaml", help="Path to YAML config file (default: ./automation.yaml)")
    args = parser.parse_args()

    global SKIP_ANNOTATION
    SKIP_ANNOTATION = args.skip_annotation

    # Load YAML config and override defaults if present
    cfg = load_config(args.config)
    if cfg:
        if 'SHARPCAP_BASE' in cfg:
            globals()['SHARPCAP_BASE'] = Path(cfg['SHARPCAP_BASE'])
        if 'OUTPUT_DIR' in cfg:
            globals()['OUTPUT_DIR'] = Path(cfg['OUTPUT_DIR'])
        if 'ASTAP_CMD' in cfg:
            # expect a list
            globals()['ASTAP_CMD'] = cfg['ASTAP_CMD']
        if 'MOVE_THRESHOLD_ARCMIN' in cfg:
            globals()['MOVE_THRESHOLD_ARCMIN'] = cfg['MOVE_THRESHOLD_ARCMIN']
    # Extract annotation options (optional)
    annotation_cfg = cfg.get('ANNOTATION', {}) if cfg else {}

    today = datetime.utcnow().strftime("%Y-%m-%d")
    watch_dir = SHARPCAP_BASE / today
    print(f"Watching {watch_dir}")

    coord = Coordinator(watch_dir, plate_scale=args.plate_scale, skip_annotation=args.skip_annotation, annotation_options=annotation_cfg)
    coord.start()
    try:
        coord.run_forever()
    finally:
        coord.stop()


if __name__ == "__main__":
    main()