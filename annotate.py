import sys
from pathlib import Path
import argparse
from astropy.table import vstack

MODULE_DIR = Path(__file__).resolve().parent / "pyastrix"
sys.path.insert(0, str(MODULE_DIR.parent))

from pyastrix import (
    FitsMeta,
    GaiaCatalog,
    ImageNormalizer,
    AnnotationRenderer,
    ImageWriter,
    SimbadCatalog,
    AsteroidAnnotator
)


def parse_args():
    parser = argparse.ArgumentParser(
        description="Annotate a FITS image using astronomical catalogs"
    )

    # Exactly one input FITS file
    parser.add_argument(
        "fits",
        help="Input FITS file"
    )

    # One output image
    parser.add_argument(
        "-o", "--output",
        default="annotated.png",
        help="Output image file"
    )

    # Annotation toggles
    parser.add_argument("--gaia", action="store_true", help="Enable Gaia annotations")
    parser.add_argument("--simbad", action="store_true", help="Enable SIMBAD annotations")
    parser.add_argument("--asteroids", action="store_true", help="Enable asteroid annotations (V < 20 mag)")

    # Gaia-specific options
    parser.add_argument("--gaia-mag-limit", type=float, default=10)
    parser.add_argument("--gaia-parallax-snr", type=float, default=5)
    parser.add_argument("--cluster-arcsec", type=float, default=.1,
                        help="Angular clustering radius in arcseconds for merging labels")
    parser.add_argument("--cluster-pixels", type=float, default=1,
                        help="Pixel clustering radius as fallback for merging labels")

    return parser.parse_args()


def main():
    args = parse_args()

    fits_meta = FitsMeta(args.fits)

    normalizer = ImageNormalizer()
    rgb = normalizer.normalize(fits_meta.data)

    if hasattr(fits_meta, 'no_valid_wcs') and fits_meta.no_valid_wcs:
        print("[ERROR] No valid WCS present. Skipping annotation rendering.")
        return
    renderer = AnnotationRenderer(fits_meta.wcs)

    catalogs = []

    if args.gaia:
        catalogs.append(
            GaiaCatalog(
                mag_limit=args.gaia_mag_limit,
                parallax_snr=args.gaia_parallax_snr
            )
        )

    if args.simbad:
        catalogs.append(SimbadCatalog())

    if args.asteroids:
        catalogs.append(AsteroidAnnotator())

    writer = ImageWriter()

    fig, ax = writer.begin(rgb)

    all_tables = []

    for catalog in catalogs:
        # Pass observation time for asteroid queries
        query_kwargs = {
            "center": fits_meta.center,
            "radius": fits_meta.radius
        }
        if isinstance(catalog, AsteroidAnnotator):
            query_kwargs["obs_time"] = fits_meta.obs_time

        cat_name = getattr(catalog, "name", None) or catalog.__class__.__name__.lower()
        print(f"[DEBUG] Querying catalog: {cat_name}")
        results = catalog.query(**query_kwargs)
        print(f"[DEBUG] {cat_name} returned {len(results) if results is not None else 0} results.")

        # tag rows with source name for color mapping and labeling
        if results is not None and "source" not in results.colnames:
            results["source"] = [cat_name] * len(results)

        # ensure textual IDs are strings so vstack() doesn't fail on mixed dtypes
        if results is not None and "source_id" in results.colnames:
            try:
                results["source_id"] = [str(x) for x in results["source_id"]]
            except Exception:
                # fallback: create stringified column
                results.remove_column("source_id")
                results["source_id"] = [str(x) for x in results.iterator()]

        if results is not None:
            all_tables.append(results)

    texts = []
    if all_tables:
        combined = vstack(all_tables)
        texts = renderer.render(
            ax,
            combined,
            rgb,
            cluster_arcsec=args.cluster_arcsec,
            cluster_pixels=args.cluster_pixels,
        )

    writer.finish(fig, ax, args.output, texts)


if __name__ == "__main__":
    main()