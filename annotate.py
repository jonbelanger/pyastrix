import sys
from pathlib import Path
import argparse

MODULE_DIR = Path(__file__).resolve().parent / "astro_annotator"
sys.path.insert(0, str(MODULE_DIR.parent))

from astro_annotator import (
    FitsMeta,
    GaiaCatalog,
    ImageNormalizer,
    AnnotationRenderer,
    ImageWriter,
    SimbadCatalog
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
    parser.add_argument("--skybot", action="store_true", help="Enable SkyBot annotations")

    # Gaia-specific options
    parser.add_argument("--gaia-mag-limit", type=float, default=10)
    parser.add_argument("--gaia-parallax-snr", type=float, default=5)

    return parser.parse_args()


def main():
    args = parse_args()

    # ---- FITS ----
    fits_meta = FitsMeta(args.fits)

    # ---- normalize image ----
    normalizer = ImageNormalizer()
    rgb = normalizer.normalize(fits_meta.data)

    # ---- renderer ----
    renderer = AnnotationRenderer(fits_meta.wcs)

    # ---- collect enabled catalogs ----
    catalogs = []

    if args.gaia:
        print("doing gaia")
        catalogs.append(
            GaiaCatalog(
                mag_limit=args.gaia_mag_limit,
                parallax_snr=args.gaia_parallax_snr
            )
        )

    if args.simbad:
        catalogs.append(SimbadCatalog())

    if args.skybot:
        raise NotImplementedError("SkyBot not wired yet")

    # ---- writer ----
    writer = ImageWriter()

    # ---- render + save ----
    fig, ax = writer.begin(rgb)

    # aggregate results from all catalogs into a single table so the
    # renderer can merge nearby sources and append multiple labels
    from astropy.table import vstack
    all_tables = []

    for catalog in catalogs:
        results = catalog.query(
            center=fits_meta.center,
            radius=fits_meta.radius
        )

        # tag rows with source name for color mapping and labeling
        src_name = getattr(catalog, "name", None) or catalog.__class__.__name__.lower()
        if "source" not in results.colnames:
            results["source"] = [src_name] * len(results)

        all_tables.append(results)

    texts = []
    if all_tables:
        combined = vstack(all_tables)
        texts = renderer.render(ax, combined, rgb)

    writer.finish(fig, ax, args.output, texts)


if __name__ == "__main__":
    main()