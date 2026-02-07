import time
import os
from pathlib import Path
import subprocess
import re

fits_path = Path(r"C:\Users\jon\Desktop\M42_half_hour_2.fit")

aladin_exe = r"C:\Program Files\Aladin\aladin.exe"

output_fits_dir = Path(r"C:\Users\Jon\astro_automation\snapshots")      # Where snapshots go
output_fits_dir.mkdir(exist_ok=True)

ajs_path = Path(r".\temp_live.ajs")

#tap_url = "https://tapvizier.cds.unistra.fr/TAPVizieR/tap"
tap_url = "cds.vizier/tap"


tap_queries = [
    ("Bright_Object_Stars", '''SELECT TOP 9999
  DR3Name,
  RA_ICRS,
  DE_ICRS,
  Plx,
  e_Plx,
  RPlx,
  Gmag,
  3261.56 / Plx AS distance_ly
FROM \\"I/355/gaiadr3\\"
WHERE
  CONTAINS(
    POINT('ICRS', RA_ICRS, DE_ICRS),
    CIRCLE('ICRS', 83.82344, -5.36601, 0.72733)
  ) = 1
AND Plx BETWEEN 2.2 AND 2.7
AND RPlx > 5
AND Gmag <= 8'''),
    ("Object_Stars", '''SELECT TOP 9999
  DR3Name,
  RA_ICRS,
  DE_ICRS,
  Plx,
  e_Plx,
  RPlx,
  Gmag,
  3261.56 / Plx AS distance_ly
FROM \\"I/355/gaiadr3\\"
WHERE
  CONTAINS(
    POINT('ICRS', RA_ICRS, DE_ICRS),
    CIRCLE('ICRS', 83.82344, -5.36601, 0.72733)
  ) = 1
AND Plx BETWEEN 2.2 AND 2.7
AND RPlx > 5''')
]

last_mtime = 0

while True:
    try:
        mtime = fits_path.stat().st_mtime
        if mtime != last_mtime:
            last_mtime = mtime
            print(f"[INFO] FITS updated, generating .ajs script...")

            # Generate output FITS filename
            timestamp = int(time.time())
            output_fits = output_fits_dir / f"orion_snapshot_{timestamp}.png"

            # --- Write Aladin script ---
            with open(ajs_path, "w") as f:
                # Load FITS
                f.write(f'load "file:///{fits_path.as_posix()}"\n')

                # Load TAP queries
                layer_idx = 2
                for name, query in tap_queries:
                    query_clean = re.sub(r'\s+',' ', " ".join(query.splitlines())).strip()
                    f.write(f'{name}=get TAP({tap_url},"{query_clean}")\n')
                    layer_idx += 1

                # Save current view as FITS
                f.write(f'save "{output_fits.as_posix()}"\n')

            # --- Run Aladin script ---
            subprocess.run([aladin_exe, "-script", str(ajs_path)])
            print(f"[INFO] Saved snapshot to {output_fits}")

        time.sleep(5)  # check every 5 seconds

    except KeyboardInterrupt:
        print("Watcher stopped.")
        break