#!/usr/bin/env python3
"""
rebuild_combined_csv.py
------------------------
One-time migration + combined.csv rebuild for the 2026-07 field data drop:

  1. Converts the 3 newly added GeoJSON exports (LaElenita, MindoMirador,
     Saloya) through the existing geojson-to-csv.py -> convert_observations.py
     pipeline, producing data/csv/Processed<site>.csv for each.
  2. The legacy MiradorMindo GeoJSON export bundled observations from all
     three of those locations under one site name. Splits
     data/csv/ProcessedMiradorMindo.csv by nearest-centroid proximity to each
     new site's own GeoJSON points and appends the reassigned rows onto the
     matching Processed<site>.csv from step 1.
  3. Rebuilds data/csv/combined.csv as the union of every current
     Processed<site>.csv (Maquipucuna, Mashpi, MindoTarabita, LaElenita,
     MindoMirador, Saloya, Yanayacu), matching the file's existing formatting
     conventions (NA for missing values, ISO-8601 datetimes, R-style numeric
     formatting for the double-typed columns per COMBINED_COL_TYPES in
     helper_functions.R).

The original ProcessedMiradorMindo.csv, geojson_to_csv/raw/MiradorMindo.geojson
and geojson_to_csv/csv/MiradorMindo.csv are left untouched (superseded, not
deleted) as an audit trail.

combinedv2.csv / combinedv3.csv (manually curated species IDs) are NOT
touched -- promote the new/reassigned rows from combined.csv into those
yourself, same as every previous field data drop.

Safe to re-run: steps 1 and 2 always regenerate Processed<site>.csv fresh
from the source GeoJSON before appending, so re-running does not double up
the migrated rows.

Run from repo root:
    python scripts/data_prep/rebuild_combined_csv.py
"""
import csv
import json
import subprocess
import sys
from datetime import datetime
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent.parent
RAW_DIR = ROOT / "geojson_to_csv" / "raw"
CSV_DIR = ROOT / "geojson_to_csv" / "csv"
DATA_CSV_DIR = ROOT / "data" / "csv"
GEOJSON_TO_CSV = ROOT / "scripts" / "geojson-csv-sql-conversion-tools" / "python" / "geojson-to-csv.py"
CONVERT_OBS = ROOT / "scripts" / "data_prep" / "convert_observations.py"

NEW_SITES = ["LaElenita", "MindoMirador", "Saloya"]
LEGACY_SITE = "MiradorMindo"
# Final site order for combined.csv (alphabetical -- matches the order the
# repo's various SITES=(...) arrays already use for the other 4 sites).
ALL_SITES = ["LaElenita", "Maquipucuna", "Mashpi", "MindoMirador", "MindoTarabita", "Saloya", "Yanayacu"]

OUTPUT_COLUMNS = [
    'Source', 'Area_or_Site', 'lat', 'lon', 'Elevation_m',
    'FieldID', 'Abundance', 'Height_m', 'CanopyHeight_m',
    'pictures', 'note', 'AI_ID', 'FinalID', 'Genus', 'species', 'datetime'
]
DOUBLE_COLS = {'lat', 'lon', 'Abundance', 'Height_m', 'CanopyHeight_m'}  # COMBINED_COL_TYPES col_double()
# ~2.2km -- bigger than any real within-site GPS scatter seen in these exports
# (<0.006 deg), smaller than the real gaps between the 3 sites (>0.05 deg).
OUTLIER_DEGREES = 0.02


def run_conversion(site_name):
    geojson_path = RAW_DIR / f"{site_name}.geojson"
    raw_csv = CSV_DIR / f"{site_name}.csv"
    out_csv = DATA_CSV_DIR / f"Processed{site_name}.csv"

    print(f"[{site_name}] GeoJSON -> CSV")
    subprocess.run([sys.executable, str(GEOJSON_TO_CSV),
                     "--input", str(geojson_path), "--output", str(raw_csv)], check=True)

    print(f"[{site_name}] CSV -> Processed{site_name}.csv")
    subprocess.run([sys.executable, str(CONVERT_OBS),
                     "--input", str(raw_csv), "--output", str(out_csv),
                     "--site", site_name], check=True)
    return out_csv


def load_geojson_points(site_name):
    d = json.loads((RAW_DIR / f"{site_name}.geojson").read_text())
    return [tuple(ft["geometry"]["coordinates"][:2])
            for ft in d["features"] if ft["geometry"]["type"] == "Point"]  # (lon, lat)


def robust_centroid(points):
    """Mean of points, after dropping any farther than OUTLIER_DEGREES from
    the median -- guards against a single mislabeled/stray GPS point (seen in
    the Saloya export) skewing the site's centroid."""
    lons = sorted(p[0] for p in points)
    lats = sorted(p[1] for p in points)
    med = (lons[len(lons) // 2], lats[len(lats) // 2])
    kept = [p for p in points
            if abs(p[0] - med[0]) <= OUTLIER_DEGREES and abs(p[1] - med[1]) <= OUTLIER_DEGREES]
    if not kept:
        kept = points
    return (sum(p[0] for p in kept) / len(kept), sum(p[1] for p in kept) / len(kept))


def split_legacy_rows():
    legacy_csv = DATA_CSV_DIR / f"Processed{LEGACY_SITE}.csv"
    with legacy_csv.open(newline="") as f:
        rows = list(csv.DictReader(f))

    centroids = {s: robust_centroid(load_geojson_points(s)) for s in NEW_SITES}

    assigned = {s: [] for s in NEW_SITES}
    for row in rows:
        lon, lat = float(row["lon"]), float(row["lat"])
        best = min(NEW_SITES, key=lambda s: (lon - centroids[s][0]) ** 2 + (lat - centroids[s][1]) ** 2)
        row["Area_or_Site"] = best
        assigned[best].append(row)
        print(f"  legacy row lon={lon:.6f} lat={lat:.6f} -> {best}")
    return assigned


def append_legacy_rows(site_name, out_csv, legacy_rows):
    if not legacy_rows:
        return
    with out_csv.open(newline="") as f:
        rows = list(csv.DictReader(f))
    rows.extend(legacy_rows)
    with out_csv.open("w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=OUTPUT_COLUMNS)
        w.writeheader()
        for r in rows:
            w.writerow({c: r.get(c, "") for c in OUTPUT_COLUMNS})
    print(f"[{site_name}] +{len(legacy_rows)} rows moved from {LEGACY_SITE} -> {out_csv} ({len(rows)} total)")


def r_double_str(raw):
    if raw is None or raw == "":
        return "NA"
    f = float(raw)
    return str(int(f)) if f == int(f) else str(f)


def r_string_str(raw):
    return "NA" if raw is None or raw == "" else raw


def r_datetime_str(raw):
    if raw is None or raw == "":
        return "NA"
    dt = datetime.strptime(raw, "%Y-%m-%d %H:%M:%S")
    return dt.strftime("%Y-%m-%dT%H:%M:%SZ")


def rebuild_combined():
    out_path = DATA_CSV_DIR / "combined.csv"
    total = 0
    with out_path.open("w", newline="") as out_f:
        w = csv.writer(out_f)
        w.writerow(OUTPUT_COLUMNS)
        for site in ALL_SITES:
            site_csv = DATA_CSV_DIR / f"Processed{site}.csv"
            with site_csv.open(newline="") as f:
                for row in csv.DictReader(f):
                    out_row = []
                    for col in OUTPUT_COLUMNS:
                        raw = row.get(col, "")
                        if col == "datetime":
                            out_row.append(r_datetime_str(raw))
                        elif col in DOUBLE_COLS:
                            out_row.append(r_double_str(raw))
                        else:
                            out_row.append(r_string_str(raw))
                    w.writerow(out_row)
                    total += 1
    print(f"combined.csv rebuilt: {total} rows -> {out_path}")


def main():
    for site in NEW_SITES:
        run_conversion(site)

    print(f"\nSplitting legacy {LEGACY_SITE} observations by proximity...")
    assigned = split_legacy_rows()
    for site in NEW_SITES:
        append_legacy_rows(site, DATA_CSV_DIR / f"Processed{site}.csv", assigned[site])

    print()
    rebuild_combined()


if __name__ == "__main__":
    main()
