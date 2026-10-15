"""Import raw COMSOL cell exports into a strict cell-table dataset (docs/cell_table_contract.md).

    python analysis/import_cell_tables.py --manifest PATH_TO_EXPORT_MANIFEST \
        --output data/processed/cell_tables/DATASET_ID

Raw exports are read, hashed and never modified. Units must match the contract exactly (no silent
conversion). The output directory must not already hold a dataset (historical datasets are never
overwritten). Non-physical data roles (fixtures, analytical examples) are refused under
data/processed/cell_tables/. A qualification report is written next to the dataset; schema validity,
numerical qualification and case match are reported separately, and an import never promotes a
dataset to "qualified" by itself.

Exit codes: 0 imported and schema-valid, 1 import refused or schema-invalid, 2 usage.
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from _common import REPO

from nanofab.device.importer import ImportError_, import_from_manifest


def main(argv=None) -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True, type=Path, help="raw-export manifest (nanofab.raw_export_manifest/1)")
    ap.add_argument("--output", required=True, type=Path, help="new dataset directory")
    args = ap.parse_args(argv)
    if not args.manifest.is_file():
        print(f"import refused: manifest not found: {args.manifest}", file=sys.stderr)
        return 1
    try:
        res = import_from_manifest(args.manifest, args.output, REPO)
    except (ImportError_, FileExistsError, ValueError) as exc:
        print(f"import refused: {exc}", file=sys.stderr)
        return 1
    r = res["report"]
    print(f"dataset {r['dataset_id']} ({r['data_role']}) written: {res['manifest']}")
    print(f"  schema valid:          {r['schema_valid']}")
    print(f"  verification only:     {r['verification_only']}")
    print(f"  numerically qualified: {r['numerically_qualified']}")
    for e in r["errors"]:
        print(f"  ERROR: {e}")
    for q in r["qualification_failures"]:
        print(f"  not qualified: {q}")
    print(f"  qualification report:  {Path(args.output) / 'qualification_report.json'}")
    return 0 if r["schema_valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
