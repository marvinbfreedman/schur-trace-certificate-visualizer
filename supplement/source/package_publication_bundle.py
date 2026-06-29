#!/usr/bin/env python3
"""Package the current RH publication bundle into a frozen archive.

This copies the manuscript source, audit outputs, and theorem/certificate
artifacts into a dated directory and then creates a tar.gz archive of that
directory.  The bundle is source-only in this environment because no TeX
toolchain is installed to render the manuscript PDF.
"""

from __future__ import annotations

import json
import shutil
import tarfile
from datetime import datetime
from pathlib import Path


ROOT = Path(__file__).resolve().parent
STAMP = datetime.now().strftime("%Y%m%d-%H%M%S")
BUNDLE_DIR = ROOT / f"publication_bundle_{STAMP}"
ARCHIVE_PATH = ROOT / f"{BUNDLE_DIR.name}.tar.gz"


def should_include(path: Path) -> bool:
    name = path.name
    if path.is_dir():
        return False
    if name in {
        "rh_weyl_positive_draft.tex",
        "rh_klm_notes.md",
        "publication_bundle_manifest.md",
        "publication_bundle_manifest.json",
    }:
        return True
    if name.startswith("publication_audit_") and path.suffix in {".py", ".json", ".md"}:
        return True
    if name.endswith(("_theorem.py", "_theorem.json", "_consequence_theorem.py", "_consequence_theorem.json", "_certificate.py", "_certificate.json")):
        return True
    return False


def main() -> None:
    if BUNDLE_DIR.exists():
        raise SystemExit(f"bundle dir already exists: {BUNDLE_DIR}")
    BUNDLE_DIR.mkdir()

    files: list[str] = []
    for path in sorted(ROOT.iterdir()):
        if should_include(path):
            target = BUNDLE_DIR / path.name
            shutil.copy2(path, target)
            files.append(path.name)

    manifest = {
        "stamp": STAMP,
        "bundle_dir": BUNDLE_DIR.name,
        "archive": ARCHIVE_PATH.name,
        "note": "Source-only bundle; PDF build requires a TeX toolchain not present in this environment.",
        "included_files": files,
    }
    (BUNDLE_DIR / "bundle_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    (ROOT / "publication_bundle_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )

    with tarfile.open(ARCHIVE_PATH, "w:gz") as tf:
        tf.add(BUNDLE_DIR, arcname=BUNDLE_DIR.name)

    print(f"bundle_dir={BUNDLE_DIR}")
    print(f"archive={ARCHIVE_PATH}")
    print(f"file_count={len(files)}")


if __name__ == "__main__":
    main()
