#!/usr/bin/env python3
"""Copy the editable Figure 1 artwork and render its publication PNG."""
from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess

ROOT = Path(__file__).resolve().parents[3]
SOURCE = ROOT / "figures/released/fig1/recomb_figure1_editable_2.svg"
OUT = Path(os.environ["RECOMB_OUTPUT_DIR"]).resolve()
if not OUT.is_relative_to(ROOT):
    raise ValueError("Figure 1 outputs must be written inside the publication repository.")
OUT.mkdir(parents=True, exist_ok=True)
if not SOURCE.is_file() or SOURCE.is_symlink():
    raise FileNotFoundError(f"Editable Figure 1 SVG is unavailable: {SOURCE}")
inkscape = shutil.which("inkscape")
if not inkscape:
    raise RuntimeError("Inkscape is required to render Figure 1 from the editable SVG.")

svg_output = OUT / SOURCE.name
png_output = OUT / "recomb_figure_1.png"
shutil.copyfile(SOURCE, svg_output)
subprocess.run(
    [
        inkscape, str(SOURCE), "--export-type=png", f"--export-filename={png_output}",
        "--export-width=2400", "--export-background=#ffffff", "--export-background-opacity=255",
    ],
    check=True,
)
