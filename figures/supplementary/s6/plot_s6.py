#!/usr/bin/env python3
"""Render the two-panel S6 figure from publication-frozen Figure 4 geometry."""
from __future__ import annotations

import io
import hashlib
import json
import lzma
import os
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[3]
OUT = Path(os.environ["RECOMB_OUTPUT_DIR"]).resolve()
if not OUT.is_relative_to(ROOT):
    raise ValueError("S6 outputs must be written inside the publication repository.")
OUT.mkdir(parents=True, exist_ok=True)
DATA = ROOT / "data/figure_inputs/fig4"
METHODS = ["iMAT", "GIMME", "CORDA", "RIPTiDe"]
COLORS = {"iMAT": "#4477AA", "GIMME": "#EE6677", "CORDA": "#228833", "RIPTiDe": "#AA3377"}


def load_geometry() -> pd.DataFrame:
    frames = []
    input_records = {
        record["destination_path"]: record
        for record in json.loads((ROOT / "manifests/FIGURE_INPUTS.json").read_text())["inputs"]
    }
    manifest_rel = "data/figure_inputs/fig4/fig4_data_manifest.json"
    required_paths = [manifest_rel] + [
        f"data/figure_inputs/fig4/fig4_geometry_paired/algorithm={method}/part.parquet.xz"
        for method in ["CORDA", "GIMME", "iMAT", "RIPTiDe"]
    ]
    for relative in required_paths:
        record = input_records[relative]
        path = ROOT / relative
        if not record["included"] or record["rights_status"] != "APPROVED":
            raise ValueError(f"S6 input is not approved for release: {relative}")
        if hashlib.sha256(path.read_bytes()).hexdigest() != record["sha256"]:
            raise ValueError(f"S6 frozen input checksum mismatch: {relative}")
    for method in ["CORDA", "GIMME", "iMAT", "RIPTiDe"]:
        path = DATA / f"fig4_geometry_paired/algorithm={method}/part.parquet.xz"
        with lzma.open(path, "rb") as stream:
            frame = pd.read_parquet(io.BytesIO(stream.read()), engine="pyarrow")
        frames.append(frame.assign(algorithm=method))
    paired = pd.concat(frames, ignore_index=True)
    keys = ["algorithm", "evaluation_id", "reaction_id"]
    if len(paired) != 1_671_600 or paired.duplicated(keys).any():
        raise ValueError("S6 requires 1,671,600 unique method/evaluation/reaction keys.")
    if paired.groupby("algorithm").size().to_dict() != {m: 417_900 for m in METHODS}:
        raise ValueError("S6 method populations do not match the frozen 417,900 keys each.")
    required = [
        "A1_direction_explained_magnitude_variance",
        "A2_direction_explained_magnitude_variance",
        "A1_supported_sign_state_count",
        "A2_supported_sign_state_count",
    ]
    missing = sorted(set(required) - set(paired.columns))
    if missing:
        raise ValueError(f"Frozen Figure 4 inputs lack S6 columns: {missing}")
    return paired


def main() -> None:
    paired = load_geometry()
    eta_delta = (
        paired["A2_direction_explained_magnitude_variance"].to_numpy(float)
        - paired["A1_direction_explained_magnitude_variance"].to_numpy(float)
    )
    finite_eta = np.isfinite(
        paired["A1_direction_explained_magnitude_variance"].to_numpy(float)
    ) & np.isfinite(paired["A2_direction_explained_magnitude_variance"].to_numpy(float))
    paired["delta_eta2"] = eta_delta
    eta_counts = paired.loc[finite_eta].groupby("algorithm").size().to_dict()
    expected_eta_counts = {"CORDA": 212_164, "GIMME": 141_365, "iMAT": 207_930, "RIPTiDe": 12_423}
    if eta_counts != expected_eta_counts or sum(eta_counts.values()) != 573_882:
        raise ValueError(f"Unexpected finite paired eta-squared population: {eta_counts}")

    transition_counts = pd.crosstab(
        paired["A1_supported_sign_state_count"].astype(int),
        paired["A2_supported_sign_state_count"].astype(int),
    ).reindex(index=[1, 2, 3], columns=[1, 2, 3], fill_value=0)
    expected_transitions = np.array([[1_293_850, 1_055, 0], [0, 376_470, 0], [0, 0, 225]])
    if not np.array_equal(transition_counts.to_numpy(), expected_transitions):
        raise ValueError("S6 supported-sign-state transition counts differ from frozen reference.")

    fig, axes = plt.subplots(
        1, 2, figsize=(10.2, 4.5), constrained_layout=True,
        gridspec_kw={"width_ratios": [1.35, 1.0]},
    )

    ax = axes[0]
    for method in METHODS:
        values = np.sort(paired.loc[finite_eta & paired.algorithm.eq(method), "delta_eta2"].to_numpy(float))
        y = np.arange(1, len(values) + 1) / len(values)
        ax.plot(values, y, lw=1.7, color=COLORS[method], label=f"{method} (n={len(values):,})")
    ax.axvline(0, color="#444444", lw=1, ls="--")
    ax.set(xlim=(-1.02, 1.02), ylim=(0, 1), xlabel=r"Paired change in $\eta^2$ (A2 − A1)", ylabel="Cumulative fraction")
    ax.set_title(r"A  Direction-explained magnitude variance ($\eta^2$)", loc="left")
    ax.grid(axis="y", color="#dddddd", lw=0.6)
    ax.legend(frameon=False, fontsize=8, loc="lower right")

    ax = axes[1]
    fractions = transition_counts / len(paired)
    im = ax.imshow(
        fractions.to_numpy(), origin="lower", cmap="Blues",
        norm=Normalize(vmin=0, vmax=max(float(fractions.to_numpy().max()), 1e-12)), aspect="equal",
    )
    for i in range(3):
        for j in range(3):
            fraction = float(fractions.iloc[i, j])
            count = int(transition_counts.iloc[i, j])
            ax.text(j, i, f"{fraction:.1%}\n({count:,})", ha="center", va="center", fontsize=8,
                    color="white" if fraction > float(fractions.to_numpy().max()) * 0.55 else "#222222")
    ax.set_xticks([0, 1, 2], ["1", "2", "3"])
    ax.set_yticks([0, 1, 2], ["1", "2", "3"])
    ax.set_xlabel("Number of sign states with\nnonzero probability mass under A2\n" + r"(among $\pi_{A2,-1}$, $\pi_{A2,0}$, $\pi_{A2,+1}$)")
    ax.set_ylabel("Number of sign states with\nnonzero probability mass under A1\n" + r"(among $\pi_{A1,-1}$, $\pi_{A1,0}$, $\pi_{A1,+1}$)")
    ax.set_title(f"B  Supported sign-state transitions\n(n={len(paired):,})", loc="left")
    fig.colorbar(im, ax=ax, shrink=0.78, pad=0.04, label="Fraction of paired rows")

    fig.suptitle("Supplementary Figure S6 | Complementary flux-difference distribution descriptors",
                 fontsize=13, fontweight="bold")
    svg_path = OUT / "supp_fig6.svg"
    fig.savefig(svg_path, bbox_inches="tight")
    # Matplotlib's SVG path serialization pads several path-command lines; trim
    # that irrelevant whitespace so the generated publication SVG passes diff checks.
    svg_path.write_text("\n".join(line.rstrip() for line in svg_path.read_text().splitlines()) + "\n")
    fig.savefig(OUT / "supp_fig6.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()
