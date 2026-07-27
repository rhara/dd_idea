"""Pandas DataFrame construction for the Streamlit Overview / Active-site
comparison tabs (`view.py`)."""
from __future__ import annotations

from typing import Tuple

import pandas as pd


def overview_dataframe(report: dict) -> pd.DataFrame:
    """`"% identity to reference"` is `None` (not the string `"(reference)"`)
    for the reference row itself -- mixing a string into an otherwise-float
    column makes pandas fall back to an `object` dtype, which Streamlit's
    Arrow serialization can't always convert cleanly (surfaced as a "Expected
    bytes, got a 'float' object" error, silently "fixed" by re-running the
    whole script -- including rebuilding the Structure overlay tab's py3Dmol
    scene from scratch, which defeats that scene's own resize-after-visible
    fix in `viewer3d`). `"(reference)"` is folded into `Note` instead, kept
    as a plain string column throughout."""
    reference = report["reference"]
    rows = []
    for p in report["proteins"]:
        is_reference = p["accession"] == reference
        note = "(reference)" if is_reference else (p["align_error"] or "")
        rows.append(
            {
                "Accession": p["accession"],
                "Name": p["name"],
                "Length": p["length"],
                "% identity to reference": None if is_reference else round(p["pct_identity"], 1),
                "RMSD to reference (Å)": round(p["rmsd"], 3) if p["rmsd"] is not None else None,
                "Note": note,
            }
        )
    return pd.DataFrame(rows)


def pocket_comparison_frames(report: dict) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """Two same-shape DataFrames (one row per protein, one column per
    reference pocket residue -- proteins as rows rather than columns reads
    each protein's own mapped residues left-to-right like a short sequence,
    instead of top-to-bottom): `values` (display string, e.g. `"Y81"` or
    `"-"` for no counterpart) and `conservation` (the raw category, used to
    color `values`'s cells in the Streamlit app -- see `view.py`'s
    `CONSERVATION_COLORS`)."""
    proteins = report["proteins"]
    col_labels = [f"{c['reference_residue']}{c['reference_position']}" for c in proteins[0]["pocket_comparison"]]

    row_labels = []
    values_rows = []
    conservation_rows = []
    for p in proteins:
        row_labels.append(f"{p['accession']} ({p['name']})")
        values_rows.append([
            f"{c['target_residue']}{c['target_position']}" if c["target_residue"] else "-"
            for c in p["pocket_comparison"]
        ])
        conservation_rows.append([c["conservation"] for c in p["pocket_comparison"]])

    return (
        pd.DataFrame(values_rows, index=row_labels, columns=col_labels),
        pd.DataFrame(conservation_rows, index=row_labels, columns=col_labels),
    )


def candidates_dataframe(candidates_report: dict) -> pd.DataFrame:
    rows = [
        {"Accession": c["accession"], "Name": c["name"], "Length": c["length"], "% identity to seed": round(c["pct_identity"], 1)}
        for c in candidates_report["candidates"]
    ]
    return pd.DataFrame(rows)
