"""py3Dmol multi-protein overlay scene, built from `structalign`'s
superposed AlphaFold-model coordinate files, plus zero or more independent
real-RCSB-structure overlays per protein (see `pdbstruct.py`/`pipeline.py`'s
`pdb_structures` report field). Adapted from
`dd_seqalign.scene.build_overlay_view` for the cross-protein case: that
function highlights the *same* canonical positions across multiple
structures of one protein; here every protein has its own numbering, so
each structure gets its own list of residues to highlight (translated
position-by-position by `pocketmap.py`).
"""
from __future__ import annotations

from pathlib import Path
from typing import Dict, Optional, Sequence, Tuple, Union

import py3Dmol

PALETTE = [
    # No black/gray: an achromatic entry here (the original list had
    # "#7f7f7f", matplotlib's tab10 gray) is hard to tell apart from
    # other structures once rendered semi-transparent, and easy to
    # mistake for "no color at all" -- every entry here should read as a
    # clear, distinct hue at a glance.
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#20c997", "#bcbd22", "#17becf",
]
REFERENCE_COLOR = "#264653"  # dark teal -- distinct from PALETTE but, per the
# achromatic-avoidance rule above, deliberately not gray/black either (the
# original choice, "#444444", was)
SITE_COLOR = "yellow"
LIGAND_COLOR = "magenta"  # fixed, not cycled per-structure like PDB_PALETTE:
# a bound ligand is the actual payoff of the real-structure overlay and
# should read as "the ligand" at a glance regardless of which structure it
# came from -- sharing its parent structure's own cartoon color (the
# original behavior) made it blend in instead of standing out. Distinct
# from SITE_COLOR (used for pocket-residue sticks) so the two don't get
# confused when both are shown at once.

# Distinct hues from PALETTE (protein/AFDB cartoon colors) -- also no
# black/gray, same reasoning as above -- so a real structure's own color
# (used for *both* its solid cartoon and its bound ligand's sticks, see
# assign_pdb_colors) never visually collides with any protein's AFDB
# cartoon, and multiple real structures shown together (even across
# different proteins) stay distinguishable from each other.
PDB_PALETTE = [
    "#e6194b", "#3cb44b", "#ffe119", "#4363d8", "#f58231",
    "#911eb4", "#46f0f0", "#f032e6", "#bcf60c", "#fabebe",
]

# 3Dmol.js's built-in "*Carbon" colorschemes (e.g. "yellowCarbon") tint carbon
# and leave every other element at its RasMol default, but only for a fixed
# set of named CSS colors -- this is the same RasMol table, copied here so
# any per-protein hex can be used as the carbon tint via a
# `{"prop": "elem", "map": {...}}` colorscheme instead of a scheme name.
_HETERO_ELEMENT_COLORS = {
    "H": "#ffffff", "He": "#ffc0cb", "Li": "#b22222", "B": "#00ff00",
    "N": "#8f8fff", "O": "#f00000", "F": "#daa520", "Na": "#0000ff",
    "Mg": "#228b22", "Al": "#808090", "Si": "#daa520", "P": "#ffa500",
    "S": "#ffc832", "Cl": "#00ff00", "Ca": "#808090",
}


def _carbon_tint_scheme(carbon_color: str) -> dict:
    return {"prop": "elem", "map": {**_HETERO_ELEMENT_COLORS, "C": carbon_color}}


def _readable_font_color(background_hex: str) -> str:
    """"black" or "white", whichever reads more clearly on top of
    `background_hex` -- some palette entries (e.g. `PDB_PALETTE`'s bright
    yellow `#ffe119`) are light enough that white label text is nearly
    invisible. Perceived brightness via the standard ITU-R BT.601
    luma weighting (0.299R + 0.587G + 0.114B), same formula commonly used
    for this exact "pick readable text color" problem."""
    h = background_hex.lstrip("#")
    if len(h) != 6:
        return "white"
    r, g, b = (int(h[i:i + 2], 16) for i in (0, 2, 4))
    luma = 0.299 * r + 0.587 * g + 0.114 * b
    return "black" if luma > 150 else "white"


def assign_colors(labels: Sequence[str], reference_label: str) -> Dict[str, str]:
    """One color per label, cycling through `PALETTE`; the reference always
    gets the same fixed `REFERENCE_COLOR` so it reads as "the reference",
    not just another protein in the cycle. Duplicate labels (e.g. a protein's
    AlphaFold entry and its own real-PDB-overlay entries, which
    deliberately share a label so they render in the same color) collapse
    to one color, not one each."""
    colors: Dict[str, str] = {}
    i = 0
    for label in labels:
        if label in colors:
            continue
        if label == reference_label:
            colors[label] = REFERENCE_COLOR
        else:
            colors[label] = PALETTE[i % len(PALETTE)]
            i += 1
    return colors


def assign_pdb_colors(keys: Sequence[str]) -> Dict[str, str]:
    """One color per key (e.g. an `f"{accession}:{pdb_id}"` string),
    cycling through `PDB_PALETTE` -- used for a real structure's cartoon
    *and* its bound ligand's sticks alike, so the two read as one visual
    unit. Assign this over the *full* set of available real-structure
    entries (not just the ones currently checked to show) so a given
    structure's color stays stable as the user toggles others on/off."""
    return {key: PDB_PALETTE[i % len(PDB_PALETTE)] for i, key in enumerate(keys)}


def _ca_coords(pdb_path: str, chain_id: str, resnums: Sequence[int]) -> Dict[int, Tuple[float, float, float]]:
    """CA atom position for each of `resnums` in the given chain, read
    directly from the PDB file text (fixed-column parsing, same convention
    as `pdbio.py`) -- py3Dmol/3Dmol.js has no synchronous "read back a
    loaded model's coordinates" API from Python, so label placement needs
    its own tiny parse of the same file `addModel` was given."""
    wanted = set(resnums)
    coords: Dict[int, Tuple[float, float, float]] = {}
    for line in Path(pdb_path).read_text().splitlines():
        if line[:6] != "ATOM  " or line[12:16].strip() != "CA" or line[21] != chain_id:
            continue
        try:
            resseq = int(line[22:26])
        except ValueError:
            continue
        if resseq in wanted and resseq not in coords:
            coords[resseq] = (float(line[30:38]), float(line[38:46]), float(line[46:54]))
    return coords


def _add_site_labels(view: py3Dmol.view, pdb_path: str, chain_id: str, site_labels: Dict[int, str], color: str) -> None:
    font_color = _readable_font_color(color)
    for resnum, coord in _ca_coords(pdb_path, chain_id, list(site_labels)).items():
        view.addLabel(site_labels.get(resnum, str(resnum)), {
            "position": {"x": coord[0], "y": coord[1], "z": coord[2]},
            "backgroundColor": color, "backgroundOpacity": 0.8,
            "fontColor": font_color, "fontSize": 11, "showBackground": True, "borderThickness": 0.4,
        })


def build_overlay_view(
    structures: Sequence[dict], reference_label: str, *,
    colors: Optional[Dict[str, str]] = None, width: Union[int, str] = "100%", height: Union[int, str] = 600,
    label_residues: bool = False,
) -> py3Dmol.view:
    """`structures`: a flat list of independent entries to draw -- normally
    one `kind="afdb"` entry per shown protein (its AlphaFold model) plus
    zero or more `kind="pdb"` entries (real RCSB structures, already
    superposed onto the reference frame -- see `pdbstruct.py`/
    `pipeline.py`'s `pdb_structures` report field). Keeping every real
    structure as its own list entry, rather than nesting them under a
    parent protein entry, lets the caller show/hide a protein's AlphaFold
    model and each of its real structures independently of one another.

    Each entry: `label` (used to look up a shared per-protein `color` from
    `colors` when `color` isn't given directly on the entry -- an afdb
    entry normally leaves `color` unset and inherits its protein's shared
    color this way, while a pdb entry normally passes its own explicit
    `color`, e.g. from `assign_pdb_colors`, so distinct real structures of
    the same protein don't all render identically), `pdb_path`, `chain_id`,
    `site_resseqs`/`site_labels` (that structure's *own* numbering -- for
    `kind="pdb"` see `pipeline.py`'s `pocket_resseq` mapping, since a real
    PDB entry's residue numbers don't equal canonical UniProt position,
    unlike an AlphaFold model; `site_labels` is an optional
    `{resnum: "K33"}` map, used only when `label_residues` is True). `kind`
    (`"afdb"`, the default, vs. `"pdb"` -- both render as a solid, fully
    opaque cartoon; `kind` only affects pocket-residue stick thickness and
    whether a bound ligand gets drawn, see below). For `kind="pdb"`
    entries only: `ligand_resname` (drawn as sticks in the fixed
    `LIGAND_COLOR`, not the entry's own cartoon `color`, if given -- the
    actual payoff of overlaying a real structure at all; a fixed color
    makes "where's the ligand" readable at a glance across every shown
    structure, rather than blending in with whichever cartoon color that
    structure happened to get)."""
    view = py3Dmol.view(width=width, height=height)
    colors = colors or assign_colors([s["label"] for s in structures], reference_label)

    for model_index, s in enumerate(structures):
        pdb_text = Path(s["pdb_path"]).read_text()
        view.addModel(pdb_text, "pdb")
        color = s.get("color") or colors.get(s["label"], "gray")
        kind = s.get("kind", "afdb")
        chain_sel = {"model": model_index, "chain": s["chain_id"]}

        # Every chain is explicitly styled to nothing first, since 3Dmol.js
        # falls back to a default line/wireframe rendering for any atom
        # left unstyled rather than hiding it.
        view.setStyle({"model": model_index}, {})
        view.setStyle(chain_sel, {"cartoon": {"color": color}})

        site = s.get("site_resseqs")
        if site:
            stick_radius = 0.12 if kind == "pdb" else 0.17
            view.addStyle(
                {"model": model_index, "chain": s["chain_id"], "resi": list(site)},
                {"stick": {"colorscheme": _carbon_tint_scheme(SITE_COLOR), "radius": stick_radius}},
            )
            if label_residues:
                site_labels = {r: (s.get("site_labels") or {}).get(r, str(r)) for r in site}
                _add_site_labels(view, s["pdb_path"], s["chain_id"], site_labels, color)

        if kind == "pdb":
            ligand_resname = s.get("ligand_resname")
            if ligand_resname:
                view.addStyle(
                    {"model": model_index, "resn": ligand_resname},
                    {"stick": {"colorscheme": _carbon_tint_scheme(LIGAND_COLOR), "radius": 0.35}},
                )

    view.zoomTo()
    return view
