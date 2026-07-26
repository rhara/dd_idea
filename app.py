"""Streamlit UI for dd_idea.

Run with `streamlit run app.py -- --report-dir data` (after `dd_idea-run`/
`dd_idea-fetch`+`dd_idea-align` have populated that directory with
`report.json` and `aligned/*_aligned.pdb`), or just `streamlit run app.py`
and enter the directory in the sidebar.
"""
import argparse
import json
import sys
from pathlib import Path

import streamlit as st

from dd_idea import dashboard, scene
from dd_idea.viewer3d import html_fill_container, html_with_camera_events, html_with_initial_view, view3d

st.set_page_config(page_title="dd_idea", layout="wide")

CONSERVATION_COLORS = {
    "identical": "#c8e6c9",
    "conservative": "#fff59d",
    "non-conservative": "#ef9a9a",
    "gap": "#e0e0e0",
}


def _parse_cli_defaults() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--report-dir")
    args, _ = parser.parse_known_args(sys.argv[1:])
    return args


@st.cache_data(show_spinner=False)
def _load_json(path: str, mtime: float) -> dict:
    # `mtime` deliberately has no leading underscore: Streamlit's
    # cache_data excludes underscore-prefixed params from the cache key
    # (meant for unhashable args like DB connections) -- naming it `_mtime`
    # would silently defeat the whole point of passing it at all.
    return json.loads(Path(path).read_text())


def _load_json_fresh(path: Path) -> dict:
    """`_load_json`, but keyed on the file's current mtime too -- a plain
    `@st.cache_data` keyed on `path` alone keeps serving the *first*
    version ever loaded for that path, forever, even after re-running
    `dd_idea-fetch`/`-align` writes a new `report.json` to the same
    file while this app's own process keeps running (Streamlit has no way
    to know the file changed otherwise). Re-fetching more real structures
    into an already-open report directory is an explicitly supported
    workflow (see README "Real-structure overlay"), so this cache must
    actually invalidate when the file does."""
    return _load_json(str(path), path.stat().st_mtime)


def main() -> None:
    defaults = _parse_cli_defaults()
    st.title("dd_idea -- cross-protein structure & active-site comparison")

    with st.sidebar:
        st.header("Report")
        report_dir = st.text_input("Report directory (dd_idea-run/-align output)", value=defaults.report_dir or "data")

    report_path = Path(report_dir) / "report.json"
    if not report_path.exists():
        st.info(f"No report.json found in {report_dir!r}. Run `dd_idea-run ACC1 ACC2 [...] -o {report_dir}` first.")
        st.stop()

    report = _load_json_fresh(report_path)
    proteins = report["proteins"]
    labels = [p["accession"] for p in proteins]
    reference = report["reference"]

    st.caption(
        f"Reference: {reference} -- {len(labels)} proteins -- "
        f"pocket: {len(report['pocket']['residues'])} residue(s) (druggability={report['pocket']['druggability_score']:.2f})"
    )

    by_label = {p["accession"]: p for p in proteins}

    with st.sidebar:
        st.header("Overlay display")
        st.write("Proteins to show")
        select_all_col, deselect_all_col = st.columns(2)
        if select_all_col.button("Select all", width="stretch"):
            for label in labels:
                st.session_state[f"show_prot_{label}"] = True
                for i in range(len(by_label[label].get("pdb_structures") or [])):
                    st.session_state[f"show_pdb_{label}_{i}"] = True
        if deselect_all_col.button("Deselect all", width="stretch"):
            for label in labels:
                st.session_state[f"show_prot_{label}"] = False
                for i in range(len(by_label[label].get("pdb_structures") or [])):
                    st.session_state[f"show_pdb_{label}_{i}"] = False
        for label in labels:
            st.checkbox(label, value=True, key=f"show_prot_{label}")
            pdb_structures = by_label[label].get("pdb_structures") or []
            if not pdb_structures:
                continue
            # Collapsed by default and behind an expander, not a flat
            # indented checkbox list -- a well-studied protein (e.g. CDK2
            # with hundreds of RCSB entries) can have dozens of these kept
            # even after --pdb-max-structures/--resolution-cutoff filtering,
            # and a flat list that long buries every other protein's own
            # controls below the fold.
            expander = st.expander(f"{len(pdb_structures)} real structure(s)", expanded=False)
            for i, pdb in enumerate(pdb_structures):
                ligand_note = pdb.get("ligand_resname") or "apo"
                res_note = f"{pdb['resolution']:.2f}Å" if pdb.get("resolution") is not None else "n/a"
                expander.checkbox(
                    f"{pdb['pdb_id']} ({ligand_note}, {res_note})", value=False, key=f"show_pdb_{label}_{i}",
                )
        selected = [label for label in labels if st.session_state[f"show_prot_{label}"]]
        show_pocket = st.checkbox("Highlight reference pocket residues", value=True)
        label_residues = st.checkbox("Show residue labels", value=True, disabled=not show_pocket)

        if "camera_generation" not in st.session_state:
            st.session_state.camera_generation = 0
        if st.button("Reset view"):
            st.session_state.camera_generation += 1
            st.session_state.camera_view = None

    candidates_path = Path(report_dir) / "candidates.json"
    tab_names = ["Overview", "Active-site comparison", "Structure overlay"]
    if candidates_path.exists():
        tab_names.append("Candidates")
    tabs = st.tabs(tab_names)

    with tabs[0]:
        summary = dashboard.overview_dataframe(report)
        st.dataframe(summary, width="stretch", height=(len(summary) + 1) * 35 + 3, hide_index=True)

    with tabs[1]:
        st.caption(
            "Every reference pocket residue, mapped onto each protein's own numbering. "
            "Color: green = identical, yellow = conservative substitution, red = non-conservative, gray = no counterpart (gap)."
        )
        values_df, cons_df = dashboard.pocket_comparison_frames(report)

        def _apply_colors(_df):
            return cons_df.map(lambda v: f"background-color: {CONSERVATION_COLORS.get(v, '')}")

        st.dataframe(values_df.style.apply(_apply_colors, axis=None), width="stretch", height=(len(values_df) + 1) * 35 + 3)

    with tabs[2]:
        # Canonical-position pocket labels ("K33") per protein, computed once
        # regardless of AFDB/PDB checkbox state -- a real-PDB checkbox can be
        # checked independently of that protein's own AlphaFold checkbox.
        canon_site_labels_by_label = {}
        if show_pocket:
            for label in labels:
                covered = [c for c in by_label[label]["pocket_comparison"] if c["target_position"] is not None]
                canon_site_labels_by_label[label] = {
                    c["target_position"]: f"{c['target_residue']}{c['target_position']}" for c in covered
                }

        all_pdb_entries = [
            (label, i, pdb) for label in labels for i, pdb in enumerate(by_label[label].get("pdb_structures") or [])
        ]
        pdb_colors = scene.assign_pdb_colors([f"{label}:{pdb['pdb_id']}" for label, i, pdb in all_pdb_entries])

        scene_structures = []
        for label in selected:
            p = by_label[label]
            if not p.get("aligned_pdb"):
                continue  # skipped during structural alignment (see p.get("align_error"))
            site_labels = canon_site_labels_by_label.get(label)
            site = list(site_labels) if site_labels is not None else None
            scene_structures.append({
                "label": label, "pdb_path": p["aligned_pdb"], "chain_id": p["chain"],
                "site_resseqs": site, "site_labels": site_labels, "kind": "afdb",
            })

        pdb_caption_lines = []
        for label, i, pdb in all_pdb_entries:
            if not st.session_state.get(f"show_pdb_{label}_{i}", True) or not pdb.get("aligned_pdb"):
                continue
            pdb_site, pdb_site_labels = None, None
            if show_pocket:
                site_labels = canon_site_labels_by_label.get(label) or {}
                pdb_site, pdb_site_labels = [], {}
                for canon_str, pdb_resseq in (pdb.get("pocket_resseq") or {}).items():
                    if pdb_resseq is None:
                        continue
                    pdb_site.append(pdb_resseq)
                    pdb_site_labels[pdb_resseq] = site_labels.get(int(canon_str), canon_str)
            structure_color = pdb_colors.get(f"{label}:{pdb['pdb_id']}")
            scene_structures.append({
                "label": label, "pdb_path": pdb["aligned_pdb"], "chain_id": pdb["chain"], "kind": "pdb",
                "site_resseqs": pdb_site, "site_labels": pdb_site_labels,
                "ligand_resname": pdb.get("ligand_resname"), "color": structure_color,
            })
            ligand_note = f"ligand {pdb['ligand_resname']}" if pdb.get("ligand_resname") else "apo"
            res_note = f"{pdb['resolution']:.2f}Å" if pdb.get("resolution") is not None else "resolution n/a"
            swatch = (
                f'<span style="display:inline-block;width:0.9em;height:0.9em;'
                f'background:{structure_color};border-radius:2px;margin-right:0.3em;'
                f'vertical-align:middle;"></span>'
            )
            pdb_caption_lines.append(f"{swatch}<b>{label}</b> {pdb['pdb_id']} ({res_note}, {ligand_note})")

        no_pdb_labels = [label for label in labels if not by_label[label].get("pdb_structures")]

        if pdb_caption_lines:
            st.markdown(
                "Real PDB structure overlay (thinner sticks): "
                + " &nbsp;|&nbsp; ".join(pdb_caption_lines),
                unsafe_allow_html=True,
            )
        if no_pdb_labels:
            st.caption(f"No RCSB structure found for: {', '.join(no_pdb_labels)} (AlphaFold model only)")

        if not scene_structures:
            st.info("No selected protein has a superposed coordinate file to show (all skipped during structural alignment?).")
        else:
            # Rebuilt -- and re-baked with the latest known camera view --
            # only when the *displayed content* actually changes, not on
            # every rerun. This matters because `view3d`'s return value
            # below feeds back into `st.session_state.camera_view`, which
            # this same block bakes into the HTML: rebuilding unconditionally
            # every run means a rerun triggered purely by the component
            # reporting its own camera position sends back a *fresh* HTML
            # string (`_make_html()` picks a new random viewer id every
            # call regardless of content -- see `htmlpatch.get_viewer_
            # variable`'s docstring), which the component always treats as
            # a new scene to load, whose own initial-paint self-report
            # triggers another rerun, forever -- confirmed live as a
            # permanent back-and-forth between the correct camera position
            # and the scene's own default fit that never settled. Reusing
            # the exact same HTML string when nothing structural changed
            # means the component's own `args.html !== currentHtml` check
            # sees no change and does nothing, breaking that loop.
            # `camera_generation` is part of the signature so "Reset view"
            # (which bumps it and clears `camera_view`) still forces a
            # fresh, un-rotated build rather than being cached away.
            scene_signature = (
                reference, label_residues, st.session_state.camera_generation,
                tuple(
                    (s["label"], s["pdb_path"], s["kind"], s.get("ligand_resname"), s.get("color"),
                     tuple(s.get("site_resseqs") or ()))
                    for s in scene_structures
                ),
            )
            if st.session_state.get("scene_signature") != scene_signature:
                view = scene.build_overlay_view(scene_structures, reference_label=reference, label_residues=label_residues)
                html = html_with_camera_events(
                    html_with_initial_view(html_fill_container(view._make_html()), st.session_state.get("camera_view"))
                )
                st.session_state.scene_signature = scene_signature
                st.session_state.scene_html = html
            else:
                html = st.session_state.scene_html
            # The component's return value is its own live-updated camera
            # view (see `component.view3d`'s docstring for why this, not
            # just the component's in-memory `lastView`, is what actually
            # survives a widget-triggered rerun that tears the component's
            # iframe down and recreates it). Stashed in session_state so the
            # *next rebuild* above can bake it into a fresh scene's HTML
            # before it ever reaches the browser.
            reported_view = view3d(html, height=650, reset_camera_token=st.session_state.camera_generation)
            if reported_view:
                st.session_state.camera_view = reported_view

        skipped = [p for p in proteins if p["accession"] in selected and not p.get("aligned_pdb")]
        if skipped:
            st.warning(
                "Not shown (couldn't be structurally superposed): "
                + ", ".join(f"{p['accession']} ({p['align_error']})" for p in skipped)
            )

    if candidates_path.exists():
        with tabs[3]:
            candidates_report = _load_json_fresh(candidates_path)
            st.caption(
                f"Seed: {candidates_report['seed_accession']} ({candidates_report['seed_name']}) -- "
                f"family {candidates_report['family_database']}:{candidates_report['family_id']} "
                f"({candidates_report['family_entry_name']}, {candidates_report['family_member_count']} member(s))"
            )
            cand_df = dashboard.candidates_dataframe(candidates_report)
            st.dataframe(cand_df, width="stretch", height=(len(cand_df) + 1) * 35 + 3, hide_index=True)


if __name__ == "__main__":
    main()
