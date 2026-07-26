"""Picking a small, illustrative set of real RCSB structures for a
protein already fetched as an AlphaFold model (see README "Why always
AlphaFold, never a PDB structure" and its "real-structure overlay" caveat)
-- an *additive* visualization layer only; pocket detection and the
cross-protein sequence/pocket mapping in `pocketmap.py`/`sequence.py` stay
anchored on the AlphaFold model unconditionally. Builds on `rcsb.fetch`
(lookup/download) and `rcsb.chain` (canonical-sequence alignment).
"""
from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple, Union

from .. import pdbio
from .chain import ChainAlignment, align_to_canonical, extract_chain_sequences, pick_target_chain
from .fetch import EntryMetadata, download_pdb, fetch_entry_metadata, list_pdb_ids_for_uniprot


@dataclass
class SelectedPdbStructure:
    accession: str
    pdb_id: str
    resolution: Optional[float]
    chain_id: str
    pdb_path: str
    ligand_resname: Optional[str]  # None if apo
    chain_alignment: ChainAlignment
    n_candidates_scanned: int


def select_pdb_structures(
    accession: str, canonical_seq: str, out_dir: Union[str, Path], *,
    scan_cap: int = 25, min_ligand_atoms: int = 5, max_structures: int = 3,
    resolution_cutoff: float = 2.0, dedupe_ligands: bool = True, show_progress: bool = True,
) -> List["SelectedPdbStructure"]:
    """Pick up to `max_structures` real PDB structures for `accession`,
    preferring ones with a genuine bound ligand (not water/cryoprotectant/
    cofactor -- see `pdbio.classify_hetero_groups`). By default
    (`dedupe_ligands=True`) picks one per *distinct* ligand (deduplicated by
    resname, resolution-ranked-first-seen wins) so a target crystallized
    many times with the same fragment doesn't crowd out other chemical
    matter -- set `dedupe_ligands=False` to keep every ligand-bound entry
    that meets `resolution_cutoff`, including repeat entries for the same
    ligand (e.g. wanting every sub-1.5Å structure regardless of how many
    share a ligand). Falls back to a single best-resolution entry overall if
    none of the (resolution-sorted) first `scan_cap` candidates has a ligand
    at all. Returns `[]` if the accession has no RCSB structures. Candidates
    are scanned cheapest-metadata-first and downloaded/parsed for ligand
    content only as needed, stopping once `max_structures` entries are kept
    or `scan_cap` candidates have been scanned -- a well-studied target
    (e.g. CDK2's 512 entries) would otherwise mean downloading hundreds of
    structures; pass a large `max_structures`/`scan_cap` explicitly to
    actually collect that many.

    `resolution_cutoff` (Angstrom, lower is better -- 2.0 by default)
    excludes any candidate whose resolution is worse than this or entirely
    unreported (e.g. NMR structures don't have one) *before* downloading
    its coordinates at all, for both the ligand-bound and apo-fallback
    paths -- this is a quality floor, not just a ranking tiebreaker, so a
    low-resolution entry is never picked just because it happened to be
    the first ligand-bound one scanned."""
    out_dir = Path(out_dir)
    pdb_ids = list_pdb_ids_for_uniprot(accession)  # already best-resolution-first (server-side sort)
    if not pdb_ids:
        if show_progress:
            print(f"[rcsb.select] {accession}: no RCSB structures found", flush=True)
        return []

    selections: List[SelectedPdbStructure] = []
    seen_ligands: set = set()
    best_apo: Optional[Tuple[EntryMetadata, str]] = None
    scanned = 0
    for pdb_id in pdb_ids[:scan_cap]:
        if len(selections) >= max_structures:
            break
        scanned += 1
        try:
            meta = fetch_entry_metadata(pdb_id)
        except Exception:
            continue  # a single entry's metadata being unfetchable shouldn't abort the whole scan
        if meta.resolution is None or meta.resolution > resolution_cutoff:
            continue  # doesn't meet the quality bar; skip without downloading coordinates
        dest = out_dir / f"{meta.pdb_id}.pdb"
        try:
            text = download_pdb(meta.pdb_id, dest)
        except Exception:
            continue
        groups = pdbio.classify_hetero_groups(pdbio.collect_hetero_groups(text))
        ligand = pdbio.pick_ligand_of_interest(groups, min_atoms=min_ligand_atoms)
        if ligand is not None:
            if dedupe_ligands and ligand.resname in seen_ligands:
                continue  # already have this ligand from a better-resolution entry; prefer diversity
            seen_ligands.add(ligand.resname)
            selections.append(_finalize_selection(accession, meta, dest, canonical_seq, ligand.resname, scanned))
            if show_progress:
                print(
                    f"[rcsb.select] {accession}: picked {meta.pdb_id} (resolution={meta.resolution}, "
                    f"ligand={ligand.resname}) [{len(selections)}/{max_structures}] after scanning "
                    f"{scanned} candidate(s)", flush=True,
                )
            continue
        if best_apo is None:
            best_apo = (meta, str(dest))

    if selections:
        return selections
    if best_apo is None:
        if show_progress:
            print(
                f"[rcsb.select] {accession}: {scanned} candidate(s) scanned, none met the "
                f"resolution_cutoff={resolution_cutoff}Å bar (or were downloadable)", flush=True,
            )
        return []
    meta, dest = best_apo
    if show_progress:
        print(
            f"[rcsb.select] {accession}: no ligand-bound entry among {scanned} scanned candidate(s), "
            f"falling back to best resolution {meta.pdb_id} ({meta.resolution})", flush=True,
        )
    return [_finalize_selection(accession, meta, Path(dest), canonical_seq, None, scanned)]


def _finalize_selection(
    accession: str, meta: EntryMetadata, pdb_path: Path, canonical_seq: str,
    ligand_resname: Optional[str], scanned: int,
) -> SelectedPdbStructure:
    chains = extract_chain_sequences(pdb_path)
    alignments = {cid: align_to_canonical(cs, canonical_seq) for cid, cs in chains.items()}
    target_chain = pick_target_chain(alignments)
    return SelectedPdbStructure(
        accession=accession, pdb_id=meta.pdb_id, resolution=meta.resolution, chain_id=target_chain,
        pdb_path=str(pdb_path), ligand_resname=ligand_resname, chain_alignment=alignments[target_chain],
        n_candidates_scanned=scanned,
    )


def rehydrate_selection(
    accession: str, canonical_seq: str, *, pdb_id: str, resolution: Optional[float], chain_id: str,
    pdb_path: str, ligand_resname: Optional[str],
) -> SelectedPdbStructure:
    """Reconstruct a `SelectedPdbStructure` from a previously-cached pick
    (`fetch_all`'s manifest.json entry) without any network access -- just
    re-reads the already-downloaded `pdb_path` and re-runs the local
    canonical-sequence alignment for the (already-known) `chain_id`. Used
    by `pipeline.analyze` so re-running the align step never re-hits RCSB:
    the real network lookup/selection only happens once, in `fetch_all`,
    exactly like canonical-sequence/AlphaFold-model fetching already
    works -- `select_pdb_structures` itself is for that fetch-time
    lookup, not for align-time reuse."""
    chains = extract_chain_sequences(pdb_path)
    alignment = align_to_canonical(chains[chain_id], canonical_seq)
    return SelectedPdbStructure(
        accession=accession, pdb_id=pdb_id, resolution=resolution, chain_id=chain_id, pdb_path=pdb_path,
        ligand_resname=ligand_resname, chain_alignment=alignment, n_candidates_scanned=0,
    )


@dataclass
class PdbOverlayResult:
    accession: str
    pdb_id: str
    resolution: Optional[float]
    ligand_resname: Optional[str]
    chain_id: str
    aligned_pdb: str
    rmsd: Optional[float]
    n_aligned_atoms: int
    error: Optional[str] = None


def align_pdb_overlays(
    reference_afdb_path: Union[str, Path], reference_label: str,
    selections: Dict[str, List[SelectedPdbStructure]], out_dir: Union[str, Path], *, show_progress: bool = True,
) -> Dict[str, List[PdbOverlayResult]]:
    """Superpose every protein's selected real PDB structure(s) (`selections`,
    keyed by accession -- proteins with no RCSB structure are simply absent
    or map to an empty list) onto the reference protein's AlphaFold-model
    frame via whole-chain `cealign`, reusing `structalign.align_structures`
    unmodified (it doesn't care whether an input is an AlphaFold model or a
    real structure). The reference's own AFDB file is the fixed anchor; if
    the reference protein itself has selected PDB structure(s) too, those
    get cealigned onto its own AFDB model just like every other protein's
    does -- a real PDB entry's coordinate frame has no relation to
    AlphaFold's, even for the same protein, so each still needs its own
    fit. Each structure gets a unique `structalign` label
    (`f"{acc}_pdb{i}"`) so multiple entries for the same protein don't
    collide, and so `structalign`'s output filenames stay distinct."""
    from .. import structalign

    if not any(selections.values()):
        return {acc: [] for acc in selections}

    anchor = structalign.StructureInput(label=reference_label, pdb_path=str(reference_afdb_path), chain_id="A")
    mobiles = [
        structalign.StructureInput(label=f"{acc}_pdb{i}", pdb_path=sel.pdb_path, chain_id=sel.chain_id)
        for acc, sels in selections.items() for i, sel in enumerate(sels)
    ]
    results = structalign.align_structures(
        [anchor] + mobiles, reference_label=reference_label, out_dir=out_dir, show_progress=show_progress,
    )
    by_label = {r.label: r for r in results}

    out: Dict[str, List[PdbOverlayResult]] = {}
    for acc, sels in selections.items():
        entries = []
        for i, sel in enumerate(sels):
            r = by_label.get(f"{acc}_pdb{i}")
            if r is None:
                continue
            entries.append(PdbOverlayResult(
                accession=acc, pdb_id=sel.pdb_id, resolution=sel.resolution, ligand_resname=sel.ligand_resname,
                chain_id=sel.chain_id, aligned_pdb=r.aligned_pdb, rmsd=r.rmsd, n_aligned_atoms=r.n_atoms,
                error=r.error,
            ))
        out[acc] = entries
    return out
