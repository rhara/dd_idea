"""End-to-end orchestration in three phases:

- `fetch_all`: given a list of UniProt accessions, download each protein's
  canonical sequence and AlphaFold DB model, and (unless disabled) look up
  and select its real RCSB structures too -- all genuinely fetch-time
  network work, all cached on disk -- writing `manifest.json`.
- `discover`: given one seed accession, propose candidate similar proteins
  (`similarity.py`), writing `candidates.json` -- a human-reviewed proposal,
  not automatically fed into `fetch_all`.
- `analyze`: pick a reference protein, detect its druggable pocket, align
  every other protein's canonical sequence to the reference and map the
  pocket onto each, then superpose every AlphaFold model -- and every
  already-selected real structure, rehydrated from `manifest.json` with no
  network access -- onto the reference via whole-chain `cealign`, writing
  `report.json` plus the superposed coordinate files under `aligned/`.
  Real-structure *selection* is deliberately not repeated here: it's
  fetch-time work, done once in `fetch_all`, exactly like canonical-
  sequence/AlphaFold-model fetching -- re-running `analyze` (e.g. to try a
  different `--reference`/`--pocket-rank`) never re-hits RCSB.
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Dict, List, Optional, Sequence, Union

from . import alphafold, pocketmap, rcsb, sequence, similarity, structalign, uniprot

AFDB_CHAIN = "A"  # every AlphaFold DB model is a single chain, always "A"


def _read_fasta(path: Union[str, Path]) -> str:
    return "".join(Path(path).read_text().splitlines()[1:])


def fetch_all(
    accessions: Sequence[str], out_dir: Union[str, Path], *, show_progress: bool = True,
    pdb_overlay: bool = True, pdb_scan_cap: int = 25, pdb_max_structures: int = 3,
    pdb_resolution_cutoff: float = 2.0, pdb_dedupe_ligands: bool = True,
) -> dict:
    """Download every protein's canonical sequence + AlphaFold DB model,
    and (unless `pdb_overlay` is False) look up and select its real RCSB
    structures too (see `rcsb.select_pdb_structures`). Cached:
    re-running against the same `out_dir` skips anything already on disk
    -- including individual real-structure downloads, so raising
    `pdb_max_structures` or loosening `pdb_resolution_cutoff` on a later
    call cheaply reuses whatever was already fetched and only spends new
    network time on the additional candidates needed to satisfy the wider
    request. `pdb_dedupe_ligands=False` keeps every ligand-bound entry
    meeting `pdb_resolution_cutoff` instead of one per distinct ligand --
    combine with a large `pdb_max_structures`/`pdb_scan_cap` to actually
    collect e.g. "every sub-1.5Å structure, duplicate ligands included"."""
    out_dir = Path(out_dir).resolve()
    raw_dir = out_dir / "raw"
    raw_dir.mkdir(parents=True, exist_ok=True)

    proteins: List[dict] = []
    for i, accession in enumerate(accessions, start=1):
        acc = accession.upper()
        entry = uniprot.fetch_uniprot_entry(acc)
        name = uniprot.protein_name(entry)

        fasta_dest = out_dir / f"{acc}.fasta"
        if fasta_dest.exists():
            canonical = _read_fasta(fasta_dest)
            if show_progress:
                print(f"[fetch] ({i}/{len(accessions)}) {acc}: canonical sequence already downloaded, skipping", flush=True)
        else:
            canonical = uniprot.fetch_uniprot_fasta(acc)
            fasta_dest.write_text(f">{acc}\n{canonical}\n")
            if show_progress:
                print(f"[fetch] ({i}/{len(accessions)}) {acc} ({name}): canonical sequence ({len(canonical)} aa)", flush=True)

        afdb_dest = raw_dir / f"{acc}_AFDB.pdb"
        already_had_it = afdb_dest.exists()
        alphafold.download_afdb(acc, afdb_dest)
        if show_progress:
            status = "already downloaded, skipping" if already_had_it else "-> " + afdb_dest.name
            print(f"[fetch] ({i}/{len(accessions)}) {acc}: AlphaFold DB model {status}", flush=True)

        pdb_structures: List[dict] = []
        if pdb_overlay:
            try:
                selections = rcsb.select_pdb_structures(
                    acc, canonical, out_dir / "raw_pdb", scan_cap=pdb_scan_cap,
                    max_structures=pdb_max_structures, resolution_cutoff=pdb_resolution_cutoff,
                    dedupe_ligands=pdb_dedupe_ligands, show_progress=show_progress,
                )
            except Exception as e:
                if show_progress:
                    print(f"[fetch] ({i}/{len(accessions)}) {acc}: real-structure lookup failed ({e}), skipping", flush=True)
                selections = []
            pdb_structures = [
                {
                    "pdb_id": s.pdb_id, "resolution": s.resolution, "chain": s.chain_id,
                    "pdb_path": s.pdb_path, "ligand_resname": s.ligand_resname,
                }
                for s in selections
            ]

        proteins.append({
            "accession": acc, "name": name, "length": len(canonical), "fasta_path": str(fasta_dest),
            "afdb_path": str(afdb_dest), "pdb_structures": pdb_structures,
        })

    manifest = {"proteins": proteins}
    (out_dir / "manifest.json").write_text(json.dumps(manifest, indent=2))
    return manifest


def discover(
    seed_accession: str, out_dir: Union[str, Path], *,
    max_candidates: int = 20, any_organism: bool = False, show_progress: bool = True,
) -> dict:
    """Propose candidate similar proteins for `seed_accession`, writing
    `candidates.json`. A proposal step only -- pass the accessions you want
    to actually keep to `fetch_all` yourself."""
    out_dir = Path(out_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)
    result = similarity.discover(
        seed_accession, max_candidates=max_candidates, any_organism=any_organism, show_progress=show_progress,
    )
    (out_dir / "candidates.json").write_text(json.dumps(result, indent=2))
    return result


def analyze(
    out_dir: Union[str, Path], *, reference: Optional[str] = None, pocket_rank: int = 1, show_progress: bool = True,
) -> dict:
    """Run pocket detection + cross-protein alignment + structural overlay
    across every protein `fetch_all` downloaded into `out_dir`. `reference`
    defaults to the first accession in `manifest.json`."""
    out_dir = Path(out_dir).resolve()
    manifest = json.loads((out_dir / "manifest.json").read_text())
    proteins: Dict[str, dict] = {p["accession"]: p for p in manifest["proteins"]}
    accessions = list(proteins)

    reference = (reference or accessions[0]).upper()
    if reference not in proteins:
        raise ValueError(f"reference {reference!r} not found among fetched proteins: {accessions}")
    if len(proteins) < 2:
        raise ValueError(f"need at least 2 proteins to compare, found {len(proteins)}: {accessions}")

    ref_info = proteins[reference]
    ref_seq = _read_fasta(ref_info["fasta_path"])

    if show_progress:
        print(f"[align] detecting a druggable pocket on reference {reference} ({ref_info['name']})...", flush=True)
    pocket = pocketmap.detect_reference_pocket(
        ref_info["afdb_path"], out_dir / "pocket", pocket_rank=pocket_rank, show_progress=show_progress,
    )
    if show_progress:
        print(f"[align] reference pocket: {len(pocket.residues)} lining residue(s)", flush=True)

    comparisons_by_acc: Dict[str, list] = {}
    identity_by_acc: Dict[str, float] = {reference: 100.0}
    canonical_by_acc: Dict[str, str] = {reference: ref_seq}
    for acc, info in proteins.items():
        if acc == reference:
            # Trivial self-comparison, so the dashboard can treat every
            # protein (including the reference) uniformly instead of
            # special-casing one column.
            comparisons_by_acc[acc] = [
                pocketmap.PocketResidueComparison(
                    reference_position=r.resnum, reference_residue=ref_seq[r.resnum - 1],
                    target_position=r.resnum, target_residue=ref_seq[r.resnum - 1], conservation="identical",
                )
                for r in pocket.residues
            ]
            continue
        target_seq = _read_fasta(info["fasta_path"])
        canonical_by_acc[acc] = target_seq
        mapping = sequence.align_to_reference(ref_seq, reference, target_seq, acc)
        comparisons_by_acc[acc] = pocketmap.compare_pocket(pocket, mapping)
        identity_by_acc[acc] = mapping.pct_identity
        if show_progress:
            n_noncons = sum(1 for c in comparisons_by_acc[acc] if c.conservation == "non-conservative")
            n_gap = sum(1 for c in comparisons_by_acc[acc] if c.conservation == "gap")
            print(
                f"[align] {acc}: {mapping.pct_identity:.1f}% identity to {reference}; "
                f"pocket: {n_noncons} non-conservative, {n_gap} gap (of {len(pocket.residues)})", flush=True,
            )

    structures = [
        structalign.StructureInput(label=acc, pdb_path=info["afdb_path"], chain_id=AFDB_CHAIN)
        for acc, info in proteins.items()
    ]
    align_results = {
        r.label: r for r in structalign.align_structures(
            structures, reference_label=reference, out_dir=out_dir / "aligned", show_progress=show_progress,
        )
    }

    pdb_selections: Dict[str, List[rcsb.SelectedPdbStructure]] = {}
    for acc, info in proteins.items():
        sels = []
        for cached in info.get("pdb_structures") or []:
            try:
                sels.append(rcsb.rehydrate_selection(
                    acc, canonical_by_acc[acc], pdb_id=cached["pdb_id"], resolution=cached["resolution"],
                    chain_id=cached["chain"], pdb_path=cached["pdb_path"], ligand_resname=cached["ligand_resname"],
                ))
            except Exception as e:
                if show_progress:
                    print(f"[align] {acc}: couldn't reload cached real structure {cached.get('pdb_id')} ({e}), skipping", flush=True)
        pdb_selections[acc] = sels
    pdb_overlay_results = rcsb.align_pdb_overlays(
        proteins[reference]["afdb_path"], reference, pdb_selections,
        out_dir / "aligned_pdb", show_progress=show_progress,
    )

    def _pdb_report(acc: str) -> list:
        out = []
        for sel, r in zip(pdb_selections.get(acc, []), pdb_overlay_results.get(acc, [])):
            pocket_resseq = {
                c.target_position: sel.chain_alignment.resseq_for_canonical(c.target_position)
                for c in comparisons_by_acc[acc] if c.target_position is not None
            }
            out.append({
                "pdb_id": r.pdb_id, "resolution": r.resolution, "ligand_resname": r.ligand_resname,
                "chain": r.chain_id, "aligned_pdb": r.aligned_pdb, "rmsd": r.rmsd,
                "n_aligned_atoms": r.n_aligned_atoms, "error": r.error, "pocket_resseq": pocket_resseq,
            })
        return out

    report = {
        "reference": reference,
        "pocket": pocket.to_report_dict(),
        "proteins": [
            {
                "accession": acc,
                "name": info["name"],
                "length": info["length"],
                "pct_identity": identity_by_acc[acc],
                "rmsd": align_results[acc].rmsd,
                "n_aligned_atoms": align_results[acc].n_atoms,
                "chain": AFDB_CHAIN,
                "aligned_pdb": align_results[acc].aligned_pdb,
                "align_error": align_results[acc].error,
                "pocket_comparison": [c.__dict__ for c in comparisons_by_acc[acc]],
                "pdb_structures": _pdb_report(acc),
            }
            for acc, info in proteins.items()
        ],
    }
    (out_dir / "report.json").write_text(json.dumps(report, indent=2))
    return report
