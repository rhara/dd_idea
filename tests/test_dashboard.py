"""Fast, offline unit tests for dashboard.py's DataFrame construction --
no network access."""
import pandas as pd
import pyarrow as pa

from dd_idea.dashboard import overview_dataframe, pocket_comparison_frames

REPORT = {
    "reference": "Q8IZL9",
    "proteins": [
        {"accession": "Q8IZL9", "name": "Cyclin-dependent kinase 20", "length": 346,
         "pct_identity": 100.0, "rmsd": None, "align_error": None},
        {"accession": "Q00535", "name": "Cyclin-dependent kinase 5", "length": 292,
         "pct_identity": 46.4, "rmsd": 1.971, "align_error": None},
    ],
}

POCKET_REPORT = {
    "proteins": [
        {
            "accession": "Q8IZL9", "name": "Cyclin-dependent kinase 20",
            "pocket_comparison": [
                {"reference_residue": "Y", "reference_position": 15, "target_residue": "Y", "target_position": 15, "conservation": "identical"},
                {"reference_residue": "K", "reference_position": 33, "target_residue": "K", "target_position": 33, "conservation": "identical"},
            ],
        },
        {
            "accession": "P24941", "name": "Cyclin-dependent kinase 2",
            "pocket_comparison": [
                {"reference_residue": "Y", "reference_position": 15, "target_residue": "F", "target_position": 14, "conservation": "conservative"},
                {"reference_residue": "K", "reference_position": 33, "target_residue": None, "target_position": None, "conservation": "gap"},
            ],
        },
    ],
}


def test_overview_dataframe_reference_row_has_no_identity_string():
    df = overview_dataframe(REPORT)
    ref_row = df[df["Accession"] == "Q8IZL9"].iloc[0]
    assert pd.isna(ref_row["% identity to reference"])
    assert ref_row["Note"] == "(reference)"


def test_overview_dataframe_non_reference_row_keeps_identity_value():
    df = overview_dataframe(REPORT)
    row = df[df["Accession"] == "Q00535"].iloc[0]
    assert row["% identity to reference"] == 46.4
    assert row["Note"] == ""


def test_overview_dataframe_identity_column_is_arrow_serializable():
    # Regression test: mixing the string "(reference)" into this column
    # (the previous behavior) made pandas fall back to an `object` dtype
    # that pyarrow's Table.from_pandas rejects outright for a None/float
    # mix under the hood -- this is the exact conversion Streamlit's
    # dataframe rendering performs, and its failure was silently "fixed" by
    # rerunning the whole app script (see overview_dataframe's docstring).
    df = overview_dataframe(REPORT)
    pa.Table.from_pandas(df)  # raises pyarrow.lib.ArrowTypeError if broken


def test_pocket_comparison_frames_rows_are_proteins_columns_are_residues():
    # Transposed from the original row=residue/column=protein layout -- a
    # user-requested change so each protein's own mapped residues read
    # left-to-right like a short sequence instead of top-to-bottom.
    values_df, cons_df = pocket_comparison_frames(POCKET_REPORT)
    assert list(values_df.index) == ["Q8IZL9 (Cyclin-dependent kinase 20)", "P24941 (Cyclin-dependent kinase 2)"]
    assert list(values_df.columns) == ["Y15", "K33"]
    assert values_df.shape == (2, 2)
    assert cons_df.shape == (2, 2)


def test_pocket_comparison_frames_values_and_conservation_match_source():
    values_df, cons_df = pocket_comparison_frames(POCKET_REPORT)
    cdk2_row = values_df.loc["P24941 (Cyclin-dependent kinase 2)"]
    assert cdk2_row["Y15"] == "F14"
    assert cdk2_row["K33"] == "-"  # gap: no target residue
    cons_row = cons_df.loc["P24941 (Cyclin-dependent kinase 2)"]
    assert cons_row["Y15"] == "conservative"
    assert cons_row["K33"] == "gap"
