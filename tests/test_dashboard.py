"""Fast, offline unit tests for dashboard.py's DataFrame construction --
no network access."""
import pandas as pd
import pyarrow as pa

from dd_idea.dashboard import overview_dataframe

REPORT = {
    "reference": "Q8IZL9",
    "proteins": [
        {"accession": "Q8IZL9", "name": "Cyclin-dependent kinase 20", "length": 346,
         "pct_identity": 100.0, "rmsd": None, "align_error": None},
        {"accession": "Q00535", "name": "Cyclin-dependent kinase 5", "length": 292,
         "pct_identity": 46.4, "rmsd": 1.971, "align_error": None},
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
