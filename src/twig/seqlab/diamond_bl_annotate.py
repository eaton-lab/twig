#!/usr/bin/env python

"""Annotation loading and expansion helpers for diamond-bl."""

from __future__ import annotations

import csv
import gzip
import re
from pathlib import Path

GO_DEFAULT_COLS = [
    "go_id",
    "go_aspect",
    "go_qualifier",
    "go_evidence",
    "go_reference",
    "go_with_from",
    "go_assigned_by",
]
GAF_EXTRA_COLS = [
    "db_object_id",
    "db_object_symbol",
    "db_object_name",
]
GAF_DEFAULT_COLS = GO_DEFAULT_COLS + GAF_EXTRA_COLS

ISOFORM_SUFFIX_RE = re.compile(r"\.\d+$")


def _open_text(path: Path):
    """Open plain text or gzip-compressed text files."""
    if path.suffix == ".gz":
        return gzip.open(path, "rt")
    return path.open()


def detect_annotation_format(path: Path, annotation_format: str) -> str:
    """Resolve annotation format from explicit choice or path suffix."""
    fmt = annotation_format.lower()
    if fmt != "auto":
        if fmt not in {"tsv", "gaf", "gpa"}:
            raise ValueError(f"unsupported --annotation-format: {annotation_format!r}")
        return fmt

    name = path.name.lower()
    if name.endswith(".gaf") or name.endswith(".gaf.gz"):
        return "gaf"
    if name.endswith(".gpa") or name.endswith(".gpa.gz"):
        return "gpa"
    return "tsv"


def _select_go_cols(wanted_cols: list[str] | None) -> list[str]:
    """Validate and return requested GO output columns."""
    if wanted_cols:
        missing = [col for col in wanted_cols if col not in GO_DEFAULT_COLS]
        if missing:
            raise ValueError(
                f"invalid GO annotation columns: {missing}. "
                f"Valid columns are: {GO_DEFAULT_COLS}"
            )
        return list(wanted_cols)
    return list(GO_DEFAULT_COLS)


def _select_gaf_cols(wanted_cols: list[str] | None) -> list[str]:
    """Validate and return requested GAF output columns."""
    if wanted_cols:
        missing = [col for col in wanted_cols if col not in GAF_DEFAULT_COLS]
        if missing:
            raise ValueError(
                f"invalid GAF annotation columns: {missing}. "
                f"Valid columns are: {GAF_DEFAULT_COLS}"
            )
        return list(wanted_cols)
    return list(GAF_DEFAULT_COLS)


def _load_tsv_annotations(
    path: Path,
    key_col: str,
    wanted_cols: list[str] | None,
) -> tuple[dict[str, list[dict[str, str]]], list[str]]:
    """Load generic TSV annotations keyed by a user-defined column."""
    table: dict[str, list[dict[str, str]]] = {}
    with _open_text(path) as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if reader.fieldnames is None or key_col not in reader.fieldnames:
            raise ValueError(f"annotation file must contain a '{key_col}' column: {path}")

        if wanted_cols:
            for col in wanted_cols:
                if col not in reader.fieldnames:
                    raise ValueError(f"annotation column '{col}' not found in {path}")
            cols = list(wanted_cols)
        else:
            cols = [name for name in reader.fieldnames if name != key_col]

        for row in reader:
            key = row.get(key_col)
            if not key:
                continue
            rec = {col: row.get(col, "") for col in cols}
            table.setdefault(key, []).append(rec)

    return table, cols


def _load_gaf_annotations(
    path: Path,
    wanted_cols: list[str] | None,
) -> tuple[dict[str, list[dict[str, str]]], list[str]]:
    """Load GO annotations from GAF keyed by DB_Object_ID."""
    cols = _select_gaf_cols(wanted_cols)
    table: dict[str, list[dict[str, str]]] = {}

    with _open_text(path) as handle:
        for i, line in enumerate(handle, start=1):
            if not line or line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 15:
                raise ValueError(f"malformed GAF row at line {i}: expected >=15 columns, found {len(parts)}")

            key = parts[1].strip()  # DB_Object_ID
            if not key:
                continue

            rec = {
                "db_object_id": parts[1].strip(),
                "db_object_symbol": parts[2].strip(),
                "db_object_name": parts[9].strip(),
                "go_id": parts[4].strip(),
                "go_aspect": parts[8].strip(),
                "go_qualifier": parts[3].strip(),
                "go_evidence": parts[6].strip(),
                "go_reference": parts[5].strip(),
                "go_with_from": parts[7].strip(),
                "go_assigned_by": parts[14].strip(),
            }
            table.setdefault(key, []).append({col: rec.get(col, "") for col in cols})

    return table, cols


def _load_gpa_annotations(
    path: Path,
    wanted_cols: list[str] | None,
) -> tuple[dict[str, list[dict[str, str]]], list[str]]:
    """Load GO annotations from GPA keyed by DB_Object_ID."""
    cols = _select_go_cols(wanted_cols)
    table: dict[str, list[dict[str, str]]] = {}

    with _open_text(path) as handle:
        for i, line in enumerate(handle, start=1):
            if not line or line.startswith("!"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 10:
                raise ValueError(f"malformed GPA row at line {i}: expected >=10 columns, found {len(parts)}")

            key = parts[1].strip()  # DB_Object_ID
            if not key:
                continue

            rec = {
                "go_id": parts[3].strip(),
                "go_aspect": "",
                "go_qualifier": parts[2].strip(),
                "go_evidence": parts[5].strip(),
                "go_reference": parts[4].strip(),
                "go_with_from": parts[6].strip(),
                "go_assigned_by": parts[9].strip(),
            }
            table.setdefault(key, []).append({col: rec.get(col, "") for col in cols})

    return table, cols


def load_annotations(
    path: Path,
    annotation_format: str,
    key_col: str,
    wanted_cols: list[str] | None,
) -> tuple[dict[str, list[dict[str, str]]], list[str]]:
    """Load annotations table in TSV/GAF/GPA formats."""
    fmt = detect_annotation_format(path, annotation_format)

    if fmt == "tsv":
        return _load_tsv_annotations(path, key_col, wanted_cols)

    if key_col != "sseqid":
        raise ValueError("--annotation-key is only supported for TSV annotations; use default key for GAF/GPA")

    if fmt == "gaf":
        return _load_gaf_annotations(path, wanted_cols)
    if fmt == "gpa":
        return _load_gpa_annotations(path, wanted_cols)

    raise ValueError(f"unsupported annotation format: {fmt!r}")


def _fallback_gene_key_from_isoform_id(key: str) -> str:
    """Return a gene-level fallback key by stripping a trailing '.<digits>' suffix."""
    return ISOFORM_SUFFIX_RE.sub("", key)


def expand_rows_with_annotations(
    rows: list[dict[str, str]],
    annots: dict[str, list[dict[str, str]]] | None,
    annot_cols: list[str] | None,
    subject_key: str = "sseqid",
    match_mode: str = "exact",
) -> list[dict[str, str]]:
    """Expand hit rows to one row per matching annotation record."""
    if not rows:
        return []
    if not annots or not annot_cols:
        return [dict(row) for row in rows]
    if match_mode not in {"exact", "exact_then_gene_fallback"}:
        raise ValueError(f"unsupported annotation match mode: {match_mode!r}")

    expanded: list[dict[str, str]] = []
    empty_annot = {col: "" for col in annot_cols}

    for row in rows:
        sid = row.get(subject_key, "")
        matches = annots.get(sid, [])
        if not matches and match_mode == "exact_then_gene_fallback":
            fallback = _fallback_gene_key_from_isoform_id(sid)
            if fallback != sid:
                matches = annots.get(fallback, [])
        if not matches:
            rec = dict(row)
            rec.update(empty_annot)
            expanded.append(rec)
            continue

        for match in matches:
            rec = dict(row)
            rec.update({col: match.get(col, "") for col in annot_cols})
            expanded.append(rec)

    return expanded
