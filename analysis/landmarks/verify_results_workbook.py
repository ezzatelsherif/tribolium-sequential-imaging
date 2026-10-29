#!/usr/bin/env python3
"""Verify the submitted S1 Data workbook against its three source CSV tables.

Uses only the Python standard library. This is a read-only verification tool.
It does not rerun biological analyses or modify the workbook.
"""

from __future__ import annotations

import argparse
import collections
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import xml.etree.ElementTree as ET
import zipfile


NS = {
    "m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "r": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}
NUMBER = re.compile(r"^[+-]?(?:\d+\.?\d*|\.\d+)(?:[eE][+-]?\d+)?$")
SHEETS = [
    "Overall summary", "By landmark", "Landmark inventory", "Original peak audit",
    "Summary source", "ReadMe", "Data dictionary",
]


def typed(value: str):
    if value == "":
        return None
    if value in ("True", "False"):
        return value == "True"
    if NUMBER.fullmatch(value):
        return int(value) if re.fullmatch(r"[+-]?\d+", value) else float(value)
    return value


def read_csv(path: Path):
    with path.open(newline="", encoding="utf-8-sig") as stream:
        rows = list(csv.reader(stream))
    return rows[0], [[typed(value) for value in row] for row in rows[1:]]


def address(row: int, column: int):
    label = ""
    while column >= 0:
        label = chr(65 + column % 26) + label
        column = column // 26 - 1
    return f"{label}{row + 1}"


def read_workbook(path: Path):
    sheets, formulas = {}, {}
    with zipfile.ZipFile(path) as archive:
        strings = []
        if "xl/sharedStrings.xml" in archive.namelist():
            strings = [
                "".join(node.text or "" for node in item.findall(".//m:t", NS))
                for item in ET.fromstring(archive.read("xl/sharedStrings.xml")).findall("m:si", NS)
            ]
        relationships = {
            node.attrib["Id"]: node.attrib["Target"]
            for node in ET.fromstring(archive.read("xl/_rels/workbook.xml.rels"))
        }
        for sheet in ET.fromstring(archive.read("xl/workbook.xml")).findall("m:sheets/m:sheet", NS):
            name = sheet.attrib["name"]
            target = relationships[sheet.attrib["{" + NS["r"] + "}id"]]
            member = target.lstrip("/") if target.startswith("/") else "xl/" + target
            xml = ET.fromstring(archive.read(member))
            values, formulas[name] = {}, []
            for cell in xml.findall("m:sheetData/m:row/m:c", NS):
                ref, kind = cell.attrib["r"], cell.attrib.get("t")
                node = cell.find("m:v", NS)
                raw = node.text if node is not None else None
                if kind == "s":
                    value = strings[int(raw)]
                elif kind == "inlineStr":
                    value = "".join(node.text or "" for node in cell.findall(".//m:t", NS))
                elif kind == "b":
                    value = raw == "1"
                elif kind == "e":
                    raise AssertionError(f"Excel error in {name}!{ref}: {raw}")
                elif kind == "str":
                    value = raw
                else:
                    value = float(raw) if raw is not None else None
                values[ref] = value
                formula = cell.find("m:f", NS)
                if formula is not None:
                    formulas[name].append((ref, formula.text, value))
            sheets[name] = values
    return sheets, formulas


def verify(args):
    sources = {
        "Landmark inventory": args.landmark_inventory,
        "Original peak audit": args.peak_audit,
        "Summary source": args.summary_source,
    }
    datasets = {name: read_csv(path) for name, path in sources.items()}
    sheets, formulas = read_workbook(args.workbook)
    assert list(sheets) == SHEETS, f"Unexpected worksheet structure: {list(sheets)}"
    source_checks = {}
    for name, (headers, rows) in datasets.items():
        expected = [headers, *rows]
        nonempty_addresses = set()
        for row_number, row in enumerate(expected):
            assert len(row) == len(headers), f"CSV row width changed: {name}, row {row_number + 1}"
            for column_number, value in enumerate(row):
                ref = address(row_number, column_number)
                actual = sheets[name].get(ref)
                if value is not None:
                    nonempty_addresses.add(ref)
                if isinstance(value, (int, float)) and not isinstance(value, bool):
                    equal = actual == value
                else:
                    equal = type(actual) is type(value) and actual == value
                assert equal, f"Value mismatch at {name}!{ref}: expected {value!r}, found {actual!r}"
        assert not formulas[name], f"Source sheet contains unexpected formulas: {name}"
        actual_nonempty = {ref for ref, value in sheets[name].items() if value is not None}
        assert actual_nonempty == nonempty_addresses, f"Unexpected populated source cells: {name}"
        source_checks[name] = {
            "rows": len(rows), "columns": len(headers),
            "exact_source_cells": len(expected) * len(headers),
            "sha256": hashlib.sha256(sources[name].read_bytes()).hexdigest(),
        }

    inventory_headers, inventory_rows = datasets["Landmark inventory"]
    inventory = [dict(zip(inventory_headers, row)) for row in inventory_rows]
    summary_headers, summary_rows = datasets["Summary source"]
    summary = [dict(zip(summary_headers, row)) for row in summary_rows]
    for row in summary:
        eligible = [
            item for item in inventory
            if item["eligible"] and item["group"] == row["group"]
            and (
                row["panel"] == "all_landmarks"
                or (item["landmark"] == "peak" if row["panel"] == "peak"
                    else item["landmark"] in ("rise", "fall"))
            )
        ]
        categories = collections.Counter(item["summary_category"] for item in eligible)
        assert categories[row["category"]] == row["count"], f"Summary count mismatch: {row}"
        assert len(eligible) == row["total"], f"Summary denominator mismatch: {row}"
        assert math.isclose(row["percentage"], 100 * row["count"] / row["total"], abs_tol=1e-12)

    group_codes = ["pair-rule", "gap-gap", "gap-pair-rule"]
    category_codes = [row["category"] for row in summary[:5]]
    lookups = {(row["panel"], row["group"], row["category"]): row for row in summary}
    blocks = [("Overall summary", "all_landmarks", 7), ("By landmark", "peak", 7), ("By landmark", "half_height", 18)]
    for sheet_name, panel, first in blocks:
        totals = [lookups[panel, group, category_codes[0]]["total"] for group in group_codes]
        for offset, category in enumerate(category_codes):
            row_number = first + offset
            counts = [lookups[panel, group, category]["count"] for group in group_codes]
            counts.append(sum(counts))
            denominators = totals + [sum(totals)]
            for count_column, percent_column, count, total in zip(("B", "D", "F", "H"), ("C", "E", "G", "I"), counts, denominators):
                assert sheets[sheet_name][f"{count_column}{row_number}"] == count
                assert math.isclose(sheets[sheet_name][f"{percent_column}{row_number}"], count / total, abs_tol=1e-12)
        total_row = first + 5
        for column, total in zip(("B", "D", "F", "H"), totals + [sum(totals)]):
            assert sheets[sheet_name][f"{column}{total_row}"] == total
        for ref, formula, cache in formulas[sheet_name]:
            assert cache is not None, f"Missing formula cache: {sheet_name}!{ref}"

    return {
        "status": "passed",
        "workbook_sha256": hashlib.sha256(args.workbook.read_bytes()).hexdigest(),
        "sheets": list(sheets),
        "source_tables": source_checks,
        "all_source_values_exact": True,
        "source_summaries_match_inventory": True,
        "displayed_counts_and_percentages_match_sources": True,
        "excel_formula_errors": 0,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workbook", type=Path, required=True)
    parser.add_argument("--peak-audit", type=Path, required=True)
    parser.add_argument("--landmark-inventory", type=Path, required=True)
    parser.add_argument("--summary-source", type=Path, required=True)
    parser.add_argument("--report", type=Path)
    args = parser.parse_args()
    report = verify(args)
    output = json.dumps(report, indent=2)
    if args.report:
        args.report.parent.mkdir(parents=True, exist_ok=True)
        args.report.write_text(output + "\n", encoding="utf-8")
    print(output)


if __name__ == "__main__":
    main()
