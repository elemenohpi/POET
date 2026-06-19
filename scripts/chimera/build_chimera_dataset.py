"""Build tokenized POET datasets from the chimera MultiHance workbook."""

import argparse
import os
import re
import sys
from pathlib import Path
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))

import chimera_codec as CC
import sequence_codec as SC


NS = {
    "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
    "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
}


def _col_to_idx(ref):
    match = re.match(r"([A-Z]+)", ref or "")
    col = match.group(1) if match else ""
    value = 0
    for char in col:
        value = value * 26 + ord(char) - 64
    return value - 1


def _shared_strings(zip_file):
    root = ET.fromstring(zip_file.read("xl/sharedStrings.xml"))
    values = []
    for item in root.findall("main:si", NS):
        values.append("".join(t.text or "" for t in item.findall(".//main:t", NS)))
    return values


def _cell_value(cell, shared):
    cell_type = cell.attrib.get("t")
    value = cell.find("main:v", NS)
    if cell_type == "s":
        return shared[int(value.text)] if value is not None and value.text else ""
    if cell_type == "inlineStr":
        return "".join(t.text or "" for t in cell.findall(".//main:t", NS))
    return value.text if value is not None and value.text else ""


def _sheet_paths(zip_file):
    workbook = ET.fromstring(zip_file.read("xl/workbook.xml"))
    rels = ET.fromstring(zip_file.read("xl/_rels/workbook.xml.rels"))
    rid_to_target = {rel.attrib["Id"]: rel.attrib["Target"] for rel in rels}
    paths = {}
    rel_key = "{%s}id" % NS["rel"]
    for sheet in workbook.findall(".//main:sheets/main:sheet", NS):
        target = rid_to_target[sheet.attrib[rel_key]]
        paths[sheet.attrib["name"]] = (
            "xl/" + target.lstrip("/") if not target.startswith("xl/") else target
        )
    return paths


def _rows(zip_file, path, shared):
    root = ET.fromstring(zip_file.read(path))
    for row in root.findall(".//main:sheetData/main:row", NS):
        values = []
        for cell in row.findall("main:c", NS):
            idx = _col_to_idx(cell.attrib.get("r", ""))
            while len(values) <= idx:
                values.append("")
            values[idx] = _cell_value(cell, shared)
        while values and values[-1] == "":
            values.pop()
        if values:
            yield int(row.attrib.get("r", "0")), values


def _float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def _read_multihance_replicates(workbook_path):
    records = []
    with ZipFile(workbook_path) as zip_file:
        shared = _shared_strings(zip_file)
        paths = _sheet_paths(zip_file)

        for row_no, values in _rows(
            zip_file, paths["B2chimeraMultihance"], shared
        ):
            if len(values) < 6:
                continue
            value = _float(values[4])
            if value is None or not values[5]:
                continue
            records.append(
                {
                    "sheet": "B2chimeraMultihance",
                    "row": row_no,
                    "replicate_label": values[3],
                    "construct": values[5],
                    "raw_value": value,
                }
            )

        for row_no, values in _rows(
            zip_file, paths["B3chimeraMultihance"], shared
        ):
            if len(values) < 15:
                continue
            value = _float(values[14])
            if value is None or not values[13]:
                continue
            records.append(
                {
                    "sheet": "B3chimeraMultihance",
                    "row": row_no,
                    "replicate_label": values[12],
                    "construct": values[13],
                    "raw_value": value,
                }
            )
    return pd.DataFrame(records)


def _construct_to_sequence(sheet, construct, region_count):
    construct = str(construct).strip()

    if construct == "B2":
        return CC.encode_construct("B2", region_count=region_count), "native_B2"
    if construct == "B3":
        return CC.encode_construct("B3", region_count=region_count), "native_B3"

    if sheet == "B2chimeraMultihance":
        match = re.fullmatch(r"TM(\d+)", construct)
        if match:
            slot = int(match.group(1))
            swaps = [(slot, "B3", slot)]
            return (
                CC.encode_construct("B2", swaps, region_count),
                "B2_backbone_B3_TM{:02d}".format(slot),
            )

    if sheet == "B3chimeraMultihance":
        match = re.fullmatch(r"T(\d+)", construct)
        if match:
            slot = int(match.group(1))
            swaps = [(slot, "B2", slot)]
            return (
                CC.encode_construct("B3", swaps, region_count),
                "B3_backbone_B2_TM{:02d}".format(slot),
            )

    return None, "skipped"


def build_dataset(
    workbook_path,
    output_data,
    alphabet_data,
    manifest_data,
    raw_replicates,
    skipped_groups,
    region_count=CC.DEFAULT_REGION_COUNT,
    normalize="b2_ratio",
):
    replicates = _read_multihance_replicates(workbook_path)
    raw_replicates.parent.mkdir(parents=True, exist_ok=True)
    replicates.to_csv(raw_replicates, index=False)

    summary = (
        replicates.groupby(["sheet", "construct"], as_index=False)
        .agg(mean_raw=("raw_value", "mean"), sd_raw=("raw_value", "std"), n=("raw_value", "count"))
        .fillna({"sd_raw": 0.0})
    )

    b2_baseline = {}
    for sheet, rows in summary.groupby("sheet"):
        b2 = rows[rows["construct"] == "B2"]
        if not b2.empty:
            b2_baseline[sheet] = float(b2.iloc[0]["mean_raw"])

    data_rows = []
    skipped = []
    for _, row in summary.iterrows():
        tokens, construct_type = _construct_to_sequence(
            row["sheet"], row["construct"], region_count
        )
        if tokens is None:
            skipped.append(row.to_dict())
            continue

        if normalize == "b2_ratio":
            baseline = b2_baseline.get(row["sheet"])
            if not baseline:
                skipped.append(row.to_dict())
                continue
            fitness = float(row["mean_raw"]) / baseline
        elif normalize == "none":
            fitness = float(row["mean_raw"])
        else:
            raise ValueError("Unknown normalize mode '{}'".format(normalize))

        data_rows.append(
            {
                "sequence": SC.join_tokens(tokens, {"sequence_mode": "token"}),
                "fitness": fitness,
                "sheet": row["sheet"],
                "construct": row["construct"],
                "construct_type": construct_type,
                "mean_raw": row["mean_raw"],
                "sd_raw": row["sd_raw"],
                "n": row["n"],
                "normalization": normalize,
            }
        )

    output_data.parent.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(data_rows).to_csv(output_data, index=False)
    pd.DataFrame(skipped).to_csv(skipped_groups, index=False)

    manifest = pd.DataFrame(CC.manifest_rows(region_count))
    manifest_data.parent.mkdir(parents=True, exist_ok=True)
    manifest.to_csv(manifest_data, index=False)
    manifest[["code"]].to_csv(alphabet_data, index=False)

    return len(data_rows), len(skipped)


def main():
    parser = argparse.ArgumentParser(
        description="Build tokenized chimera MultiHance POET training data."
    )
    parser.add_argument(
        "--workbook",
        default="chimera_init_data/ChimeraICPdata.xlsx",
        help="Input workbook path.",
    )
    parser.add_argument(
        "--output",
        default="data/chimera_multihance_relative.csv",
        help="POET training CSV output.",
    )
    parser.add_argument(
        "--alphabet",
        default="chimera_init_data/chimera_token_alphabet.csv",
        help="Token alphabet CSV output.",
    )
    parser.add_argument(
        "--manifest",
        default="chimera_init_data/chimera_token_manifest.csv",
        help="Token manifest CSV output.",
    )
    parser.add_argument(
        "--raw-replicates",
        default="data/chimera_multihance_raw_replicates.csv",
        help="Raw replicate table output.",
    )
    parser.add_argument(
        "--skipped",
        default="chimera_init_data/chimera_skipped_groups.csv",
        help="Skipped construct summary output.",
    )
    parser.add_argument(
        "--normalize",
        choices=["b2_ratio", "none"],
        default="b2_ratio",
        help="Fitness scaling mode.",
    )
    parser.add_argument("--regions", type=int, default=CC.DEFAULT_REGION_COUNT)
    args = parser.parse_args()

    kept, skipped = build_dataset(
        Path(args.workbook),
        Path(args.output),
        Path(args.alphabet),
        Path(args.manifest),
        Path(args.raw_replicates),
        Path(args.skipped),
        args.regions,
        args.normalize,
    )
    print("Wrote {} training rows; skipped {} groups.".format(kept, skipped))


if __name__ == "__main__":
    main()
