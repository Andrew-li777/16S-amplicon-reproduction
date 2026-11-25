#!/usr/bin/env python3
"""
verify_outputs.py
=================
Zero-dependency Ground-Truth verification gate for reproduction/16S
(Bio-Repro Unified Pipeline Blueprint v1.0).

Validates output artifacts declared in `ground_truth_contract.json`:
  - File existence & minimum byte size
  - Minimum data row count for TSV tables
  - Exact string or numeric [min, max] assertions on key-value TSV (`tsv_kv`)
    or `Key: Value` TXT (`txt_kv`) summaries

Usage:
  python repro_scripts/verify_outputs.py --mode smoke
  python repro_scripts/verify_outputs.py --contract ground_truth_contract.json --mode full --out repro_results/verification_report.json
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys
from typing import Any, Dict, List, Tuple


def resolve_project_root() -> pathlib.Path:
    """Resolve project root dynamically regardless of current working directory."""
    here = pathlib.Path(__file__).resolve().parent
    return here.parent if here.name == "repro_scripts" else here


def parse_tsv_kv(path: pathlib.Path) -> Tuple[int, Dict[str, str]]:
    """Parse a 2-column key-value TSV with header row."""
    kv: Dict[str, str] = {}
    rows = 0
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        header = next(reader, None)
        if header is None:
            return 0, kv
        for row in reader:
            if not row or all(not c.strip() for c in row):
                continue
            rows += 1
            if len(row) >= 2:
                kv[row[0].strip()] = row[1].strip()
    return rows, kv


def parse_tsv_row_count(path: pathlib.Path) -> int:
    """Count non-empty data rows (excluding header) in a TSV table."""
    rows = 0
    with path.open("r", encoding="utf-8", errors="replace", newline="") as fh:
        reader = csv.reader(fh, delimiter="\t")
        next(reader, None)
        for row in reader:
            if row and any(c.strip() for c in row):
                rows += 1
    return rows


def parse_txt_kv(path: pathlib.Path) -> Dict[str, str]:
    """Parse `Key: Value` lines from plain-text verification summaries."""
    kv: Dict[str, str] = {}
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if ":" in line:
            k, v = line.split(":", 1)
            k_clean = k.strip()
            v_clean = v.strip()
            if k_clean:
                kv[k_clean] = v_clean
    return kv


def check_assertion(kv: Dict[str, str], rule: Dict[str, Any]) -> Tuple[bool, str]:
    """Evaluate a single key assertion (equals / min / max)."""
    key = str(rule["key"])
    if key not in kv:
        return False, f"Missing key '{key}'"
    raw_val = kv[key]
    if "equals" in rule:
        expected = str(rule["equals"])
        ok = raw_val == expected
        return ok, f"{key} == '{expected}' (got '{raw_val}')"
    try:
        num = float(raw_val)
    except ValueError:
        return False, f"Key '{key}' expected numeric value, got '{raw_val}'"
    if "min" in rule and num < float(rule["min"]):
        return False, f"{key}={num} < min {rule['min']}"
    if "max" in rule and num > float(rule["max"]):
        return False, f"{key}={num} > max {rule['max']}"
    return True, f"{key}={num} within bounds"


def verify_contract(root: pathlib.Path, contract: Dict[str, Any], mode: str) -> Dict[str, Any]:
    """Run all file and metric checks for a given execution mode."""
    modes_cfg = contract.get("modes", {})
    if mode not in modes_cfg:
        raise ValueError(f"Mode '{mode}' not defined in contract (available: {list(modes_cfg.keys())})")

    checks: List[Dict[str, Any]] = []
    failures: List[str] = []

    for spec in modes_cfg[mode].get("required_files", []):
        rel_path = spec["path"]
        target = root / rel_path
        item_res: Dict[str, Any] = {"path": rel_path, "status": "PASS", "details": []}

        if not target.exists():
            item_res["status"] = "FAIL"
            msg = f"Missing required file: {rel_path}"
            item_res["details"].append(msg)
            failures.append(msg)
            checks.append(item_res)
            continue

        size = target.stat().st_size
        min_bytes = int(spec.get("min_bytes", 1))
        if size < min_bytes:
            item_res["status"] = "FAIL"
            msg = f"File too small ({size} B < {min_bytes} B): {rel_path}"
            item_res["details"].append(msg)
            failures.append(msg)

        fmt = spec.get("format")
        kv: Dict[str, str] = {}
        data_rows: int | None = None

        if fmt == "tsv_kv":
            data_rows, kv = parse_tsv_kv(target)
        elif fmt == "tsv_table":
            data_rows = parse_tsv_row_count(target)
        elif fmt == "txt_kv":
            kv = parse_txt_kv(target)

        if "min_data_rows" in spec and data_rows is not None:
            min_rows = int(spec["min_data_rows"])
            if data_rows < min_rows:
                item_res["status"] = "FAIL"
                msg = f"Row count {data_rows} < {min_rows} in {rel_path}"
                item_res["details"].append(msg)
                failures.append(msg)

        for rule in spec.get("assertions", []):
            ok, detail = check_assertion(kv, rule)
            item_res["details"].append(detail)
            if not ok:
                item_res["status"] = "FAIL"
                failures.append(f"{rel_path}: {detail}")

        checks.append(item_res)

    return {
        "project": contract.get("project_name", root.name),
        "mode": mode,
        "status": "PASS" if not failures else "FAIL",
        "total_files": len(checks),
        "failure_count": len(failures),
        "failures": failures,
        "checks": checks,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="16S Amplicon Ground-Truth Verification Gate")
    parser.add_argument("--root", type=pathlib.Path, default=None, help="Project root directory")
    parser.add_argument(
        "--contract",
        type=pathlib.Path,
        default=pathlib.Path("ground_truth_contract.json"),
        help="Path to JSON contract (relative to root if not absolute)",
    )
    parser.add_argument("--mode", required=True, help="Execution mode to verify (smoke / quick / full)")
    parser.add_argument("--out", type=pathlib.Path, default=None, help="Optional JSON report output path")
    args = parser.parse_args()

    root = args.root.resolve() if args.root else resolve_project_root()
    contract_path = args.contract if args.contract.is_absolute() else (root / args.contract)
    if not contract_path.exists():
        print(f"[ERROR] Contract file not found: {contract_path}", file=sys.stderr)
        return 2

    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    report = verify_contract(root, contract, args.mode)

    if args.out:
        out_path = args.out if args.out.is_absolute() else (root / args.out)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"[VERIFY] Project={report['project']} | Mode={report['mode']} | Status={report['status']} | Files={report['total_files']}")
    for fail in report["failures"]:
        print(f"  [FAIL] {fail}", file=sys.stderr)
    return 0 if report["status"] == "PASS" else 1


if __name__ == "__main__":
    sys.exit(main())
