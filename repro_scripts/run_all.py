#!/usr/bin/env python3
"""
run_all.py
==========
Unified orchestrator for the 16S rRNA Amplicon Dual-Track Reproduction Pipeline
(Bio-Repro Unified Pipeline Blueprint v1.0).

Modes:
  - smoke : 10-second deterministic compositionality & decontamination smoke test
  - quick : Skip raw FASTQ DADA2 denoising & tree building; reuse 01_ps_miseq.rds
  - full  : End-to-end dual-track execution (00_smoke -> 01_dada2 -> 02_diversity -> 03_daa -> 04_network_func)
"""

from __future__ import annotations

import argparse
import csv
import pathlib
import shutil
import subprocess
import sys
import time
from typing import Dict, List, Sequence, Tuple


def get_project_root() -> pathlib.Path:
    """Resolve project root dynamically regardless of current working directory."""
    here = pathlib.Path(__file__).resolve().parent
    return here.parent if here.name == "repro_scripts" else here


def resolve_rscript() -> str:
    """Locate Rscript on PATH or common Windows R installation directories."""
    found = shutil.which("Rscript")
    if found:
        return found
    return "Rscript"


def get_mode_steps(mode: str) -> List[Tuple[str, List[str]]]:
    """Return ordered subprocess execution list for the selected pipeline mode."""
    rscript = resolve_rscript()
    if mode == "smoke":
        return [
            ("00_smoke_test", [rscript, "repro_scripts/00_smoke_test.R"]),
        ]
    if mode == "quick":
        return [
            ("02_decontam_and_diversity", [rscript, "repro_scripts/02_decontam_and_diversity.R"]),
            ("03_daa_ancombc2_maaslin3_lefse", [rscript, "repro_scripts/03_daa_ancombc2_maaslin3_lefse.R"]),
            ("04_microeco_network_and_func", [rscript, "repro_scripts/04_microeco_network_and_func.R"]),
        ]
    return [
        ("00_smoke_test", [rscript, "repro_scripts/00_smoke_test.R"]),
        ("01_dada2_fastq_to_phyloseq", [rscript, "repro_scripts/01_dada2_fastq_to_phyloseq.R"]),
        ("02_decontam_and_diversity", [rscript, "repro_scripts/02_decontam_and_diversity.R"]),
        ("03_daa_ancombc2_maaslin3_lefse", [rscript, "repro_scripts/03_daa_ancombc2_maaslin3_lefse.R"]),
        ("04_microeco_network_and_func", [rscript, "repro_scripts/04_microeco_network_and_func.R"]),
    ]


def check_preflight_inputs(root: pathlib.Path, required_relpaths: Sequence[str]) -> None:
    """Fail early with actionable hints if gitignored checkpoint inputs are missing."""
    missing = [p for p in required_relpaths if not (root / p).exists()]
    if missing:
        joined = "\n  - ".join(missing)
        raise FileNotFoundError(
            f"Missing gitignored checkpoint input(s) for selected mode:\n  - {joined}\n"
            "Run `--mode full` first or see data/README.md for exact acquisition details."
        )


def run_step(name: str, cmd: List[str], root: pathlib.Path, dry_run: bool) -> float:
    """Execute one pipeline step in a child process and persist stdout/stderr to logs/."""
    print(f">>> [{name}] CMD: {' '.join(cmd)}")
    if dry_run:
        print(f"<<< [{name}] DRY-RUN OK\n")
        return 0.0

    log_dir = root / "logs"
    log_dir.mkdir(parents=True, exist_ok=True)
    log_path = log_dir / f"{name}.log"

    t0 = time.perf_counter()
    proc = subprocess.run(
        cmd,
        cwd=str(root),
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
    )
    elapsed = round(time.perf_counter() - t0, 2)

    lines = [
        f"=== STEP: {name} | EXIT={proc.returncode} | ELAPSED={elapsed}s ===",
        f"CMD: {' '.join(cmd)}",
        proc.stdout.rstrip(),
    ]
    if proc.stderr.strip():
        lines.append("[STDERR]\n" + proc.stderr.rstrip())
    log_path.write_text("\n".join(lines) + "\n", encoding="utf-8")

    if proc.stdout.strip():
        print(proc.stdout.rstrip())
    if proc.stderr.strip():
        print(proc.stderr.rstrip(), file=sys.stderr)

    if proc.returncode != 0:
        raise RuntimeError(f"Step '{name}' failed with exit code {proc.returncode} (see {log_path})")

    print(f"<<< [{name}] SUCCESS in {elapsed}s (log: logs/{name}.log)\n")
    return elapsed


def ensure_python_in_versions_tsv(root: pathlib.Path) -> pathlib.Path:
    """Ensure Python runtime version is recorded alongside child R package versions."""
    out_tsv = root / "repro_results" / "environment_versions.tsv"
    out_tsv.parent.mkdir(parents=True, exist_ok=True)
    py_ver = sys.version.split()[0]

    existing_rows: List[Dict[str, str]] = []
    has_python = False
    if out_tsv.exists():
        with out_tsv.open("r", encoding="utf-8", errors="replace", newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                if row.get("Package") == "Python":
                    row["Version"] = py_ver
                    row["Repository"] = "Host Runtime (Stdlib)"
                    has_python = True
                existing_rows.append(row)

    if not has_python:
        existing_rows.insert(0, {"Package": "Python", "Version": py_ver, "Repository": "Host Runtime (Stdlib)"})

    with out_tsv.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.DictWriter(fh, fieldnames=["Package", "Version", "Repository"], delimiter="\t")
        writer.writeheader()
        writer.writerows(existing_rows)
    return out_tsv


def main() -> int:
    parser = argparse.ArgumentParser(description="16S Amplicon Dual-Track Pipeline Orchestrator")
    parser.add_argument("--mode", choices=["smoke", "quick", "full"], default="smoke")
    parser.add_argument("--dry-run", action="store_true", help="Print resolved commands without executing")
    parser.add_argument("--skip-verify", action="store_true", help="Skip Ground-Truth verification gate")
    args = parser.parse_args()

    root = get_project_root()
    print("===================================================================")
    print("  16S rRNA Amplicon Dual-Track Pipeline Orchestrator (Blueprint v1.0)")
    print(f"  Workspace : {root}")
    print(f"  Mode      : {args.mode.upper()} (dry_run={args.dry_run})")
    print("===================================================================\n")

    preflight_map: Dict[str, List[str]] = {
        "smoke": [],
        "quick": ["data/processed/01_ps_miseq.rds"],
        "full": [],
    }
    if not args.dry_run:
        check_preflight_inputs(root, preflight_map[args.mode])

    steps = get_mode_steps(args.mode)
    timings: List[Tuple[str, float]] = []
    for name, cmd in steps:
        elapsed = run_step(name, cmd, root, dry_run=args.dry_run)
        timings.append((name, elapsed))

    if not args.dry_run:
        env_tsv = ensure_python_in_versions_tsv(root)
        print(f"Environment snapshot verified: {env_tsv}")

        verifier = root / "repro_scripts" / "verify_outputs.py"
        contract = root / "ground_truth_contract.json"
        if not args.skip_verify and verifier.exists() and contract.exists():
            run_step(
                "verify_outputs",
                [
                    sys.executable,
                    str(verifier),
                    "--root",
                    str(root),
                    "--contract",
                    str(contract),
                    "--mode",
                    args.mode,
                    "--out",
                    "repro_results/verification_report.json",
                ],
                root,
                dry_run=False,
            )

    print("=== Pipeline Execution Summary ===")
    for name, sec in timings:
        print(f"  - {name}: {sec}s")
    return 0


if __name__ == "__main__":
    sys.exit(main())
