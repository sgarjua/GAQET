#!/usr/bin/env python
"""GAQET.py
Simple command-line tool that runs four external pipelines (AGAT, BUSCO,
LTR_retriever / LAI and RNA-seq with StringTie + GFFcompare) for each sample
listed in a *file-of-files* (FOF) and outputs a single tab-separated summary.

Usage
-----
GAQET.py -i samples.fof -o results/ -t 8
"""

# === Standard library imports ===
import argparse
import sys
from csv import DictReader
from pathlib import Path

# === Project-specific imports ===
from src.agat import run_agat, get_agat_stats
from src.busco import run_busco, run_gffread, get_busco_results
from src.LTR_retriever import (
    create_outdir, run_suffixerator, run_harvest, run_finder,
    concatenate_outputs, run_LTR_retriever, run_LAI, get_LAI
)
from src.stringtie import run_stringtie, run_gffcompare, calculate_annotation_scores
from src.table import AGAT_COLS, RNASEQ_COLS


# ---------------------------------------------------------------------------
# CLI helpers
# ---------------------------------------------------------------------------
def parse_arguments() -> argparse.Namespace:
    """Return parsed command-line arguments."""
    description = "Run quality metrics on genome annotations."
    parser = argparse.ArgumentParser(description=description)

    parser.add_argument("-i", "--input",  required=True, help="Input FOF")
    parser.add_argument("-o", "--output", required=True, help="Output folder")
    parser.add_argument("-t", "--threads", type=int, default=1,
                        help="Threads to use (default 1)")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit()

    return parser.parse_args()

# ---------------------------------------------------------------------------
# Load and validate input
# ---------------------------------------------------------------------------
def get_arguments():
    """Load the FOF file and return a plain dict with paths and settings."""
    parser = parse_arguments()
    fof_fpath = Path(parser.input)
    if not fof_fpath.exists():
        raise RuntimeError("FOF does not exist")
    else:
        samples = {}
        with open(fof_fpath) as fof_fhand:
            samples = {line["name"]: line for line in DictReader(fof_fhand, delimiter="\t")}
                       
    return {"input": samples,
            "threads": parser.threads,
            "output": Path(parser.output)}

# ---------------------------------------------------------------------------
# Main workflow
# ---------------------------------------------------------------------------
def main():
    """Run all analyses and write `summary.tsv`."""
    arguments = get_arguments()
    # Output directory
    out_dir =  arguments["output"]
    if not out_dir.exists():
        out_dir.mkdir(parents=True, exist_ok=True)
    
    # Dictionary for results  
    stats = {}   

    # For each sample: create a folder, run the 4 pipelines and save results in "stats"
    for name, values in arguments["input"].items():
        stats[name] = {}
        name_dir = out_dir / name
        values["output"] = name_dir
        values["threads"] = arguments["threads"]
        if not name_dir.exists():
            name_dir.mkdir(parents=True, exist_ok=True)

        # AGAT
        agat_statistics = run_agat(values)
        print(agat_statistics)
        stats[name]["agat_statistics"] = get_agat_stats(agat_statistics)

        # BUSCO
        gffread_results = run_gffread(values)
        print(gffread_results)
        busco_results = run_busco(values)
        print(busco_results)
        stats[name]["busco_results"] = get_busco_results(busco_results, lineage=values["lineage"])

        # LAI
        LAI_out_dir =  create_outdir(values)
        print(LAI_out_dir)
        values["LAI_dir"] = LAI_out_dir["out_fpath"]
        suffixerator =  run_suffixerator(values)
        if "returncode" in suffixerator:
            if suffixerator["returncode"] == 1:
                raise RuntimeError("Suffixerator has failed")
        print(suffixerator)
        harvest = run_harvest(values)
        print(harvest)
        finder = run_finder(values)
        print(finder)
        cat = concatenate_outputs(values)
        print(cat)
        LTR = run_LTR_retriever(values)
        print(LTR)
        LAI = run_LAI(values)
        print(LAI)
        stats[name]["LAI"] = get_LAI(LAI)

        # RNA-seq support
        stringtie = run_stringtie(values)
        print(stringtie)
        gffcompare = run_gffcompare(values)
        print(gffcompare)
        annotation_scores = calculate_annotation_scores(values)
        print(annotation_scores)
        stats[name]["annotation_scores"] = annotation_scores

    # Write summary as a table
    with open("summary.tsv", "w") as summary:
        header = ["Name"]
        header += AGAT_COLS
        header += ["Busco results"]
        header += ["LAI"]
        header += RNASEQ_COLS
        summary.write("\t".join(header)+"\n")
        
        for name in arguments["input"]:
            results = [name]
            results += [stats[name]["agat_statistics"][stat] for stat in AGAT_COLS]
            results += [stats[name]["busco_results"]]
            results += [stats[name]["LAI"]]
            results += [stats[name]["annotation_scores"][score] for score in RNASEQ_COLS]
            summary.write("\t".join(results)+"\n")
        
if __name__ == "__main__":
    main()