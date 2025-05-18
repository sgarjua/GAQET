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
import time
from csv import DictReader
from pathlib import Path
from datetime import datetime
import os

# === Project-specific imports ===
from src.agat import run_agat, get_agat_stats
from src.busco import run_busco, run_gffread, get_busco_results
from src.LTR_retriever import (
    create_outdir, run_suffixerator, run_harvest, run_finder,
    concatenate_outputs, run_LTR_retriever, run_LAI, get_LAI
)
from src.stringtie import run_stringtie, run_gffcompare, calculate_annotation_scores
from src.table import AGAT_COLS, RNASEQ_COLS


def write_time_in_file(file: str, text: str):
    with open(file, 'a') as f:
        f.write(text+"\n")


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

    name_time_file: str = "time_file.txt"
    route_time_file = out_dir / name_time_file

    # For each sample: create a folder, run the 4 pipelines and save results in "stats"
    for name, values in arguments["input"].items():
        write_time_in_file(route_time_file, name)
        stats[name] = {}
        name_dir = out_dir / name
        values["output"] = name_dir
        values["threads"] = arguments["threads"]
        if not name_dir.exists():
            name_dir.mkdir(parents=True, exist_ok=True)
        start = time.ctime()

        # # AGAT
        # start_time = time.time()
        # agat_statistics = run_agat(values)
        # end_time = time.time()
        # print(agat_statistics)
        # write_time_in_file(route_time_file, "   Time consumed by GenomeAnnStats: {}s\n\n".format(round(end_time-start_time, 2)))
        # stats[name]["agat_statistics"] = get_agat_stats(agat_statistics)

        # # BUSCO
        # start_time = time.time()
        # gffread_results = run_gffread(values)
        # print(gffread_results)
        # busco_results = run_busco(values)
        # print(busco_results)
        # end_time = time.time()
        # write_time_in_file(route_time_file, "   Time consumed by BUSCOCompleteness: {}s\n\n".format(round(end_time-start_time, 2)))
        # print("\nTime consumed by BUSCOCompleteness: {}s\n\n".format(round(end_time-start_time, 2)))
        # stats[name]["busco_results"] = get_busco_results(busco_results, lineage=values["lineage"])

        # # LAI
        # start_time = time.time()
        # LAI_out_dir =  create_outdir(values)
        # print(LAI_out_dir)
        # values["LAI_dir"] = LAI_out_dir["out_fpath"]
        # suffixerator =  run_suffixerator(values)
        # if "returncode" in suffixerator:
        #     if suffixerator["returncode"] == 1:
        #         raise RuntimeError("Suffixerator has failed")
        # print(suffixerator)
        # harvest = run_harvest(values)
        # print(harvest)
        # finder = run_finder(values)
        # print(finder)
        # cat = concatenate_outputs(values)
        # print(cat)
        # LTR = run_LTR_retriever(values)
        # print(LTR)
        # LAI = run_LAI(values)
        # print(LAI)
        # end_time = time.time()
        # write_time_in_file(route_time_file, "   Time consumed by LAICompleteness: {}s\n\n".format(round(end_time-start_time, 2)))
        # print("\nTime consumed by LAICompleteness: {}s\n\n".format(round(end_time-start_time, 2)))
        # stats[name]["LAI"] = get_LAI(LAI)

        dir_bam: str = values["alignments"]
        list_bam_files = [file for file in dir_bam if file.endswith(".bam")]
        start_time = time.time()

        for bam in list_bam_files:
            values["alignments"] = dir_bam + bam
            # RNA-seq support
            stringtie = run_stringtie(values)
            print(stringtie)
            gffcompare = run_gffcompare(values)
            print(gffcompare)
            annotation_scores = calculate_annotation_scores(values)
            print(annotation_scores)
            stats[name][bam]["annotation_scores"] = annotation_scores
        end_time = time.time()
        print("\nTime consumed by RNASeqCheck: {}s\n\n".format(round(end_time-start_time, 2)))
        write_time_in_file(route_time_file, "   Time consumed by RNASeqCheck: {}s\n\n".format(round(end_time-start_time, 2)))

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
            # results += [stats[name]["agat_statistics"][stat] for stat in AGAT_COLS]
            # results += [stats[name]["busco_results"]]
            # results += [stats[name]["LAI"]]
            results += [stats[name][bam]["annotation_scores"][score] for score in RNASEQ_COLS for bam in stats[name].keys()]
            summary.write("\t".join(results)+"\n")
        
if __name__ == "__main__":
    main()