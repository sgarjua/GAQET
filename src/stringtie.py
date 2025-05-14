"""
stringtie.py
============
Helpers that run **StringTie** and **GFFcompare** on RNA-seq alignments
and derive simple F1-like support scores for each genome annotation.
"""


import subprocess
from pathlib import Path
from typing import Any, Dict

# ---------------------------------------------------------------------------
# 1.  Assemble transcripts with StringTie
# ---------------------------------------------------------------------------
def run_stringtie(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Run *stringtie* on the BAM file (or skip if the GTF already exists)."""
    
    # Output directory
    outdir = arguments["output"] / "RNASeqCheck"
    if not outdir.exists():
        outdir.mkdir(parents=True, exist_ok=True)
    
    # sample‑specific prefix for outputs
    output_name = outdir / Path(arguments["alignments"]).stem
    
    # StringTie command line
    cmd = "stringtie -o {}.gtf -p {} {}".format(output_name,
                                                arguments["threads"],
                                                arguments["alignments"])

    outfile = outdir / "{}.gtf".format(Path(arguments["alignments"]).stem)
    if outfile.exists():
        return {"command": cmd, 
                "msg": "stringtie already done",
                "out_fpath": outdir,
                "returncode": 99}
    else:
        # Run stringtie
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        if run_.returncode == 0:
            msg = "stringtie ran successfully"
        else:
            msg = "stringtie Failed: \n {}".format(run_.stderr)
        return {"command": cmd,
                "msg": msg,
                "out_fpath": outdir,
                "returncode": run_.returncode}


# ---------------------------------------------------------------------------
# 2.  Compare transcripts with GFFcompare
# ---------------------------------------------------------------------------
def run_gffcompare(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Run *gffcompare* against the reference annotation (or skip if done)."""
    
    outdir = arguments["output"] / "RNASeqCheck"
    gtffile = outdir / "{}.gtf".format(Path(arguments["alignments"]).stem)
    output_name = outdir / Path(arguments["alignments"]).stem
    
    # GFFcompare command line
    cmd = "gffcompare -r {} {} -o {}.stats".format(arguments["ref_annotation"],
                                            gtffile,
                                            output_name)
    
    outfile = outdir / "{}.stats".format(Path(arguments["alignments"]).stem)
    if outfile.exists():
        return {"command": cmd, 
                "msg": "gffcompare already done",
                "out_fpath": outfile,
                "returncode": 99}
    else:
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        if run_.returncode == 0:
            msg = "gffcompare ran successfully"
        else:
            msg = "gffcompare Failed: \n {}".format(run_.stderr)
        return {"command": cmd, 
                "msg": msg,
                "out_fpath": outfile, 
                "returncode": run_.returncode}


# ---------------------------------------------------------------------------
# 3.  Derive simple F1-style support scores from *.stats*
# ---------------------------------------------------------------------------
def calculate_annotation_scores(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Return F1 scores and matching counts parsed from the *.stats* file."""
    
    # Locate the stats file produced by GFFcompare
    statsfile = arguments["output"] / "RNASeqCheck"/ "{}.stats".format(Path(arguments["alignments"]).stem)

    # Headings we look for in the .stats report
    f1_checks = ["Transcript level:", "Locus level:"]
    number_check = ["Matching transcripts:", "Matching loci:"]
    
    f1_scores = {} # dictionary for F1 scores
   
    with open(statsfile) as stats_fhand:
        for line in stats_fhand:
            # Extract F1‑like scores
            for check in f1_checks:
                if check in line:
                    line = line.strip()
                    line = line.split()
                    sensitivity = float(line[2])
                    precision = float(line[4])
                    f1_calc = 2*(sensitivity*precision)/(sensitivity+precision)
                    f1_scores[check[:-1]+"_f1"] = f1_calc
            # Matching numbers
            for check in number_check:
                if check in line:
                    line = line.strip()
                    line = line.split()
                    matching_number = line[-1]
                    f1_scores[check] = matching_number
    return f1_scores

