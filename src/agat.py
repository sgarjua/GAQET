"""
agat.py
=======
Wrapper helper for a single external tool:

* **AGAT** - computes basic statistics on a genome annotation (GFF/GTF).

Each runner returns a dictionary with the executed command, an informational message,
the main output path, and a ``returncode`` (99 means “already done”).
"""

import subprocess
from pathlib import Path
from typing import Dict, Any


# ---------------------------------------------------------------------------
# 1. Run AGAT statistics (extracts metrics from a GFF/GTF file)
# ---------------------------------------------------------------------------
def run_agat(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Run AGAT statistics (or skip if already done)."""

    # Output directory and file
    outdir = arguments["output"] / "GenomeAnnStats"
    outdir.mkdir(parents=True, exist_ok=True)
    out_fpath = outdir / "ResultAgat.txt"

    # AGAT command
    cmd_list = ["agat_sp_statistics.pl", "--gff", "{}".format(arguments["annotation"]), "-o", "{}".format(out_fpath)]
    command = ' '.join(cmd_list)

# Check if AGAT is already done
    if out_fpath.exists():
        return {"command": command,
                "msg": "AGAT already done",
                "out_fpath": out_fpath, 
                "returncode": 99}
    else:
        # Run AGAT
        run_ = subprocess.run(command, shell=True, stdout=subprocess.PIPE)

        if run_.returncode == 0:
            msg = "AGAT run successfully"
        else:
            msg = "AGAT Failed: \n {}".format(run_.stdout)

        return {"command": command,
                "msg": msg,
                "out_fpath": out_fpath,
                "returncode": run_.returncode}


# ---------------------------------------------------------------------------
# 2. Parse the AGAT statistics file and keep selected metrics
# ---------------------------------------------------------------------------
def get_agat_stats(agat_statistics: Dict[str, Any]) -> Dict[str, Any]:
    """Parse the AGAT result file into a dict of metrics."""

    # Pre‑fill all expected keys with 0
    results = {
        "Gene_Models (N)": 0,
        "Transcript_Models (N)": 0,
        "CDS_Models (N)": 0,
        "Exons (N)": 0,
        "UTR5' (N)": 0,
        "UTR3' (N)": 0,
        "Overlapping_Gene_Models (N)": 0,
        "Single Exon Gene Models (N)": 0,
        "Single Exon Transcripts (N)": 0,
        "Total Gene Space (Mb)": 0,
        "Mean Gene Model Length (bp)": 0,
        "Mean CDS Model Length (bp)": 0,
        "Mean Exon Length (bp)": 0,
        "Mean Intron Length (bp)": 0,
        "Longest Gene Model Length (bp)": 0,
        "Longest CDS Model Length (bp)": 0,
        "Longest Intron Length (bp)": 0,
        "Shortest Gene Model Length (bp)": 0,
        "Shortest CDS Model Length (bp)": 0,
        "Shortest Intron Length (bp)": 0
    }
    
    # Mapping dictionaries (unchanged from original version)
    mapping_mrna = {
        "Number of gene": "Gene_Models (N)",
        "Number of mrna": "Transcript_Models (N)",
        "Number of cds": "CDS_Models (N)",
        "Number of exon": "Exons (N)",
        "Number of five_prime_utr": "UTR5' (N)",
        "Number of three_prime_utr": "UTR3' (N)",
        "Number gene overlapping": "Overlapping_Gene_Models (N)",
        "Number of single exon gene": "Single Exon Gene Models (N)",
        "Number of single exon mrna": "Single Exon Transcripts (N)",
        "Total gene length (bp)": "Total Gene Space (Mb)",
        "mean gene length (bp)": "Mean Gene Model Length (bp)",
        "mean cds length (bp)": "Mean CDS Model Length (bp)",
        "mean exon length (bp)": "Mean Exon Length (bp)",
        "mean intron in cds length (bp)": "Mean Intron Length (bp)",
        "Longest gene (bp)": "Longest Gene Model Length (bp)",
        "Longest cds (bp)": "Longest CDS Model Length (bp)",
        "Longest intron into cds part (bp)": "Longest Intron Length (bp)",
        "Shortest gene (bp)": "Shortest Gene Model Length (bp)",
        "Shortest cds piece (bp)": "Shortest CDS Model Length (bp)",
        "Shortest intron into cds part (bp)": "Shortest Intron Length (bp)"
    }
    mapping_transcript = {
        "Number of gene": "Gene_Models (N)",
        "Number of transcript": "Transcript_Models (N)",
        # "Number of cds": "CDS_Models (N)",
        "Number of exon": "Exons (N)",
        # "Number of five_prime_utr": "UTR5' (N)",
        # "Number of three_prime_utr": "UTR3' (N)",
        "Number gene overlapping": "Overlapping_Gene_Models (N)",
        "Number of single exon gene": "Single Exon Gene Models (N)",
        "Number of single exon transcript": "Single Exon Transcripts (N)",
        "Total gene length (bp)": "Total Gene Space (Mb)",
        "mean gene length (bp)": "Mean Gene Model Length (bp)",
        # "mean cds length (bp)": "Mean CDS Model Length (bp)",
        "mean exon length (bp)": "Mean Exon Length (bp)",
        # "mean intron in cds length (bp)": "Mean Intron Length (bp)",
        "Longest gene (bp)": "Longest Gene Model Length (bp)",
        # "Longest cds (bp)": "Longest CDS Model Length (bp)",
        # "Longest intron into cds part (bp)": "Longest Intron Length (bp)",
        "Shortest gene (bp)": "Shortest Gene Model Length (bp)",
        # "Shortest cds piece (bp)": "Shortest CDS Model Length (bp)",
        "Shortest intron into exon part (bp)": "Shortest Intron Length (bp)"
    }
    
    # Read the statistics file produced by AGAT
    with open(agat_statistics["out_fpath"], 'r') as stats_fhand:
        for line in stats_fhand:
            # Switch mapping depending on the section header
            if "--- transcript ---" in line:
                mapping = mapping_transcript
            elif "--- mrna ---" in line:
                mapping = mapping_mrna

            if not line.rstrip():   # empty line → skip it
                continue
            if ':' in line:         # end-of-block marker → stop reading
                break

            # Save metrics we care about
            try:
                key, val = line.rsplit(maxsplit=1)
                key = key.strip()
                val = int(val.strip())
                if key in mapping:
                    result_key = mapping[key]
                    if result_key == "Total Gene Space (Mb)":
                        results[result_key] = round(val / 1_000_000, 2) # convert bp → Mb
                    else:
                        results[result_key] = val
            except ValueError:
                continue           # skip unparsable line
    return results