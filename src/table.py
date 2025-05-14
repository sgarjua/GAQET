"""
table.py
========
Defines the column orders for the final summary TSV:

- AGAT_COLS: metrics produced by AGAT statistics.
- RNASEQ_COLS: support scores and counts from the RNA-seq pipeline.
"""

from typing import List

# ---------------------------------------------------------------------------
# 1. Columns from AGAT statistics
# ---------------------------------------------------------------------------
# These are the exact headers (and order) expected in the AGAT output,
# and used to extract values in the summary table.
AGAT_COLS: List[str] = ["Gene_Models (N)",
                        "Transcript_Models (N)",
                        "CDS_Models (N)",
                        "Exons (N)",
                        "UTR5' (N)",
                        "UTR3' (N)",
                        "Overlapping_Gene_Models (N)",
                        "Single Exon Gene Models (N)",
                        "Single Exon Transcripts (N)",
                        "Total Gene Space (Mb)",
                        "Mean Gene Model Length (bp)",
                        "Mean CDS Model Length (bp)",
                        "Mean Exon Length (bp)",
                        "Mean Intron Length (bp)",
                        "Longest Gene Model Length (bp)",
                        "Longest CDS Model Length (bp)",
                        "Longest Intron Length (bp)",
                        "Shortest Gene Model Length (bp)",
                        "Shortest CDS Model Length (bp)",
                        "Shortest Intron Length (bp)"]

# ---------------------------------------------------------------------------
# 2. Columns from RNA-seq support metrics
# ---------------------------------------------------------------------------
# These headers correspond to the F1-style scores and matching counts
# produced by the StringTie + GFFcompare step.
RNASEQ_COLS: List[str] = ["Transcript level_f1", 
                        "Locus level_f1",
                        "Matching transcripts:",
                        "Matching loci:"]
