"""
LTR_retriever.py
================
Wrappers for the GeneTools pipeline used to compute LAI (LTR Assembly Index):

* **suffixerator**   - builds a suffix-array index of the genome.
* **ltrharvest**     - detects LTR retrotransposons (harvest step).
* **LTR_FINDER**     - complementary detector (finder step).
* **cat**            - concatenates harvest + finder output.
* **LTR_retriever**  - filters and refines candidates.
* **LAI**            - computes the final LTR Assembly Index score.

Each runner returns a dictionary with the executed command, an informational message,
the main output path, and a ``returncode`` (99 means “already done”).
"""


import os
import subprocess
import shutil
from pathlib import Path
from typing import Any, Dict


# ---------------------------------------------------------------------------
# 0. Prepare working directory for LAI
# ---------------------------------------------------------------------------
def create_outdir(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Ensure the LAI working directory exists and symlink the genome."""

    # Output directory (to save LTR_retriever input and output files)
    outdir = arguments["output"] / "LAICompleteness"
    outfile = outdir / Path(arguments["ref_assembly"]).name

    if not outdir.exists():
        outdir.mkdir(parents=True)
    if not outfile.exists():
        cmd = f"ln -s {str(arguments['ref_assembly'])} {str(outfile)}"
        run_ = subprocess.run(cmd, shell=True)
    msg = "The output directory for LAICompleteness has been created"

    return {"msg": msg, 
            "out_fpath": outdir}

# ---------------------------------------------------------------------------
# 1. Build suffix-array index
# ---------------------------------------------------------------------------
def run_suffixerator(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Run **gt suffixerator** or skip if the index already exists."""

    # Output index path
    index = arguments["LAI_dir"] / Path(arguments["ref_assembly"]).name

    # Suffixerator command
    cmd = "gt suffixerator -db {} -indexname {} -tis -suf -lcp -des -ssp -sds -dna".format(arguments["ref_assembly"], index)

    #Check if suffixerator is already done
    md5 = arguments["LAI_dir"] / "{}.md5".format(Path(arguments["ref_assembly"]).name)
    if md5.exists():
        return {"command": cmd,
                "msg": "suffixerator already done",
                "out_fpath": index,
                "returncode": 99}
    else:
        #Run suffixerator
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)

        if run_.returncode == 0:
            msg = "suffixerator ran successfully"
        else:
            msg = "suffixerator Failed: \n {}".format(run_.stderr)

        return {"command": cmd,
                "msg": msg,
                "out_fpath": index, 
                "returncode": run_.returncode}

# ---------------------------------------------------------------------------
# 2. ltrharvest
# ---------------------------------------------------------------------------
def run_harvest(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Run *gt ltrharvest* or skip if done."""

    # Take suffixerator index to use it as input
    index = arguments["LAI_dir"] / Path(arguments["ref_assembly"]).name

    # Create output: output path + file .harvest.scn /w ref_assembly file name
    out = arguments["LAI_dir"] / "{}.harvest.scn".format(Path(arguments["ref_assembly"]).name)

    # HARVEST command
    cmd = "gt ltrharvest -index {} -minlenltr 100 -maxlenltr 7000 -mintsd 4 -maxtsd 6 -motif TGCA -motifmis 1 -similar 85 -vic 10 -seed 20 -seqids yes > {}".format(index, out)

    # Check if HARVEST is already done
    if out.exists():
        return {"command": cmd, 
                "msg": "harvest already done",
                "out_fpath": out,
                "returncode": 99}
 
    else:
        #Run harvest
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
 
        if run_.returncode == 0:
            msg = "HARVEST ran successfully"
        else:
            msg = "HARVEST Failed: \n {}".format(run_.stderr)

        return {"command": cmd, 
                "msg": msg,
                "out_fpath": out, 
                "returncode": run_.returncode}

# ---------------------------------------------------------------------------
# 3. LTR_FINDER_parallel
# ---------------------------------------------------------------------------
def run_finder(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Run *LTR_FINDER_parallel* or skip if done."""
    cwd = Path(os.getcwd())

    # FINDER command
    cmd = "LTR_FINDER_parallel -seq {} -threads {} -harvest_out -size 1000000 -time 300".format(arguments["ref_assembly"],
                                                                                                arguments["threads"])

    # Check if FINDER is already done
    out_file = arguments["LAI_dir"] / "{}.finder.combine.scn".format(Path(arguments["ref_assembly"]).name)
    if out_file.exists():
        return {"command": cmd,
                "msg": "harvest already done",
                "out_fpath": out_file,
                "returncode": 99}

    else:
        # Change the working directory to the "output" path
        os.chdir(arguments["LAI_dir"])
        # Run FINDER
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)

        if run_.returncode == 0:
            msg = "FINDER ran successfully"
        else:
            msg = " FINDER Failed: \n {}".format(run_.stderr)

        #Restore the original working directory
        os.chdir(cwd)

        return {"command": cmd, 
                "msg": msg,
                "out_fpath": out_file, 
                "returncode": run_.returncode}

# ---------------------------------------------------------------------------
# 4. Concatenate harvest + finder output
# ---------------------------------------------------------------------------
def concatenate_outputs(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Concatenate harvest & finder results or skip if done."""

    # "cat" command
    outpath = arguments["LAI_dir"] / Path(arguments["ref_assembly"]).name
    cmd = "cat {}.harvest.scn {}.finder.combine.scn > {}.rawLTR.scn".format(outpath, 
                                                                            outpath, 
                                                                            outpath)

    # Check if "cat" is already done
    out_file = arguments["LAI_dir"] / "{}.rawLTR.scn".format(Path(arguments["ref_assembly"]).name)
    if out_file.exists():
        return {"command": cmd, 
                "msg": "Concatenation of the output files from Harvest and Finder is already done.",
                "out_fpath": out_file,
                "returncode": 99}

    else:
        # Run "cat"
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)

        if run_.returncode == 0:
            msg = "Concatenation of the output files from Harvest and Finder successfully completed."
        else:
            msg = "Failed: \n {}".format(run_.stderr)

        return {"command": cmd,
                "msg": msg,
                "out_fpath": out_file, 
                "returncode": run_.returncode}

# ---------------------------------------------------------------------------
# 5. LTR_retriever refinement
# ---------------------------------------------------------------------------
def run_LTR_retriever(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Run *LTR_retriever* or skip if done."""
    cwd = Path(os.getcwd())

    # LTR_retriever command
    cmd = "LTR_retriever -genome {} -inharvest {}.rawLTR.scn -threads {}".format(Path(arguments["ref_assembly"]).name,
                                                                                Path(arguments["ref_assembly"]).name,
                                                                                arguments["threads"])

    # Check if LTR_retriever is already done
    outfile = arguments["LAI_dir"] / "{}.mod.pass.list".format(Path(arguments["ref_assembly"]).name)
    if outfile.exists():
        return {"command": cmd, 
                "msg": "LTR_retriever already done",
                "out_fpath": outfile,
                "returncode": 99}

    else:
        # Change the working directory to the "output" path
        os.chdir(arguments["LAI_dir"])
        # Run LTR_retriever
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)

        if run_.returncode == 0:
            msg = "LTR_retriever ran successfully"
        else:
            msg = "LTR_retriever Failed: \n {}".format(run_.stderr)

        #Restore the original working directory
        os.chdir(cwd)

        return {"command": cmd, 
                "msg": msg,
                "out_fpath": outfile, 
                "returncode": run_.returncode}

# ---------------------------------------------------------------------------
# 6. Compute LAI
# ---------------------------------------------------------------------------
def run_LAI(arguments: Dict[str, Any]) -> Dict[str, Any]:
    """Run *LAI* tool or skip if done."""
    cwd = Path(os.getcwd())

    # LAI command
    cmd = "LAI -genome {} -intact {}.mod.pass.list -all {}.mod.out".format(Path(arguments["ref_assembly"]).name,
                                                                            Path(arguments["ref_assembly"]).name,
                                                                            Path(arguments["ref_assembly"]).name)


    outfile = arguments["LAI_dir"] / "{}.mod.out.LAI".format(Path(arguments["ref_assembly"]).name)
    if outfile.exists():
        return {"command": cmd, 
                "msg": "LAI already done",
                "out_fpath": outfile,
                "returncode": 99}

    else:
        # Change the working directory to the "output" path
        os.chdir(arguments["LAI_dir"])
        # Run
        run_ = subprocess.run(cmd, shell=True, stderr=subprocess.PIPE)
        if run_.returncode == 0:
            msg = "LAI ran successfully"
        else:
            msg = "LAI Failed: \n {}".format(run_.stderr)

        #Restore the original working directory
        os.chdir(cwd)

        return {"command": cmd, 
                "msg": msg,
                "out_fpath": outfile, 
                "returncode": run_.returncode}

# ---------------------------------------------------------------------------
# 7. Retrieve LAI value
# ---------------------------------------------------------------------------
def get_LAI(lai_run: Dict[str, Any]) -> str:
    """Return ``"LAI:<value>"`` from the final *.LAI* file."""

    with open(lai_run["out_fpath"], encoding="utf-8") as fh:
        for line in fh:
            if "whole_genome" in line:
                line = line.strip()
                line = line.split()
                LAI_value = line[6]
                return f"LAI:{LAI_value}"