# GAQET - Genome Annotation Quality Evaluation Tool

GAQET is a Python-based tool designed to evaluate the quality of genome annotations. Using GFF and FASTA files, GAQET generates statistical reports to help identify common errors and artifacts in gene structural annotations.


## Requirements

### Python Dependencies (??????)

- Python == 3.10
- ete3
- PyYAML

### Software dependencies
- AGAT == 1.4.2 (https://github.com/NBISweden/AGAT)
- GFFread == 0.12.7 (https://github.com/gpertea/gffread)
- BUSCO == 5.8.3 (https://github.com/WenchaoLin/BUSCO-Mod)
- LTR_retriever == 3.0.2 (https://github.com/oushujun/LTR_retriever)
- StringTie == 3.0.0 (https://github.com/gpertea/stringtie)

## Installation

Clone this repository:

```bash
git clone https://github.com/victorgcb1987/GAQET.git
```

Create a conda enviroment:

```bash
conda upgrade -c bioconda --all
conda create -c bioconda -n GAQET
conda activate GAQET
conda install -c bioconda agat
conda install -c bioconda gffread
conda install -c bioconda busco
conda install -c bioconda -c conda-forge ltr_retriever
conda install python==3.10
```


## Usage

GAQET uses as a primary input a **File of files**, a tab-separated file with a header row. Each subsequent line defines one sample:
```
| name    | ref_assembly   | ref_annotation   | annotation     | alignments   | lineage        |
|---------|----------------|------------------|----------------|--------------|----------------|
| sample1 | sample1.fa     | sample1.ref.gff3 | sample1.gff3   | sample1.bam  | eudicots_odb10 |
| sample2 | sample2.fa     | sample2.ref.gff3 | sample2.gff3   | sample2.bam  | eudicots_odb10 |

```

| Parameter     | Description                                      |
|---------------|--------------------------------------------------|
| name          | unique sample identifier                         |
| ref_assembly  | path to the reference genome FASTA               |
| ref_annotation| path to the reference annotation GFF/GTF         |
| annotation    | path to the annotation to evaluate (GFF/GTF)     |
| alignments    | path to the RNA-seq alignments (BAM)             |
| lineage       | BUSCO lineage dataset name (e.g. eudicots_odb10) |


With the FOF you can **run GAQET** as follows:

```bash
GAQET.py -i samples.fof -o results/ -t 8
```

-i, --input Path to the FOF (TSV)
-o, --output Output directory (will contain one subfolder per sample)
-t, --threads Number of CPU threads to use (default: 1)