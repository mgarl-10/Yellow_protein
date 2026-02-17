

Bioinformatics analysis pipeline for: 

García-Lozano M, Emmerich C, Henzler C, Koch I, Lanz C, Ayas A, Pons I, Buttstedt A, Hipp K, Salem H. Yellow protein co-opted to sustain obligate symbiosis in beetles.


This repository contains analysis pipelines for:

1. **Comparative transcriptomic analysis** of *Chelymorpha alternans*
2. **Phylogenetic analysis** of Yellow protein family genes
3. 

All analyses were performed on a Linux high-performance computing (HPC) cluster using `qsub`.


# 1. Comparative Transcriptomics

Navigate to the transcriptomics workflow:

```bash
cd comparative_transcriptomics
```

## Input Files

The following files are required:

- Paired-end FASTQ files  
  - `*_R1.fastq.gz`
  - `*_R2.fastq.gz`
- Reference genome (FASTA format)
- Gene annotation file (gff3 format)
- `TruSeq3-PE.fa` adapter file

## Output Files

The workflow generates:

- Trimmed FASTQ files
- Sorted and indexed BAM files
- Gene-level count matrix
- Differential expression results table
- Volcano plot

## 1.1 Quality Control

### Adapter removal and quality trimming

```bash
qsub trimmomatic.sh
```

This step:
- Removes Illumina adapter sequences
- Performs quality trimming
- Generates cleaned paired-end FASTQ files

---

## 1.2 Mapping to Reference Genome

### Step 1: Index the reference genome

```bash
hisat2-build genome.fa genome_index
```

### Step 2: Align trimmed reads to *Chelymorpha alternans* genome

```bash
qsub hisat2.sh
```

This step:
- Aligns paired-end reads to the reference genome
- Generates SAM alignment files

### Step 3: Convert, sort, and index BAM files

```bash
qsub sam_formatting.sh
```

This step:
- Converts SAM to BAM
- Sorts BAM files
- Indexes BAM files for downstream analysis

---

## 1.3 Quantification of Mapped Reads

```bash
qsub htseq_exon.sh
```

This step:
- Uses HTSeq to count reads mapped to exons
- Generates a gene-level count matrix for downstream analysis

---

## 1.4 Differential Gene Expression Analysis

```bash
Rscript DESeq2.R
```

This script:
- Performs differential expression analysis using DESeq2
- Outputs normalized counts
- Produces statistical results (log2FC, p-values, adjusted p-values)
- Generates a volcano plot

---

# 2. Phylogenetic Analysis of Yellow Proteins

Navigate to the phylogenetic analysis directory:

```bash
cd phylogenies/Yellow_protein_tree
```

---

## 2.1 Multiple Sequence Alignment

Sequences were aligned using MAFFT:

```bash
mafft --auto input_sequences.fasta > aligned_sequences.fasta
```

---

## 2.2 Phylogenetic Tree Construction

(Include your tree-building software here if applicable, e.g., IQ-TREE, RAxML, or FastTree.)

Example:

```bash
iqtree -s aligned_sequences.fasta -m TEST -bb 1000
```

This step:
- Selects the best-fit substitution model
- Performs bootstrap analysis
- Generates a phylogenetic tree file

---

# Reproducibility Notes

- All scripts were executed in a Linux HPC environment.
- File paths may need to be adjusted depending on system configuration.
- Ensure software versions are compatible.

---

# Citation

If you use this pipeline in your research, please cite:

Author et al., Year, Journal (if applicable)

---

# Contact

For questions or issues, please open a GitHub issue or contact:

Your Name  
your.email@institution.edu






