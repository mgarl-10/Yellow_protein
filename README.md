

Bioinformatics analysis pipeline for: 

García-Lozano M, Emmerich C, Henzler C, Koch I, Lanz C, Ayas A, Pons I, Buttstedt A, Hipp K, Salem H. Yellow protein co-opted to sustain obligate symbiosis in beetles.


This repository contains analysis pipelines for:

1. **Comparative transcriptomic analysis** of *Chelymorpha alternans*
2. **Phylogenetic analysis** of Yellow protein family genes
3. **Ovary-associated gland transcriptome analysis** for additional Cassidinae species 

All analyses were performed on a Linux high-performance computing (HPC) cluster using `qsub`.

---

# 1. Comparative Transcriptomics

Navigate to the phylogenetic analysis directory:

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

## 2.1 Sequence Renaming

To standardize FASTA headers and retain species names, sequences were renamed using:

- `fasta_IDS.csv`
- `rename_fasta.py`

```bash
python rename_fasta.py input_sequences.fasta fasta_IDS.csv > renamed_sequences.fasta
```

This step ensures consistent species naming across downstream analyses.

---

## 2.2 Multiple Sequence Alignment

Sequences were aligned using MAFFT:

```bash
mafft --auto renamed_sequences.fasta > aligned_sequences.fasta
```

This step:
- Performs multiple sequence alignment
- Automatically selects appropriate alignment parameters

---

## 2.3 Model Selection

The best-fit amino acid substitution model was determined using ModelTest:

```bash
qsub modeltest.sh
```

This step:
- Evaluates candidate substitution models
- Selects the optimal model based on information criteria

---

## 2.4 Phylogenetic Tree Construction

Phylogenetic reconstruction was performed using RAxML:

```bash
qsub raxml.sh
```

This step:
- Applies the best-fit substitution model
- Performs bootstrap analysis
- Generates the final maximum likelihood tree

---

# 3. Ovary-Associated Gland Transcriptome Analysis for Additional Cassidinae Species

Navigate to the analysis directory:

```bash
cd AG_transcriptomes
```

---

## Input Files

- Paired-end FASTQ files:
  - `*_R1.fastq.gz`
  - `*_R2.fastq.gz`
- `TruSeq3-PE.fa` adapter file

---

## Output Files

This workflow generates:

- Trimmed FASTQ files
- De novo transcriptome assemblies
- BUSCO completeness scores
- Transcript abundance estimates
- Predicted protein-coding sequences

---

## 3.1 Quality Control

### Adapter Removal and Quality Trimming

```bash
qsub trimmomatic.sh
```

This step removes adapter contamination and low-quality bases.

---

## 3.2 De Novo Transcriptome Assembly (Trinity)

```bash
qsub trinity.sh
```

This step:
- Assembles transcriptomes without a reference genome
- Generates assembled transcript FASTA files

---

## 3.3 Assembly Quality Assessment (BUSCO)

```bash
qsub busco.sh
```

This step:
- Evaluates assembly completeness
- Reports single-copy ortholog recovery statistics

---

## 3.4 Transcript Abundance Estimation (RSEM)

```bash
qsub trinity_misc.sh
```

This step:
- Estimates transcript abundance
- Generates normalized expression values (e.g., TPM)

---

## 3.5 Protein-Coding Gene Prediction (TransDecoder)

```bash
qsub transdecoder.sh
```

This step:
- Identifies candidate coding regions
- Predicts protein sequences from assembled transcripts


