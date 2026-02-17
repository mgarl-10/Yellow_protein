

Bioinformatics analysis pipeline for: 

García-Lozano M, Emmerich C, Henzler C, Koch I, Lanz C, Ayas A, Pons I, Buttstedt A, Hipp K, Salem H. Yellow protein co-opted to sustain obligate symbiosis in beetles.


This repository contains analysis pipelines for:

1. **Comparative transcriptomic analysis** of *Chelymorpha alternans*
2. **Phylogenetic analysis** of Yellow protein family genes
3. **Ovary-associated gland transcriptome analysis** for additional Cassidinae species
4. **Co-phylogenetic analysis** between Yellow proteins and _Stammera_ symbionts
5. **Positive selection analysis** for Yellow proteins

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
DESeq2.R
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
- Generates normalized expression values

---

## 3.5 Protein-Coding Gene Prediction (TransDecoder)

```bash
qsub transdecoder.sh
```

This step:
- Identifies candidate coding regions
- Predicts protein sequences from assembled transcripts

# 4. Co-phylogenetic Analysis Between Yellow Proteins and *Stammera* Symbionts

This section describes the reconstruction of host and symbiont phylogenies for co-phylogenetic comparison.

Navigate to the analysis directory:

```bash
cd phylogenies
```

---

## 4.1 Host Phylogeny: Cassidinae Species (7 Species + 1 Outgroup)

Phylogenetic reconstruction was performed using protein sequences from seven Cassidinae species and one outgroup species. 

### Step 1: Multiple Sequence Alignment

Alignment was performed using MUSCLE:

```bash
muscle -in input_proteins.fasta -out aligned_proteins.fasta
```

---

### Step 2: Model Selection

The best-fit amino acid substitution model was determined using ModelTest:

```bash
qsub modeltest.sh
```

---

### Step 3: Phylogenetic Tree Reconstruction

Maximum likelihood phylogeny was reconstructed using RAxML:

```bash
qsub raxml.sh
```

---

## 4.2 *Stammera* Symbiont Phylogeny (Same Cassidinae Species)

Phylogenetic reconstruction of *Stammera* symbionts was conducted using core gene sequences and corresponding outgroups.

---

### Step 1: Extraction of Core Sequences

Core gene sequences were extracted using:

```bash
bash extract_sequences.sh
```
Requires: fasomerecords.py and species.txt

This step:
- Extracts homologous sequences across species
- Generates per-gene FASTA files for alignment

---

### Step 2: Individual Gene Alignments

Each gene was aligned separately using MUSCLE:

```bash
qsub align_muscle_seq.sh
```

---

### Step 3: Concatenation of Alignments

Individual alignments were concatenated into a supermatrix:

```bash
perl catfasta2phyml.pl *.fasta > concatenated_alignment.phy
```

---

### Step 4: Partitioning Scheme and Model Selection

PartitionFinder was used to determine the best partitioning scheme and substitution models:

```bash
qsub partitionfinder.sh
```

---

### Step 5: Phylogenetic Tree Reconstruction

The concatenated alignment was used for maximum likelihood tree reconstruction with RAxML:

```bash
qsub raxml.sh
```

---

## 4.3 Co-phylogenetic Comparison

The resulting trees were visualized side-by-side in Dendroscope3 as a tanglegram using the Neighbor Net Tanglegram algorithm. Cophylogenetic congruence between yellow genes and Stammera symbionts was assessed using eMPRess GUI. 


# 5 Positive selection analysis for Yellow proteins

This section describes the detection of selection signatures in Yellow protein-coding genes.

Navigate to the analysis directory:

```bash
cd positive_selection
```

---

## 5.1 Sequence Retrieval

Protein and corresponding coding sequences (CDS) were obtained from Trinity assemblies using TransDecoder (see Section 3.5).

Outputs used:

- Predicted protein sequences (.pep)
- Predicted coding sequences (.cds)

---

## 5.2 Protein Alignment

Protein sequences were aligned using MUSCLE:

```bash
muscle -in yellow_proteins.fasta -out yellow_proteins_aligned.fasta
```

This step generates amino acid alignments for downstream codon-aware alignment.

---

## 5.3 Codon-Aware Alignment

Codon-aware nucleotide alignments were generated using PAL2NAL:

```bash
qsub pal2nal_alignments.sh
```

This step:
- Combines protein alignment and CDS sequences
- Produces codon-based alignments
- Preserves reading frame for selection analysis

---

## 5.4 Detection of Positive Selection

Positive selection was tested using the **FEL (Fixed Effects Likelihood)** method implemented in Datamonkey (HyPhy framework).

Steps:
1. Upload codon alignment to Datamonkey.
2. Select the FEL model.
3. Use the corresponding maximum likelihood phylogenetic tree.
4. Identify sites under selection based on statistical significance (p-value threshold).

# 6 Alpha-fold modelling 

cd pymol


 
