# Bioinformatics workflows and statistical analyses

This repository contains the computational and statistical analyses supporting:

García-Lozano M, Emmerich C, Henzler C, Koch I, Lanz C, Ayas A, Pons I, Buttstedt A, Hipp K, Salem H.  
**Yellow protein co-opted to sustain obligate symbiosis in beetles.**

---

## Overview


The analyses include:

1. **Assembly, decontamination, and annotation** for _Chelymorpha alternans_ genome
2. **Comparative transcriptomic analysis** of symbiotic organs in *Chelymorpha alternans*
3. **Phylogenetic reconstruction** of Yellow protein family genes
4. **Ovary-associated gland transcriptome analysis** for additional Cassidinae species
5. **Co-phylogenetic analysis** between Yellow proteins and *Stammera*
6. **Positive selection analysis** of Yellow proteins
7. **AlphaFold structural modeling and structural comparison**
8. **Symbiont transcriptome response** under low humidity conditions
9. **Statistical analyses and figure generation**

---

## Software Requirements

Most of the analyses were performed on a Linux HPC cluster using `qsub`.

Core tools:
- Hifiasm
- Redundans
- Seqkit
- BWA  
- BLAST+  
- BlobTools  
- RepeatModeler  
- RepeatMasker  
- GeneMark  
- Funannotate
- Trimmomatic
- HISAT2
- SAMtools
- HTSeq
- DESeq2
- Trinity
- BUSCO
- RSEM
- TransDecoder
- MAFFT
- MUSCLE
- RAxML
- PartitionFinder
- PAL2NAL
- Datamonkey (HyPhy FEL)
- AlphaFold Server (AF3)
- PyMOL
- Dendroscope3
- eMPRess
- Bowtie
- featureCounts
- R

## Important:

Resource specifications such as:

```bash
qsub -l h_vmem=XXG -pe parallel XX
```

are not included in the submitted scripts, as memory limits (h_vmem) and parallel environments (-pe) vary substantially across HPC systems.

Users should therefore modify job submission headers in each *.sh file according to their local cluster configuration and resource availability.


# 1. Genome assembly for *Chelymorpha alternans*


```bash
cd genome_assembly
```

---

## 1.1 Assembly

### 1.1.1 Primary assembly (HiFi reads)

Genome assembly was performed using **Hifiasm**, followed by redundancy reduction with **Redundans**.

**Input:**
- PacBio HiFi reads (`*.fasta`)

```bash
qsub hifiasm.sh
qsub redun.sh
qsub busco.sh
```

**Output:**
- Primary genome assembly (FASTA)
- Purged/redundancy-reduced assembly (FASTA)
- Assembly statistics
- Genome completeness statistics

---

## 1.2 Decontamination

```bash
cd decontamination
```

### 1.2.1 Read mapping

Reads were mapped back to the assembly using **BWA-MEM** to assess coverage and support contamination screening.

**Input:**
- Redundancy-reduced genome assembly (FASTA)
- PacBio HiFi reads (FASTA/FASTQ)

```bash
qsub mapping.sh
```

**Output:**
- Sorted and indexed BAM file
- Mapping statistics

---

### 1.2.2 Contamination screening

Taxonomic screening and decontamination were performed using **BLASTn** and **BlobTools**.

```bash
qsub blast.sh
qsub blobtools.sh
```

**Input:**
- Genome assembly (FASTA)
- BAM mapping file
- NCBI nt database

**Output:**
- BLAST tabular output
- BlobTools database (`blobDB.json`)
- Taxon-annotated coverage plots
- Filtered/decontaminated genome assembly

---

## 1.3 Annotation

```bash
cd ../annotation
```

Prior to repeat identification and gene prediction, the genome assembly was cleaned and standardized using Funannotate utilities.

---

### 1.3.1 Genome cleaning and sorting

Small contigs (< 1000 bp) were removed and contigs were renamed for consistency.

```bash
funannotate clean \
-i C_alternans_assembly_filtered.fasta \
--minlen 1000 \
-o C_alternans_genome_cleaned.fa \
--cpus 32

funannotate sort \
-i C_alternans_genome_cleaned.fa \
-b contig \
-o C_alternans_genome_cleaned_sorted.fa
```

**Input:**
- Decontaminated genome assembly (`C_alternans_assembly_filtered.fasta`)

**Output:**
- Length-filtered genome assembly (`C_alternans_genome_cleaned.fa`)
- Sorted and standardized genome assembly (`C_alternans_genome_cleaned_sorted.fa`)

---


### 1.3.2 Repeat identification and masking

Repeat annotation was performed using **RepeatModeler v2.0.3** to build a species-specific repeat library, followed by genome masking using **Funannotate (RepeatMasker backend)**.

---

#### 1.3.2.1 Build RepeatModeler database

```bash
repeatmodeler.sh
```

**Input:**
- Cleaned and sorted genome assembly (`C_alternans_genome_cleaned_sorted.fa`)

**Output:**
- RepeatModeler database files

---

#### 1.3.2.2 De novo repeat identification

```bash
repeatmodeler.sh
```

**Input:**
- RepeatModeler database

**Output:**
- Species-specific repeat library  
  (`consensi.fa.classified`)

---

#### 1.3.2.3 Genome masking

The custom repeat library was used to mask the genome via Funannotate (RepeatMasker backend).

```bash
funannotate mask \
-i C_alternans_genome_cleaned_sorted.fa \
-m repeatmasker \
-l ./RM_XXXX/consensi.fa.classified \
--cpus 64 \
-o C_alternans_assembly.fa
```

**Input:**
- Cleaned and sorted genome assembly
- Custom repeat library (`consensi.fa.classified`)

**Output:**
- Soft-masked genome assembly (`C_alternans_assembly.fa`)
- Repeat annotation files

---


## 1.4 Gene prediction

### 1.4.1 GeneMark

Ab initio gene prediction was performed using GeneMark-ES.

```bash
qsub genemark.sh
```

**Input:**
- Repeat-masked genome assembly

**Output:**
- Predicted gene models (GTF/GFF format)
- Predicted coding sequences (CDS)

---

### 1.4.2 Funannotate

Structural and functional annotation were conducted using Funannotate.

```bash
qsub fun.sh
```

**Input:**
- Repeat-masked genome assembly
- GeneMark predictions
- Optional: RNAseq evidence

**Output:**
- Final annotated gene models (GFF3)
- Predicted proteins (FAA)
- Coding sequences (FNA)
- Functional annotations (InterPro, Pfam, GO terms)

# 2. Comparative transcriptomics

Comparative transcriptomics between foregut symbiotic organs and ovary-associated glands in _Chelymorpha alternans_

```bash
cd comparative_transcriptomics
```

## Input files

The following files are required:

- Paired-end FASTQ files  
  - `*_R1.fastq.gz`
  - `*_R2.fastq.gz`
- Reference genome (FASTA format)
- Gene annotation file (gff3 format)
- `TruSeq3-PE.fa` adapter file

## Output files

The workflow generates:

- Trimmed FASTQ files
- Sorted and indexed BAM files
- Gene-level count matrix
- Differential expression results table
- Volcano plot

## 2.1 Quality control

### Adapter removal and quality trimming

```bash
qsub trim.sh
```

This step:
- Removes Illumina adapter sequences
- Performs quality trimming
- Generates cleaned paired-end FASTQ files

---

## 2.2 Mapping to reference genome

### Step 1: Index the reference genome

```bash
hisat2-build-l C_alternans_genome.fa C_alternans
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

## 2.3 Quantification of mapped Reads

```bash
qsub htseq_exon.sh
```

This step:
- Uses HTSeq to count reads mapped to exons
- Generates a gene-level count matrix for downstream analysis

---

## 2.4 Differential gene expression analysis

```bash
Rscript DESeq2.R
```

This script:
- Performs differential expression analysis using DESeq2
- Outputs normalized counts
- Produces statistical results (log2FC, p-values, adjusted p-values)
- Generates a volcano plot

---

# 3. Phylogenetic analysis of Yellow proteins


```bash
cd phylogenies/Yellow_protein_tree
```

---

## 3.1 Sequence renaming

To standardize FASTA headers and retain species names, sequences were renamed using:

- `fasta_IDS.csv`
- `rename_fasta.py`

```bash
python rename_fasta.py
```

This step ensures consistent species naming across downstream analyses.

---

## 3.2 Multiple sequence alignment

Sequences listed in the fasta IDs.csv file and Drosophila sequences obtained from NCBI were aligned using MAFFT:

```bash
mafft --auto yellow_filtered_coleopt_hymenopt_chely_bact_Drosophila_renamed.fasta > yellow_filtered_coleopt_hymenopt_chely_bact_Drosophila_renamed_alignment.fasta
```
*The alignment was converted to PHYLIP format in Geneious.


---

## 3.3 Model selection

The best-fit amino acid substitution model was determined using ModelTest:

```bash
qsub modeltest.sh
```

This step:
- Evaluates candidate substitution models
- Selects the optimal model based on information criteria

---

## 3.4 Phylogenetic tree construction

Phylogenetic reconstruction was performed using RAxML:

```bash
qsub raxml.sh
```

This step:
- Applies the best-fit substitution model
- Performs bootstrap analysis
- Generates the final maximum likelihood tree

---

# 4. Ovary-associated gland transcriptome analysis for additional Cassidinae species


```bash
cd AG_transcriptomes
```

---

## Input files

- Paired-end FASTQ files:
  - `*_R1.fastq.gz`
  - `*_R2.fastq.gz`
- `TruSeq3-PE.fa` adapter file

---

## Output files

This workflow generates:

- Trimmed FASTQ files
- De novo transcriptome assemblies
- BUSCO completeness scores
- Transcript abundance estimates
- Predicted protein-coding sequences

---

## 4.1 Quality control

### Adapter removal and quality trimming

```bash
qsub trim.sh
```

This step removes adapter contamination and low-quality bases.

---

## 4.2 De novo transcriptome assembly (Trinity)

```bash
qsub trinity.sh
```

This step:
- Assembles transcriptomes without a reference genome
- Generates assembled transcript FASTA files

---

## 4.3 Assembly quality assessment (BUSCO)

```bash
qsub busco.sh
```

This step:
- Evaluates assembly completeness

---

## 4.4 Transcript abundance estimation (RSEM)

```bash
qsub trinity_misc.sh
```

This step:
- Estimates transcript abundance
- Generates normalized expression values

---

## 4.5 Protein-coding gene prediction (TransDecoder)

```bash
qsub transdecoder.sh
```

This step:
- Identifies candidate coding regions
- Predicts protein sequences from assembled transcripts

# 5. Co-phylogenetic analysis between Yellow proteins and *Stammera* 

This section describes the reconstruction of host and symbiont phylogenies for co-phylogenetic comparison.


```bash
cd phylogenies/Yellow_protein_Cassidinae
```

---

## 5.1 Host phylogeny: Cassidinae species (7 species + 1 outgroup)

Phylogenetic reconstruction was performed using protein sequences from seven Cassidinae species and one outgroup species. 

### Step 1: Multiple sequence alignment

Alignment was performed using MUSCLE:

```bash
muscle -align yellow_protein_outgroup.fasta -output yellow_protein_alignment_outgroup.fasta
```
*The alignment was converted to PHYLIP format in Geneious.

---

### Step 2: Model selection

The best-fit amino acid substitution model was determined using ModelTest:

```bash
qsub modeltest.sh
```

---

### Step 3: Phylogenetic tree reconstruction

Maximum likelihood phylogeny was reconstructed using RAxML:

```bash
qsub raxml.sh
```

---

## 5.2 *Stammera* phylogeny 

Phylogenetic reconstruction of *Stammera* was conducted using core gene sequences and corresponding outgroups.


```bash
cd phylogenies/Stammera_phylogeny
```

---

### Step 1: Extraction of core sequences

Core gene sequences were extracted using:

```bash
bash extract_sequences.sh
```
Requires faSomeRecords.py and species.txt files 

This step:
- Extracts homologous sequences across species
- Generates per-gene FASTA files for alignment

---

### Step 2: Individual gene alignments

Each gene was aligned separately using MUSCLE:

```bash
qsub align_muscle_seq.sh
```
*The alignment was converted to PHYLIP format in Geneious.

---

### Step 3: Concatenation of alignments

Individual alignments were concatenated into a supermatrix:

```bash
perl catfasta2phyml.pl *.fasta > concatenated_aln_yellow_species.phy 2> partitions.txt
```

---

### Step 4: Partitioning scheme and model selection

PartitionFinder was used to infer the optimal partitioning scheme and substitution models. 

```bash
qsub partitionfinder.sh
```

---

### Step 5: Phylogenetic tree reconstruction

The concatenated alignment was used for maximum likelihood tree reconstruction with RAxML:

```bash
qsub raxml.sh
```

---

## 5.3 Co-phylogenetic comparison

The resulting trees were visualized side-by-side in Dendroscope3 as a tanglegram using the Neighbor Net Tanglegram algorithm. Cophylogenetic congruence between yellow genes and _Stammera_ was assessed using eMPRess GUI. 


# 6 Positive selection analysis of Yellow proteins

This section describes the detection of selection signatures in Yellow protein-coding genes.


```bash
cd positive_selection
```

---

## 6.1 Sequence retrieval

Protein and corresponding coding sequences (CDS) were obtained from Trinity assemblies using TransDecoder (see Section 3.5).

Outputs used:

- Predicted protein sequences (.pep)
- Predicted coding sequences (.cds)

---

## 6.2 Protein alignment

Protein sequences were aligned using MUSCLE:

```bash
muscle -align yellow_proteins_7species.fasta -output protein_aln_7species.fasta
```

This step generates amino acid alignments for downstream codon-aware alignment.

---

## 6.3 Codon-aware alignment

Codon-aware nucleotide alignments were generated using PAL2NAL:

```bash
qsub pal2nal_alignments.sh
```

This step:
- Combines protein alignment and CDS sequences
- Produces codon-based alignments
- Preserves reading frame for selection analysis

---

## 6.4 Detection of positive selection

Positive selection was tested using the **FEL (Fixed Effects Likelihood)** method implemented in Datamonkey (HyPhy framework).

Steps:
1. Upload codon alignment to Datamonkey.
2. Select the FEL model.
3. Use the corresponding maximum likelihood phylogenetic tree.
4. Identify sites under selection based on statistical significance (p-value threshold).

# 7. AlphaFold structural modeling and structural comparison

This section describes structural prediction and comparative analysis of Yellow proteins across seven Cassidinae species.

---

## 7.1 Structure prediction

Predicted protein structures were generated using the AlphaFold Server (AlphaFold 3).

```bash
cd pymol
```

- The seven Cassidinae Yellow protein sequences were submitted individually.
- CIF files corresponding to **model 0** were downloaded for downstream structural analyses.

---

## 7.2 Structural alignment in PyMOL

The predicted structures were imported into PyMOL for structural comparison.

The *C. alternans* Yellow protein was selected as the reference structure.

### Pairwise structural alignment

```python
align model_5, model_1
align model_5, model_2
align model_5, model_3
align model_5, model_4
align model_5, model_6
align model_5, model_7
```

---

## 7.3 Per-residue RMSD calculation

A custom PyMOL script (`pymol_aln.py`) was used to calculate per-residue RMSD values by measuring distances between Cα atoms across aligned models.

```python
run pymol_aln.py
multi_rmsd_residuewise("model_5", ["model_1", "model_2", "model_3", "model_4", "model_6", "model_7"])
```

RMSD values were:

- Stored in the **B-factor field** of the reference model
- Visualized using a red–blue color gradient:
  - **Red** → low structural variability
  - **Blue** → high structural variability

Legend generation:

```python
run spectrumbar_rmsd.py
spectrumbar_rmsd name=legend_rmsd, start="0,0,0", length=5, radius=0.5, steps=5, rmsd_min=10.0, rmsd_max=95.0
```

Structure representation:

```python
cartoon oval
```

---

## 7.4 Model confidence (pLDDT) analysis

Model confidence scores (pLDDT) were extracted from AlphaFold predictions and summarized by exon region.

Exon boundaries were defined based on the *C. alternans* Yellow gene model.

For the remaining species:
- Exon 1 and exons 2–3 were inferred by aligning transcript-derived sequences to the *C. alternans* reference
- This was necessary due to the absence of annotated exon structures

Extraction was performed using:

```bash
python extract_plddt_with_regions.py
```

---

## 7.5 Mapping functional and evolutionary sites

### Hydrophobic residues

Hydrophobic residues were identified and highlighted in PyMOL based on amino acid identity:

- Ala
- Val
- Leu
- Ile
- Met
- Phe
- Trp
- Pro

```python
load C_alternans.pdb
hide everything
show cartoon
color gray

select hydro_res, resn ALA+VAL+LEU+ILE+MET+PHE+TRP+PRO
show sticks, hydro_res
color yellow, hydro_res
```

---

### Positively selected sites

Codon sites identified as evolving under positive selection (see Section 5) were mapped onto the predicted structure.

```python
select pos_sites, resi 156+189+226+271+276+281+415+700
show spheres, pos_sites
color red, pos_sites
```

---

## 7.6 Model confidence visualization (pLDDT)

AlphaFold confidence scores (pLDDT), stored in the B-factor field of the predicted structures, were visualized in PyMOL using the standard AlphaFold color scheme (blue = high confidence).

```python
hide everything
show cartoon

spectrum b, red_orange_yellow_cyan_blue, minimum=0, maximum=100
```

Color interpretation:
- **Blue**: very high confidence (pLDDT > 90)
- **Cyan**: confident (70–90)
- **Yellow**: low confidence (50–70)
- **Orange/Red**: very low confidence (< 50)


# 8. Symbiont transcriptome response under low humidity conditions

This workflow processes *Stammera* RNA-seq data to quantify gene expression and assess differential transcriptional responses under low humidity conditions.


```bash
cd symbiont_transcriptomics
```

---

## 8.1 Input files

The following files are required:

- Single-end FASTQ files  
  - `*_R1.fastq.gz`
- Symbiont reference genome (FASTA format)
- Gene annotation file (GFF format)
- `TruSeq3-PE.fa` adapter file (for trimming)

---

## 8.2 Output files

This workflow generates:

- Trimmed FASTQ files
- Alignment files (SAM format)
- Gene-level count matrix
- Differential expression results tables
- Volcano plots (per treatment comparison)
- Heatmap for the Yellow condition

---

## 8.3 Quality control

### Adapter removal and quality trimming

```bash
qsub trim.sh
```


---

## 8.4 Mapping to reference genome

### Step 1: Index the reference genome

```bash
bowtie2-build Stammera_Chelymorpha_alternans.fasta Stammera
```

### Step 2: Align trimmed reads to the *Stammera* genome

```bash
qsub map.sh
```

This step:
- Aligns single-end reads to the symbiont reference genome using Bowtie2
- Generates SAM alignment files

---

## 8.5 Quantification of mapped reads

```bash
qsub featurecounts.sh
```

This step:
- Assigns mapped reads to annotated genes
- Generates a gene-level count matrix for downstream differential expression analysis

---

## 8.6 Differential gene expression analysis

```bash
Rscript DESeq_script_Stammera.R
```

This script:
- Performs differential expression analysis using DESeq2
- Tests treatment-specific responses under low humidity conditions
- Produces statistical results (log2FC, p-values, adjusted p-values)
- Generates volcano plots for each treatment comparison
- Generates a heatmap with annotations for the Yellow condition

# 9. Statistical analyses and figure generation

All statistical analyses and figure generation were performed in **R**.

```bash
cd stats
```

---

## Overview

This directory contains R scripts used to generate statistical analyses and figures presented in:

- Figure 2
- Figure 3
- Figure 4
- Figure 5
- Figure 6
- Figure 7
- Figure S8

Each script reproduces the statistical tests and visualizations corresponding to its respective figure.

---

## Input data

Raw data files required to execute the scripts are located in:

```
stats/raw_files/
```

---

## Scripts

- `Stats_Fig_2.R`
- `Stats_Fig_3.R`
- `Stats_Fig_4.R`
- `Stats_Fig_5.R`
- `Stats_Fig_6.R`
- `Stats_Fig_7.R`
- `Stats_Fig_S8.R`

---

## Output

Execution of these scripts generates:

- Publication-ready figures (PDF/svg)
- Statistical summaries

---


## Citation

If you use this repository, please cite:

García-Lozano et al. (2026). Yellow protein co-opted to sustain obligate symbiosis in beetles.

Zenodo DOI: 10.5281/zenodo.XXXXXXX



