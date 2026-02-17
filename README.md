# Bioinformatics workflows supporting:

García-Lozano M, Emmerich C, Henzler C, Koch I, Lanz C, Ayas A, Pons I, Buttstedt A, Hipp K, Salem H.  
**Yellow protein co-opted to sustain obligate symbiosis in beetles.**


This repository contains analysis pipelines for:

1. **Comparative transcriptomic analysis** of *Chelymorpha alternans*
2. **Phylogenetic analysis** of Yellow protein family genes
3. **Ovary-associated gland transcriptome analysis** for additional Cassidinae species
4. **Co-phylogenetic analysis** between Yellow proteins and _Stammera_ symbionts
5. **Positive selection analysis** for Yellow proteins
6. **AlphaFold modeling and comparison** of Yellow proteins

---

## Software Requirements

Most of  analyses were performed on a Linux HPC cluster using `qsub`.

Core tools:
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


# 1. Comparative Transcriptomics

Navigate to the transcriptomics analysis directory:

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
cd phylogenies/Yellow_protein_Cassidinae
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

## 4.2 *Stammera* Symbiont Phylogeny 

Phylogenetic reconstruction of *Stammera* symbionts was conducted using core gene sequences and corresponding outgroups.

Navigate to the analysis directory:

```bash
cd phylogenies/Stammera_phylogeny
```

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

PartitionFinder was used to infer the optimal partitioning scheme and substitution models.

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

# 6. AlphaFold Structural Modeling and Structural Comparison

This section describes structural prediction and comparative analysis of Yellow proteins across seven Cassidinae species.

---

## 6.1 Structure Prediction

Predicted protein structures were generated using the AlphaFold Server (AlphaFold 3).

```bash
cd pymol
```

- The seven Cassidinae Yellow protein sequences were submitted individually.
- CIF files corresponding to **model 0** were downloaded for downstream structural analyses.

---

## 6.2 Structural Alignment in PyMOL

The predicted structures were imported into PyMOL for structural comparison.

The *C. alternans* Yellow protein was selected as the reference structure.

### Pairwise Structural Alignment

```python
align model_5, model_1
align model_5, model_2
align model_5, model_3
align model_5, model_4
align model_5, model_6
align model_5, model_7
```

---

## 6.3 Per-Residue RMSD Calculation

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

## 6.4 Model Confidence (pLDDT) Analysis

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

## 6.5 Mapping Functional and Evolutionary Sites

### Hydrophobic Residues

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

### Positively Selected Sites

Codon sites identified as evolving under positive selection (see Section 5) were mapped onto the predicted structure.

```python
select pos_sites, resi 156+189+226+271+276+281+415+700
show spheres, pos_sites
color red, pos_sites
```

---

## 6.6 Model Confidence Visualization (pLDDT)

AlphaFold confidence scores (pLDDT), stored in the B-factor field of the predicted structures, were visualized in PyMOL using the standard AlphaFold color scheme (blue = high confidence).

```python
hide everything
show cartoon

# Apply AlphaFold-style confidence coloring
spectrum b, red_orange_yellow_cyan_blue, minimum=0, maximum=100
```

Color interpretation:
- **Blue**: very high confidence (pLDDT > 90)
- **Cyan**: confident (70–90)
- **Yellow**: low confidence (50–70)
- **Orange/Red**: very low confidence (< 50)

