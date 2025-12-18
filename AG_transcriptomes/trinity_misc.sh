#!/bin/bash
#$ -S /bin/bash
#$ -N trinity_abun
#$ -cwd

TRINITY_UTILS=/ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/funannotate_2/bin

for SPECIES_DIR in */ ; do
    FASTA="${SPECIES_DIR}/Trinity.fasta"
    [[ -f "$FASTA" ]] || continue

    SPECIES="${SPECIES_DIR%/}"

    LEFT="${SPECIES_DIR}/${SPECIES}_trimmed_1.fastq.gz"
    RIGHT="${SPECIES_DIR}/${SPECIES}_trimmed_2.fastq.gz"

    [[ -f "$LEFT" && -f "$RIGHT" ]] || continue

    echo "Processing ${SPECIES}"

    cd "$SPECIES_DIR" || exit 1

    # Trinity stats
    perl "$TRINITY_UTILS/TrinityStats.pl" Trinity.fasta \
        > trinityStats.log

    # Abundance estimation
    perl "$TRINITY_UTILS/align_and_estimate_abundance.pl" \
        --transcripts Trinity.fasta \
        --seqType fq \
        --left "${SPECIES}_trimmed_1.fastq.gz" \
        --right "${SPECIES}_trimmed_2.fastq.gz" \
        --SS_lib_type RF \
        --est_method RSEM \
        --aln_method bowtie \
        --trinity_mode \
        --prep_reference \
        --output_dir "rsem_outdir_${SPECIES}" \
        --coordsort_bam

    cd ..
done
