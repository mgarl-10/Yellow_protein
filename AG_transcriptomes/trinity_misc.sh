#!/bin/bash
#$ -S /bin/bash
#$ -N trinity_abun
#$ -cwd

trinity_misc=/Trinity/bin/

for species_dir in */ ; do
    fasta="${species_dir}/Trinity.fasta"
    [[ -f "$fasta" ]] || continue

    species="${species_dir%/}"

    left="${species_dir}/${species}_trimmed_1.fastq.gz"
    right="${species_dir}/${species}_trimmed_2.fastq.gz"

    [[ -f "$left" && -f "$right" ]] || continue

    echo "Processing ${species}"

    cd "$species_dir" || exit 1

    # Trinity stats
    perl "$trinity_misc/TrinityStats.pl" Trinity.fasta \
        > trinityStats.log

    # Abundance estimation
    perl "$trinity_misc/align_and_estimate_abundance.pl" \
        --transcripts Trinity.fasta \
        --seqType fq \
        --left "${species}_trimmed_1.fastq.gz" \
        --right "${species}_trimmed_2.fastq.gz" \
        --SS_lib_type RF \
        --est_method RSEM \
        --aln_method bowtie \
        --trinity_mode \
        --prep_reference \
        --output_dir "rsem_outdir_${species}" \
        --coordsort_bam

    cd ..
done
