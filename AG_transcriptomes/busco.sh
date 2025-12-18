#!/bin/bash
#$ -S /bin/bash
#$ -N busco_all
#$ -cwd

source ~/.bashrc
conda activate busco

LINEAGE=/ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/busco/busco_downloads/lineages/endopterygota_odb10
DOWNLOADS=/ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/busco/busco_downloads
THREADS=64

for SPECIES_DIR in */ ; do
    FASTA="${SPECIES_DIR}/Trinity.fasta"

    [[ -f "$FASTA" ]] || continue

    SPECIES_NAME="${SPECIES_DIR%/}"
    OUT="busco_${SPECIES_NAME}_transcriptome"

    echo "Running BUSCO for ${SPECIES_NAME}"

    busco \
      -m transcriptome \
      -i "$FASTA" \
      -o "$OUT" \
      -l "$LINEAGE" \
      -f \
      --offline \
      --download_path "$DOWNLOADS" \
      -c "$THREADS"
done

