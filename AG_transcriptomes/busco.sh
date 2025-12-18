#!/bin/bash
#$ -S /bin/bash
#$ -N busco_all
#$ -cwd

source ~/.bashrc
conda activate busco

lineage=/ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/busco/busco_downloads/lineages/endopterygota_odb10
downloads=/ebio/ag-salem/projects/CassidinaeGenomics/code/minconda3/envs/busco/busco_downloads
threads=64

for species_dir in */ ; do
    fasta="${species_dir}/Trinity.fasta"

    [[ -f "$fasta" ]] || continue

    species_name="${species_dir%/}"
    out="busco_${species_name}_transcriptome"

    echo "Running BUSCO for ${species_name}"

    busco \
      -m transcriptome \
      -i "$fasta" \
      -o "$out" \
      -l "$lineage" \
      -f \
      --offline \
      --download_path "$downloads" \
      -c "$threads"
done

