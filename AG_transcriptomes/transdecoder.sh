#!/bin/bash
#$ -S /bin/bash
#$ -N transdec
#$ -cwd

source ~/.bashrc

for SPECIES_DIR in */ ; do
  FASTA="${SPECIES_DIR}/Trinity.fasta"
  [[ -f "$FASTA" ]] || continue

  TAG="${SPECIES_DIR%/}"   

  echo "Processing ${TAG}"

  cd "$SPECIES_DIR" || exit 1

  TransDecoder.LongOrfs -t Trinity.fasta -S

  hmmscan \
    --cpu 16 \
    --domtblout "${TAG}_pfam.domtblout" \
    Pfam-A.hmm \
    "${TAG}.transdecoder_dir/longest_orfs.pep"

  TransDecoder.Predict \
    -t Trinity.fasta \
    --retain_pfam_hits "${TAG}_pfam.domtblout" \
    --cpu 16

  cd ..
done

