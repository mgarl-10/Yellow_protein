#!/bin/bash
#$ -S /bin/bash
#$ -N transdec
#$ -cwd

source ~/.bashrc

for SPECIES_DIR in */ ; do
    FASTA=$(ls "${SPECIES_DIR}"/*_trinity.fasta 2>/dev/null)

    [[ -f "$FASTA" ]] || continue

    BASENAME=$(basename "$FASTA")
    TAG=${BASENAME%_trinity.fasta}

    echo "Processing ${TAG}"

    cd "$SPECIES_DIR" || exit

    # Long ORFs
    TransDecoder.LongOrfs -t "$BASENAME" -S

    # Pfam scan
    hmmscan \
      --cpu 16 \
      --domtblout "${TAG}_pfam.domtblout" \
      Pfam-A.hmm \
      "${TAG}.transdecoder_dir/longest_orfs.pep"

    # Predict coding sequences
    TransDecoder.Predict \
      -t "$BASENAME" \
      --retain_pfam_hits "${TAG}_pfam.domtblout" \
      --cpu 16

    cd ..
done
