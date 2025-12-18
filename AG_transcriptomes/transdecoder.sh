#!/bin/bash
#$ -S /bin/bash
#$ -N transdec
#$ -cwd

source ~/.bashrc

for species_dir in */ ; do
  fasta="${species_dir}/Trinity.fasta"
  [[ -f "$fasta" ]] || continue

  tag="${species_dir%/}"   
  
  echo "Processing ${tag}"

  cd "$species_dir" || exit 1

  TransDecoder.LongOrfs -t Trinity.fasta -S

  hmmscan \
    --cpu 16 \
    --domtblout "${tag}_pfam.domtblout" \
    Pfam-A.hmm \
    "${tag}.transdecoder_dir/longest_orfs.pep"

  TransDecoder.Predict \
    -t Trinity.fasta \
    --retain_pfam_hits "${tag}_pfam.domtblout" \
    --cpu 16

  cd ..
done

