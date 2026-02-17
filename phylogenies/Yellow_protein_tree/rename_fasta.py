from Bio import SeqIO
import pandas as pd

# Load name mapping from CSV with semicolon delimiter
name_map = pd.read_csv("fasta_IDs.csv", delimiter=";")
name_dict = dict(zip(name_map["old_name"], name_map["new_name"]))

# Input and output FASTA files
input_fasta = "yellow_filtered_coleopt_hymenopt_chely_bact.fasta"
output_fasta = "yellow_filtered_coleopt_hymenopt_chely_bact_renamed.fasta"

# Process and rename sequences
with open(output_fasta, "w") as output_handle:
    for record in SeqIO.parse(input_fasta, "fasta"):
        if record.id in name_dict:
            record.id = name_dict[record.id]
            record.description = record.id  # optional: overwrite full header
        SeqIO.write(record, output_handle, "fasta")
