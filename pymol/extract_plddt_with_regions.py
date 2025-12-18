import os
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser

def extract_plddt(file_path, region_ranges=None):
    """Extracts per-residue pLDDT (from B-factor) and calculates stats."""
    parser = MMCIFParser(QUIET=True) if file_path.endswith(".cif") else PDBParser(QUIET=True)
    structure = parser.get_structure("model", file_path)

    ca_bfactors = []
    region_bfactors = {region: [] for region in region_ranges} if region_ranges else {}

    for atom in structure.get_atoms():
        if atom.get_name() == "CA":
            b = atom.get_bfactor()
            resi = atom.get_parent().get_id()[1]  # residue number

            ca_bfactors.append(b)

            if region_ranges:
                for label, (start, end) in region_ranges.items():
                    if start <= resi <= end:
                        region_bfactors[label].append(b)

    if not ca_bfactors:
        return None

    result = {
        "Mean pLDDT": round(sum(ca_bfactors) / len(ca_bfactors), 2),
        "% pLDDT < 70": f"{round(sum(1 for b in ca_bfactors if b < 70) / len(ca_bfactors) * 100, 1)}%"
    }

    for label, values in region_bfactors.items():
        if values:
            result[f"{label} mean pLDDT"] = round(sum(values) / len(values), 2)
        else:
            result[f"{label} mean pLDDT"] = "NA"

    return result

#CONFIGURATION

# Map protein name → file path
model_files = {
    "A.sparsa": "model_1.cif",
    "A.quinquefasciata": "model_2.cif",
    "C.rubiginosa": "model_3.cif",
    "C.sexpunctata": "model_4.cif",
    "C.alternans": "model_5.cif",
    "C.bullata": "model_6.cif",
    "C.gressoria": "model_7.cif"
}

# Per-protein exon definitions 
region_map = {
    "A.sparsa": {"exon_1": (1, 461), "exon_2_3": (462, 869)},
    "A.quinquefasciata": {"exon_1": (1, 533), "exon_2_3": (534, 876)},
    "C.rubiginosa": {"exon_1": (1, 528), "exon_2_3": (529, 868)},
    "C.sexpunctata": {"exon_1": (1, 631), "exon_2_3": (632, 974)},
    "C.alternans": {"exon_1": (1, 531), "exon_2_3": (532, 872)},
    "C.bullata": {"exon_1": (1, 529), "exon_2_3": (530, 869)},
    "C.gressoria": {"exon_1": (1, 530), "exon_2_3": (531, 870)}
}

#ANALYSIS

results = []
for name, path in model_files.items():
    if os.path.isfile(path):
        regions = region_map.get(name, None)
        stats = extract_plddt(path, regions)
        if stats:
            row = {"Protein": name}
            row.update(stats)
            results.append(row)
    else:
        results.append({"Protein": name, "Mean pLDDT": "File not found"})

# Save results
df = pd.DataFrame(results)
df.to_csv("model_confidence_summary_variable_regions_2.csv", index=False)
print(df)
