from pymol import cmd
import numpy as np

def multi_rmsd_residuewise(reference, others, selection="name CA"):
    ref_atoms = cmd.get_model(f"{reference} and {selection}").atom
    ref_map = {(a.chain, a.resi): a for a in ref_atoms}
    rmsd_data = {key: [] for key in ref_map}

    for model in others:
        other_atoms = cmd.get_model(f"{model} and {selection}").atom
        other_map = {(a.chain, a.resi): a for a in other_atoms}

        for key in rmsd_data:
            if key in other_map:
                a1 = ref_map[key]
                a2 = other_map[key]
                dist = np.linalg.norm(np.array(a1.coord) - np.array(a2.coord))
                rmsd_data[key].append(dist)

    # Compute average RMSD per residue and assign to B-factor
    for (chain, resi), values in rmsd_data.items():
        if values:
            avg_rmsd = np.mean(values)
            cmd.alter(f"{reference} and chain {chain} and resi {resi} and name CA", f"b={avg_rmsd}")

    cmd.rebuild()
    cmd.spectrum("b", "red_blue", reference)
    print("Residue-wise RMSD coloring complete.")
