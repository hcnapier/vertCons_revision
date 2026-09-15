"""
Register a new species with UCE by adding (or updating) a row in
model_files/new_species_protein_embeddings.csv. This is the only manual
registration step actually required — get_species_to_pe() in
data_proc/data_utils.py already reads this CSV and merges it into its
species dict automatically, so data_utils.py itself does not need editing.

Expected CSV columns (must match exactly, since data_utils.py does
pd.read_csv(...).set_index("species").to_dict()["path"]):
    species,path

Usage:
    python register_uce_species.py \
        --species_name bos_taurus \
        --protein_embeddings_path /path/to/bos_taurus_esm2_embeddings.pt \
        --csv_path /path/to/UCE/model_files/new_species_protein_embeddings.csv

Safe to call repeatedly / from an array job: if the species is already
registered, its path is updated in place rather than duplicated.
"""

import argparse
from pathlib import Path

import pandas as pd


def register_species(csv_path: str, species_name: str, protein_embeddings_path: str):
    csv_path = Path(csv_path)
    protein_embeddings_path = str(Path(protein_embeddings_path).resolve())

    if not Path(protein_embeddings_path).exists():
        raise FileNotFoundError(
            f"--protein_embeddings_path does not exist: {protein_embeddings_path}"
        )

    if csv_path.exists():
        df = pd.read_csv(csv_path)
        if not {"species", "path"}.issubset(df.columns):
            raise ValueError(
                f"{csv_path} exists but doesn't have the expected 'species','path' "
                f"columns (found: {list(df.columns)}). Refusing to modify it "
                f"automatically — check the file by hand."
            )
    else:
        df = pd.DataFrame(columns=["species", "path"])
        csv_path.parent.mkdir(parents=True, exist_ok=True)

    if species_name in df["species"].values:
        old_path = df.loc[df["species"] == species_name, "path"].iloc[0]
        if old_path == protein_embeddings_path:
            print(f"'{species_name}' already registered with the same path — no change.")
            return
        df.loc[df["species"] == species_name, "path"] = protein_embeddings_path
        print(f"Updated '{species_name}': {old_path} -> {protein_embeddings_path}")
    else:
        df = pd.concat(
            [df, pd.DataFrame([{"species": species_name, "path": protein_embeddings_path}])],
            ignore_index=True,
        )
        print(f"Added '{species_name}' -> {protein_embeddings_path}")

    df.to_csv(csv_path, index=False)
    print(f"Saved {csv_path} ({len(df)} total registered species)")


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--species_name", required=True)
    parser.add_argument("--protein_embeddings_path", required=True,
                         help="Path to this species' ESM2 embeddings .pt file "
                              "(the gene-symbol-keyed one).")
    parser.add_argument("--csv_path", required=True,
                         help="Path to UCE's model_files/new_species_protein_embeddings.csv")
    args = parser.parse_args()

    register_species(args.csv_path, args.species_name, args.protein_embeddings_path)


if __name__ == "__main__":
    main()
