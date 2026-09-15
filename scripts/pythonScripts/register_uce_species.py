```python
"""
Register a new species with UCE by writing a one-line registration file
to a temporary directory.

This script is designed to be safe to run concurrently as a batch array.
Each array task writes its own file, so no shared CSV is modified.

The temporary files can later be collapsed into the final CSV using
collapse_uce_species.py.

Expected final CSV columns:

    species,path

Usage:

    python register_uce_species.py \
        --species_name bos_taurus \
        --protein_embeddings_path /path/to/bos_taurus_esm2_embeddings.pt \
        --temp_dir /path/to/uce_species_temp

This creates:

    /path/to/uce_species_temp/bos_taurus.csv

containing:

    species,path
    bos_taurus,/path/to/bos_taurus_esm2_embeddings.pt
"""

import argparse
from pathlib import Path


def register_species(
    temp_dir: str,
    species_name: str,
    protein_embeddings_path: str,
):
    temp_dir = Path(temp_dir)
    protein_embeddings_path = str(Path(protein_embeddings_path).resolve())

    if not Path(protein_embeddings_path).exists():
        raise FileNotFoundError(
            f"--protein_embeddings_path does not exist: "
            f"{protein_embeddings_path}"
        )

    temp_dir.mkdir(parents=True, exist_ok=True)

    # Use the species name as the filename. Replace characters that could
    # cause problems in filenames.
    safe_species_name = (
        species_name.replace("/", "_")
        .replace("\\", "_")
        .replace(" ", "_")
    )

    output_file = temp_dir / f"{safe_species_name}.csv"

    # Each file contains exactly one registered species.
    with open(output_file, "w") as f:
        f.write("species,path\n")
        f.write(f"{species_name},{protein_embeddings_path}\n")

    print(f"Wrote '{species_name}' -> {output_file}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--species_name",
        required=True,
    )

    parser.add_argument(
        "--protein_embeddings_path",
        required=True,
        help=(
            "Path to this species' ESM2 embeddings .pt file "
            "(the gene-symbol-keyed one)."
        ),
    )

    parser.add_argument(
        "--temp_dir",
        required=True,
        help=(
            "Directory where the per-species registration files "
            "will be written."
        ),
    )

    args = parser.parse_args()

    register_species(
        args.temp_dir,
        args.species_name,
        args.protein_embeddings_path,
    )


if __name__ == "__main__":
    main()
```
