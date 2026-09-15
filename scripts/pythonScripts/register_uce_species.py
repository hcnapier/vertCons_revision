"""
Register a new species with UCE by writing a one-row CSV file to a
temporary directory.

This script is designed to be safe to call concurrently from a batch
array job. Each species writes its own independent CSV file, so array
tasks never modify the same file.

The individual per-species files can later be collapsed into the final
model_files/new_species_protein_embeddings.csv using
collapse_uce_species.py.

Expected CSV columns:
    species,path

Usage:
    python register_uce_species.py \
        --species_name Bos_taurus \
        --protein_embeddings_path /path/to/Bos_taurus_esm2_embeddings.pt \
        --temp_dir /path/to/uce_species_temp

This creates:
    /path/to/uce_species_temp/Bos_taurus.csv

containing:
    species,path
    Bos_taurus,/path/to/Bos_taurus_esm2_embeddings.pt

Safe to call repeatedly / from an array job: the file for this species
is simply overwritten with the current path.
"""

import argparse
from pathlib import Path


# Species UCE already ships pretrained embeddings for, per
# get_species_to_pe() in data_proc/data_utils.py.
#
# The keys are UCE's internal names and the values are accepted Latin-name
# equivalents. This prevents accidentally overriding UCE's native
# pretrained embeddings when the species is supplied using a different
# naming convention.
BUILTIN_SPECIES = {
    "human": {
        "human",
        "homo_sapiens",
    },
    "mouse": {
        "mouse",
        "mus_musculus",
    },
    "frog": {
        "frog",
    },
    "zebrafish": {
        "zebrafish",
        "danio_rerio",
    },
    "mouse_lemur": {
        "mouse_lemur",
        "microcebus_murinus",
    },
    "pig": {
        "pig",
        "sus_scrofa",
    },
    "macaca_fascicularis": {
        "macaca_fascicularis",
    },
    "macaca_mulatta": {
        "macaca_mulatta",
    },
}


def is_builtin_species(species_name: str) -> bool:
    """
    Return True if species_name corresponds to one of UCE's built-in
    species, allowing for differences in capitalization.
    """
    normalized_name = species_name.strip().lower()

    return any(
        normalized_name in accepted_names
        for accepted_names in BUILTIN_SPECIES.values()
    )


def register_species(
    temp_dir: str,
    species_name: str,
    protein_embeddings_path: str,
    force: bool = False,
):
    if is_builtin_species(species_name) and not force:
        print(
            f"'{species_name}' corresponds to one of UCE's natively "
            f"supported species — skipping registration so UCE's own "
            f"pretrained embeddings aren't overridden. "
            f"(Pass --force if you intentionally want to replace them.)"
        )
        return

    temp_dir = Path(temp_dir)
    protein_embeddings_path = str(Path(protein_embeddings_path).resolve())

    if not Path(protein_embeddings_path).exists():
        raise FileNotFoundError(
            f"--protein_embeddings_path does not exist: "
            f"{protein_embeddings_path}"
        )

    # Create the temporary directory if it doesn't already exist.
    temp_dir.mkdir(parents=True, exist_ok=True)

    # Each species gets its own file, so concurrent array jobs never
    # read/write the same file.
    output_file = temp_dir / f"{species_name}.csv"

    # Write exactly one data row plus the header.
    with open(output_file, "w") as f:
        f.write("species,path\n")
        f.write(f"{species_name},{protein_embeddings_path}\n")

    print(
        f"Registered '{species_name}' -> {protein_embeddings_path}"
    )
    print(f"Wrote {output_file}")


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
            "Directory where the per-species CSV registration file "
            "will be written."
        ),
    )

    parser.add_argument(
        "--force",
        action="store_true",
        help=(
            "Register even if --species_name matches one of UCE's "
            "built-in species, overriding its pretrained embeddings. "
            "Off by default."
        ),
    )

    args = parser.parse_args()

    register_species(
        args.temp_dir,
        args.species_name,
        args.protein_embeddings_path,
        force=args.force,
    )


if __name__ == "__main__":
    main()
