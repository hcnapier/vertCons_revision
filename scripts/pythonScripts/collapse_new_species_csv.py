"""
Collapse per-species UCE registration files into one CSV.

This script reads all *.csv files from a temporary directory and combines
them into a single CSV with the columns:

    species,path

Usage:

    python collapse_new_species_csv.py \
        --temp_dir /path/to/uce_species_temp \
        --csv_path /path/to/UCE/model_files/new_species_protein_embeddings.csv

If the same species occurs in multiple temporary files, the last occurrence
is retained.

The temporary files are NOT deleted unless --clean is specified.
"""

import argparse
from pathlib import Path

import pandas as pd


def collapse_species(temp_dir: str, csv_path: str, clean: bool = False):
    temp_dir = Path(temp_dir)
    csv_path = Path(csv_path)

    if not temp_dir.exists():
        raise FileNotFoundError(
            f"Temporary directory does not exist: {temp_dir}"
        )

    files = sorted(temp_dir.glob("*.csv"))

    if not files:
        raise FileNotFoundError(
            f"No *.csv registration files found in {temp_dir}"
        )

    dataframes = []

    for file in files:
        df = pd.read_csv(file)

        if not {"species", "path"}.issubset(df.columns):
            raise ValueError(
                f"{file} does not have the expected 'species','path' "
                f"columns. Found: {list(df.columns)}"
            )

        dataframes.append(df[["species", "path"]])

    # Combine all per-species files.
    combined = pd.concat(dataframes, ignore_index=True)

    # Remove duplicate species.
    # Since files are sorted above, keep="last" gives deterministic behavior.
    duplicate_count = combined["species"].duplicated(keep=False).sum()

    if duplicate_count:
        print(
            f"Warning: found {duplicate_count} rows belonging to "
            f"duplicate species names. Keeping the last occurrence."
        )

    combined = combined.drop_duplicates(
        subset="species",
        keep="last",
    )

    # Sort alphabetically by species for a stable output file.
    combined = combined.sort_values("species").reset_index(drop=True)

    # Make sure the output directory exists.
    csv_path.parent.mkdir(parents=True, exist_ok=True)

    combined.to_csv(csv_path, index=False)

    print(
        f"Collapsed {len(files)} temporary files into:\n"
        f"    {csv_path}"
    )
    print(
        f"Total registered species: {len(combined)}"
    )

    if clean:
        for file in files:
            file.unlink()

        print(f"Removed {len(files)} temporary files from {temp_dir}")


def main():
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )

    parser.add_argument(
        "--temp_dir",
        required=True,
        help="Directory containing the per-species CSV files.",
    )

    parser.add_argument(
        "--csv_path",
        required=True,
        help="Path for the final combined CSV.",
    )

    parser.add_argument(
        "--clean",
        action="store_true",
        help="Delete the temporary files after successful aggregation.",
    )

    args = parser.parse_args()

    collapse_species(
        args.temp_dir,
        args.csv_path,
        clean=args.clean,
    )


if __name__ == "__main__":
    main()
