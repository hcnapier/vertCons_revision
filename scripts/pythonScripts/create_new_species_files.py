"""
Generate the species-specific files UCE needs for a novel species: a
gene-to-chromosome-position CSV, a combined protein+chromosome token file,
and an offsets pickle. This is a faithful script port of UCE's own
"data_proc/Create New Species Files.ipynb" notebook — same logic, same
output files, just parameterized via argparse so it can run as a Slurm
array job (one species per task), the same way as
proteome_to_esm2_embeddings.py.

Requires as input:
  - The same proteome FASTA you used for the ESM2 embedding step.
  - The ESM2 embeddings .pt file that script produced (must be a dict
    keyed by gene SYMBOL — i.e. was generated with
    --gene_id_regex "gene_symbol:(\\S+)").
  - UCE's own model_files/all_tokens.torch (from the cloned UCE repo),
    which supplies the special tokens every species' token file starts with.

Usage (single species, matches one Slurm array task):
    python create_new_species_files.py \
        --species_name bos_taurus \
        --fasta Bos_taurus.ARS-UCD2.0.pep.all.fa \
        --protein_embeddings bos_taurus_esm2_embeddings.pt \
        --assembly_name ARS-UCD2.0 \
        --taxonomy_id 9913 \
        --all_tokens_path /path/to/UCE/model_files/all_tokens.torch \
        --out_dir ./uce_species_files

Output (all written to --out_dir):
    {species_name}_to_chrom_pos.csv
    {species_name}_pe_tokens.torch
    {species_name}_offsets.pkl

After this runs, you still need to (see UCE's README / notebook for the
current exact procedure):
  1. Add a row to UCE's model_files/new_species_protein_embeddings.csv:
         species_name,full_path_to_protein_embedding_file
  2. Add the species to the dict around line 247 of data_proc/data_utils.py
  3. Pass the printed CHROM_TOKEN_OFFSET value, plus the three output file
     paths above, as arguments to eval_single_anndata.py.
"""

import argparse
import pickle
from pathlib import Path

import pandas as pd
import torch

TOKEN_DIM = 5120  # matches UCE's token dimension (ESM2 15B model output)


def parse_gene_positions(fasta_path: str, species_name: str) -> tuple[pd.DataFrame, list]:
    """
    Parse gene_symbol -> (chromosome, start position) directly from FASTA
    headers. Mirrors the notebook's logic exactly, including its fallback
    order (chromosome: -> primary_assembly: -> scaffold:) and its handling
    of gene symbols that contain colons (e.g. some zebrafish gene names).
    """
    with open(fasta_path) as f:
        proteome_lines = f.readlines()

    gene_symbol_to_location = {}
    gene_symbol_to_chrom = {}
    missing_genes = []

    for line in proteome_lines:
        if not line.startswith(">"):
            continue

        split_line = line.split()
        gene_symbol_tok = [t for t in split_line if t.startswith("gene_symbol")]
        if len(gene_symbol_tok) == 0:
            continue  # no gene_symbol tag at all on this entry — skip, same as notebook

        parts = gene_symbol_tok[0].split(":")
        if len(parts) == 2:
            gene_symbol = parts[1]
        elif len(parts) > 2:
            gene_symbol = ":".join(parts[1:])  # fix for gene names containing colons
        else:
            raise ValueError(f"Unexpected gene_symbol token format: {gene_symbol_tok[0]!r}")

        chrom = None
        chrom_arr = [t for t in split_line if t.startswith("chromosome:")]
        if chrom_arr:
            chrom = chrom_arr[0].replace("chromosome:", "")
        else:
            chrom_arr = [t for t in split_line if t.startswith("primary_assembly:")]
            if chrom_arr:
                chrom = chrom_arr[0].replace("primary_assembly:", "")
            else:
                chrom_arr = [t for t in split_line if t.startswith("scaffold:")]
                if chrom_arr:
                    chrom = chrom_arr[0].replace("scaffold:", "")

        if chrom is not None:
            gene_symbol_to_location[gene_symbol] = chrom.split(":")[2]
            gene_symbol_to_chrom[gene_symbol] = chrom.split(":")[1]
        else:
            missing_genes.append(gene_symbol)

    positional_df = pd.DataFrame()
    positional_df["gene_symbol"] = [gn.upper() for gn in gene_symbol_to_chrom.keys()]
    positional_df["chromosome"] = list(gene_symbol_to_chrom.values())
    positional_df["start"] = list(gene_symbol_to_location.values())
    positional_df = positional_df.sort_values(["chromosome", "start"])
    positional_df["species"] = species_name

    return positional_df, missing_genes


def build_token_file(species_name: str, protein_embeddings_path: str,
                      all_tokens_path: str, taxonomy_id: int,
                      n_uniq_chrom: int) -> tuple[torch.Tensor, dict, int]:
    """
    Combine special tokens + this species' protein embeddings + seeded
    random chromosome tokens into a single token tensor, and compute the
    species' gene offset. Mirrors the notebook's token-generation cell
    exactly, including which offset gets saved where.
    """
    all_pe = torch.load(all_tokens_path)[0:4]  # special tokens at the top
    offset = len(all_pe)  # THIS is what goes into {species}_offsets.pkl —
                           # not the same as CHROM_TOKEN_OFFSET below

    PE = torch.load(protein_embeddings_path)
    pe_stacked = torch.stack(list(PE.values()))
    all_pe = torch.vstack((all_pe, pe_stacked))

    species_to_offsets = {species_name: offset}

    chrom_token_offset = all_pe.shape[0]  # where chromosome tokens start —
                                           # this IS the value you pass as
                                           # --CHROM_TOKEN_OFFSET later

    torch.manual_seed(taxonomy_id)
    chrom_tensors = torch.normal(mean=0, std=1, size=(n_uniq_chrom, TOKEN_DIM))
    all_pe = torch.vstack((all_pe, chrom_tensors))
    all_pe.requires_grad = False

    return all_pe, species_to_offsets, chrom_token_offset


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--species_name", required=True,
                         help="Short species name, e.g. 'bos_taurus'. Used "
                              "as the prefix for all output filenames.")
    parser.add_argument("--fasta", required=True,
                         help="Path to the proteome FASTA (same file used "
                              "for the ESM2 embedding step).")
    parser.add_argument("--protein_embeddings", required=True,
                         help="Path to the ESM2 embeddings .pt file, keyed "
                              "by gene symbol (produced with "
                              "--gene_id_regex 'gene_symbol:(\\S+)').")
    parser.add_argument("--assembly_name", required=True,
                         help="Primary assembly name — for your own sanity "
                              "check that it matches the FASTA; not used in "
                              "the parsing logic itself, but printed "
                              "alongside a header sample so you can confirm "
                              "it lines up before trusting the output.")
    parser.add_argument("--taxonomy_id", type=int, required=True,
                         help="NCBI Taxonomy ID. Determines the random seed "
                              "for chromosome tokens — use the SAME value "
                              "every time you regenerate this species' "
                              "files, or embeddings won't be reproducible.")
    parser.add_argument("--all_tokens_path", required=True,
                         help="Path to UCE's model_files/all_tokens.torch "
                              "from the cloned UCE repo.")
    parser.add_argument("--out_dir", default="./uce_species_files",
                         help="Directory to write the three output files.")
    args = parser.parse_args()

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    # Sanity check: show the user a header line so they can visually
    # confirm --assembly_name matches what's actually in the FASTA.
    with open(args.fasta) as f:
        first_header = next(line for line in f if line.startswith(">"))
    print(f"First FASTA header: {first_header.strip()}")
    print(f"Check that this contains your --assembly_name "
          f"('{args.assembly_name}') in a primary_assembly:/chromosome: tag.\n")

    print(f"Parsing gene positions from {args.fasta} ...")
    positional_df, missing_genes = parse_gene_positions(args.fasta, args.species_name)
    print(f"Mapped {len(positional_df)} genes to chromosome positions.")
    print(f"Genes with a gene_symbol but no chromosome/assembly/scaffold "
          f"tag (excluded): {len(missing_genes)}")

    chrom_pos_path = out_dir / f"{args.species_name}_to_chrom_pos.csv"
    positional_df.to_csv(chrom_pos_path, index=False)
    print(f"Saved {chrom_pos_path}")

    n_uniq_chrom = positional_df["chromosome"].nunique()
    print(f"N_UNIQ_CHROM: {n_uniq_chrom}")
    print(positional_df["chromosome"].value_counts().head(20))

    print(f"\nBuilding token file from {args.protein_embeddings} ...")
    all_pe, species_to_offsets, chrom_token_offset = build_token_file(
        args.species_name, args.protein_embeddings, args.all_tokens_path,
        args.taxonomy_id, n_uniq_chrom,
    )

    pe_tokens_path = out_dir / f"{args.species_name}_pe_tokens.torch"
    torch.save(all_pe, pe_tokens_path)
    print(f"Saved {pe_tokens_path}  (shape: {tuple(all_pe.shape)})")

    offsets_path = out_dir / f"{args.species_name}_offsets.pkl"
    with open(offsets_path, "wb") as f:
        pickle.dump(species_to_offsets, f)
    print(f"Saved {offsets_path}  (contents: {species_to_offsets})")

    print(
        f"\n{'=' * 60}\n"
        f"CHROM_TOKEN_OFFSET: {chrom_token_offset}\n"
        f"{'=' * 60}\n"
        f"Save this number — you need it to run eval_single_anndata.py.\n\n"
        f"Remaining manual steps before you can embed this species:\n"
        f"  1. Add a row to UCE's model_files/new_species_protein_embeddings.csv:\n"
        f"       {args.species_name},{Path(args.protein_embeddings).resolve()}\n"
        f"  2. Add '{args.species_name}' to the species dict around line 247 "
        f"of data_proc/data_utils.py in the UCE repo.\n"
        f"  3. Run:\n"
        f"       accelerate launch eval_single_anndata.py \\\n"
        f"         --adata_path your_data.h5ad \\\n"
        f"         --species {args.species_name} \\\n"
        f"         --CHROM_TOKEN_OFFSET {chrom_token_offset} \\\n"
        f"         --spec_chrom_csv_path {chrom_pos_path.resolve()} \\\n"
        f"         --token_file {pe_tokens_path.resolve()} \\\n"
        f"         --offset_pkl_path {offsets_path.resolve()} \\\n"
        f"         --dir output_directory --model_loc <model_weights.torch>\n"
    )


if __name__ == "__main__":
    main()