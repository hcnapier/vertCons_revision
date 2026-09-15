"""
Generate ESM2 protein embeddings for a species' proteome, in the format UCE
expects to add a new species (a .pt file mapping gene_id -> embedding vector).

UCE's official checkpoints use facebook/esm2_t48_15B_UR50D (embedding dim
5120) to match the model's token dimension — using a smaller ESM2 variant
(e.g. the 650M model, dim 1280) will produce embeddings of the wrong shape
for UCE's pretrained weights. The 15B model needs a GPU with ~40GB+ memory
and takes roughly an hour for a ~20k-gene proteome; if that's not available,
see the note on smaller models at the bottom of this file.

Steps:
  1. Take a proteome FASTA (one sequence per protein, ideally the longest
     isoform per gene — this matches how UCE's own species were processed).
  2. Reduce to one representative protein per gene.
  3. Run each protein through ESM2, mean-pool the per-residue embeddings.
  4. Save as a dict {gene_id: float16 tensor[5120]} — this is the exact
     structure to drop into UCE's model_files/protein_embeddings/ directory.

Usage:
    python proteome_to_esm2_embeddings.py \
        --fasta proteome.fasta \
        --species_name my_species \
        --out_dir ./protein_embeddings \
        --gene_id_regex "gene:(\\S+)"

IMPORTANT — model loading is local-only, by design:
This script loads the model with local_files_only=True, so it will NEVER
attempt a network download. Each array task therefore requires the model to
already be present in the HF cache it points to. Before submitting your
Slurm array job:

  1. Make sure $HF_HOME (or the default ~/.cache/huggingface) points to a
     path that is shared/visible from every compute node — e.g. a group
     project directory, not a node-local scratch disk.
  2. Set that same HF_HOME in your Slurm submission script, so every array
     task reads from the same cache:
         export HF_HOME=/hpc/group/vertgenlab/hailey/software/hf_cache
  3. Confirm the model is already downloaded there (see the earlier
     scan-cache / local_files_only verification steps) BEFORE launching
     the array job — if the cache is missing or incomplete, this script
     will now fail immediately and clearly instead of trying to download.
"""

import argparse
import re
from pathlib import Path

import torch
from tqdm import tqdm
from transformers import AutoTokenizer, AutoModel


ESM2_MODEL_NAME = "facebook/esm2_t48_15B_UR50D"  # dim 5120, matches UCE tokens
# Smaller fallback options (NOT dimension-compatible with UCE's pretrained
# weights as-is — only use these if you're training/fine-tuning your own UCE
# variant rather than using the released checkpoints):
#   "facebook/esm2_t36_3B_UR50D"    (dim 2560)
#   "facebook/esm2_t33_650M_UR50D"  (dim 1280)

MAX_SEQ_LEN = 1022  # ESM2 positional limit is 1024, minus BOS/EOS tokens


def parse_fasta_longest_isoform(fasta_path: str, gene_id_regex: str,
                                 collision_check_regex: str = None) -> dict:
    """
    Parse a FASTA file and keep only the longest protein sequence per gene.

    gene_id_regex should be a regex with one capture group that extracts the
    gene ID from each FASTA header. Ensembl peptide FASTA headers typically
    look like:
        >ENSP00000493376.2 pep primary_assembly:GRCh38:... gene:ENSG00000198888.2 ...
    so a regex like r"gene:(\\S+)" pulls out "ENSG00000198888.2".
    Adjust this for your proteome source (Ensembl, RefSeq, UniProt, etc.).

    collision_check_regex is optional. If provided, it should extract a
    second, more stable identifier from the same header (typically the
    Ensembl gene ID, e.g. r"gene:(\\S+)") so we can detect when multiple
    DISTINCT genes have been grouped under the same gene_id_regex key — most
    relevant when gene_id_regex extracts gene SYMBOLS, since symbols aren't
    guaranteed unique the way Ensembl gene IDs are (shared/informal symbols,
    paralogs, imperfect annotation, etc). Without this check, such
    collisions get silently collapsed to whichever isoform is longest.
    """
    gene_re = re.compile(gene_id_regex)
    collision_re = re.compile(collision_check_regex) if collision_check_regex else None

    sequences = {}          # grouping_key -> longest seq
    collision_ids = {}      # grouping_key -> set of distinct collision-check IDs seen
    gene_id, seq_lines, header = None, [], None

    def flush():
        if header is None:
            return
        seq = "".join(seq_lines)
        m = gene_re.search(header)
        if not m:
            return  # skip entries we can't map to a gene
        gid = m.group(1)
        if gid not in sequences or len(seq) > len(sequences[gid]):
            sequences[gid] = seq

        if collision_re is not None:
            cm = collision_re.search(header)
            if cm:
                collision_ids.setdefault(gid, set()).add(cm.group(1))

    with open(fasta_path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                flush()
                header = line[1:]
                seq_lines = []
            else:
                seq_lines.append(line)
        flush()

    print(f"Parsed {fasta_path}: {len(sequences)} genes "
          f"(longest isoform kept per gene).")

    if collision_re is not None:
        colliding = {k: v for k, v in collision_ids.items() if len(v) > 1}
        if colliding:
            print(f"\nWARNING: {len(colliding)} grouping key(s) matched more than "
                  f"one distinct ID under --collision_check_regex — likely a "
                  f"shared/duplicate gene symbol across different genes. Only "
                  f"the longest isoform overall was kept per key, silently "
                  f"discarding the other gene(s). First 10 examples:")
            for k in list(colliding)[:10]:
                print(f"  '{k}'  <-  {sorted(colliding[k])}")
            print(
                "If this matters for your analysis (e.g. you need every "
                "distinct gene represented, not just one per symbol), "
                "resolve these manually — e.g. by disambiguating symbols "
                "before running this script, or building the embeddings "
                "dict from the collision ID and only relabeling to symbol "
                "afterward.\n"
            )
        else:
            print("No gene-symbol collisions detected against "
                  "--collision_check_regex.")

    return sequences


def embed_proteome(gene_to_seq: dict, model_name: str, device: str,
                    batch_size: int = 1) -> dict:
    """Run each protein through ESM2 and mean-pool per-residue embeddings."""
    print(f"Loading {model_name} from local cache (no download attempted)...")
    tokenizer = AutoTokenizer.from_pretrained(model_name, local_files_only=True)
    model = AutoModel.from_pretrained(
        model_name,
        dtype=torch.float16,
        local_files_only=True,
        low_cpu_mem_usage=True,
        device_map="auto",  # places weights on GPU as loaded; offloads to
                             # CPU RAM automatically if GPU memory is tight,
                             # instead of crashing outright like model.to() does
    )
    model.eval()
    # With device_map="auto" the model may be split across devices, so send
    # inputs to wherever the first layer actually landed rather than
    # assuming everything is on a single --device string.
    input_device = next(model.parameters()).device

    gene_ids = list(gene_to_seq.keys())
    embeddings = {}

    with torch.no_grad():
        for i in tqdm(range(0, len(gene_ids), batch_size), desc="Embedding proteins"):
            batch_genes = gene_ids[i:i + batch_size]
            batch_seqs = [gene_to_seq[g][:MAX_SEQ_LEN] for g in batch_genes]

            tokens = tokenizer(batch_seqs, return_tensors="pt", padding=True,
                                truncation=True, max_length=MAX_SEQ_LEN + 2)
            tokens = {k: v.to(input_device) for k, v in tokens.items()}

            out = model(**tokens)
            hidden = out.last_hidden_state  # (batch, seq_len, embed_dim)

            # Mean-pool over real residues only (exclude padding + BOS/EOS)
            mask = tokens["attention_mask"].unsqueeze(-1).float()
            # zero out BOS (position 0) and EOS (last real token) contributions
            mask[:, 0, :] = 0
            for b, seq in enumerate(batch_seqs):
                eos_pos = min(len(seq) + 1, mask.shape[1] - 1)
                mask[b, eos_pos, :] = 0

            summed = (hidden * mask).sum(dim=1)
            counts = mask.sum(dim=1).clamp(min=1)
            pooled = (summed / counts).to(torch.float16).cpu()

            for gene_id, vec in zip(batch_genes, pooled):
                embeddings[gene_id] = vec

    return embeddings


def main():
    parser = argparse.ArgumentParser(description=__doc__,
                                      formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--fasta", required=True,
                         help="Path to proteome FASTA (protein sequences).")
    parser.add_argument("--species_name", required=True,
                         help="Species name, e.g. 'canis_lupus_familiaris'. "
                              "Used only for the output filename.")
    parser.add_argument("--out_dir", default="./protein_embeddings",
                         help="Directory to save the output .pt file.")
    parser.add_argument("--gene_id_regex", default=r"gene:(\S+)",
                         help=r"Regex with one capture group to pull the gene "
                              r"ID out of each FASTA header. Default matches "
                              r"Ensembl peptide FASTA format: 'gene:(\S+)'.")
    parser.add_argument("--collision_check_regex", default=None,
                         help=r"Optional. A second regex (one capture group) "
                              r"pulling a more stable ID from the same header "
                              r"— typically the Ensembl gene ID, e.g. "
                              r"'gene:(\S+)'. Use this whenever --gene_id_regex "
                              r"groups by something less strict than a unique "
                              r"gene ID, most commonly gene SYMBOL (e.g. "
                              r"'gene_symbol:(\S+)'), to detect and warn about "
                              r"cases where the same symbol is shared by more "
                              r"than one distinct gene.")
    parser.add_argument("--model_name", default=ESM2_MODEL_NAME,
                         help="HuggingFace ESM2 model to use. Defaults to "
                              "the 15B model UCE's released checkpoints use.")
    parser.add_argument("--device", default="cuda" if torch.cuda.is_available() else "cpu",
                         help="Kept for compatibility, but no longer controls "
                              "model placement — device_map='auto' (via the "
                              "'accelerate' package) now decides that "
                              "automatically, including offloading to CPU "
                              "if GPU memory is insufficient.")
    parser.add_argument("--batch_size", type=int, default=1,
                         help="Sequences per forward pass. Keep small (1-4) "
                              "for the 15B model unless you have a lot of GPU memory.")
    args = parser.parse_args()

    gene_to_seq = parse_fasta_longest_isoform(args.fasta, args.gene_id_regex,
                                               args.collision_check_regex)
    if not gene_to_seq:
        raise ValueError("No sequences parsed — check --gene_id_regex against "
                          "your FASTA headers.")

    embeddings = embed_proteome(gene_to_seq, args.model_name, args.device,
                                 args.batch_size)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{args.species_name}_esm2_embeddings.pt"
    torch.save(embeddings, out_path)

    dim = next(iter(embeddings.values())).shape[0]
    print(f"\nSaved {len(embeddings)} gene embeddings (dim={dim}) to {out_path}")
    print(
        "\nNext steps to register this species with UCE:\n"
        "  1. Copy this .pt file into UCE's model_files/protein_embeddings/ directory.\n"
        "  2. Add an entry for this species in species_offsets.pkl (assigns a\n"
        "     token-index range for its genes within UCE's shared vocabulary).\n"
        "  3. Add a row per gene to species_chrom.csv (gene -> chromosome/start\n"
        "     position), which UCE uses to order genes along the genome when\n"
        "     building each cell's input sequence.\n"
        "  4. Pass --species {species_name} to eval_single_anndata.py.\n"
        "  These steps are described in the UCE repo's data_proc/ scripts — "
        "check there for the exact current format, since it can change "
        "between releases."
    )


if __name__ == "__main__":
    main()