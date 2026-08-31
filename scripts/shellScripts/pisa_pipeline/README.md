# BGIseq scRNA-seq processing pipeline (PISA + STAR + DropletUtils + Souporcell)

Reproduces the method described for PRJNA1177647-style BGIseq runs:

- Read 1 (30bp): bases 1-20 = cell barcode (CB), bases 21-30 = UMI
- Read 2 (100bp): transcript sequence
- PISA -> FASTQ+, CB corrected against allow-list (Hamming distance <= 1)
- STAR alignment to species-specific reference
- PISA sam2bam / anno -> BAM with CB/UR/GN tags
- PISA corr -> UMI correction within CB+GN blocks (<= 1bp mismatch) -> UB tag
- PISA count -> raw cell x gene matrix (MEX)
- DropletUtils -> remove empty droplets -> filtered matrix
- Souporcell -> genotype-based clustering to flag ambient/cross-sample contamination

## 0. Software you need installed and on PATH

| Tool | Source |
|---|---|
| PISA | https://github.com/shiquan/PISA |
| STAR | https://github.com/alexdobin/STAR |
| samtools | https://github.com/samtools/samtools |
| sra-tools (prefetch/fasterq-dump) | https://github.com/ncbi/sra-tools |
| R + DropletUtils, Matrix, SingleCellExperiment | Bioconductor |
| souporcell | https://github.com/wheaton5/souporcell (Docker/Singularity recommended) |

## 1. Get the run list and download FASTQs

The BGIseq runs for this BioProject are listed at:
https://www.ncbi.nlm.nih.gov/sra?LinkName=bioproject_sra_all&from_uid=1177647

For each SRR accession:
```bash
prefetch SRR_ACCESSION
fasterq-dump --split-files SRR_ACCESSION -O raw/
gzip raw/SRR_ACCESSION_1.fastq raw/SRR_ACCESSION_2.fastq
```
Confirm which file is R1 (30bp, CB+UMI) vs R2 (100bp, transcript) -- fasterq-dump
ordering matches submission order, but always spot-check read lengths:
```bash
zcat raw/SRR_ACCESSION_1.fastq.gz | head -2 | tail -1 | wc -c   # expect ~31 (30bp + newline)
zcat raw/SRR_ACCESSION_2.fastq.gz | head -2 | tail -1 | wc -c   # expect ~101
```

## 2. Fill in `samples.tsv`

One row per run: `run_id  species  reference_key  fastq_R1  fastq_R2`.
`species` must match a row in `ref/species_refs.tsv` (cow, goat, pig, dog,
rabbit, guinea_pig).

## 3. Get reference genomes + build STAR indices

`ref/species_refs.tsv` lists starting points for each assembly
(ARS-UCD2.0, ARS1.2, Sscrofa11.1, UU_Cfam_GSD_1.0, UM_NZW_1.0, mCavPor4.1).
Verify the exact FASTA/GTF file names on the source site (Ensembl for cow/
goat/pig/guinea pig; NCBI RefSeq/GenBank for dog/rabbit, since those
assemblies aren't on Ensembl), then:

```bash
scripts/01_build_star_index.sh cow  ref/cow.fa  ref/cow.gtf  star_index/cow  16 99
# repeat per species; sjdbOverhang=99 because Read2 = 100bp
```

## 4. Cell barcode allow-list

The method states CB correction uses "the allow list" -- this is the
predefined whitelist of valid 20bp barcodes for the BGI/DNBelab C4 kit used
to generate the libraries. If you have the kit's whitelist file, pass its
path as the last argument to `03_run_all.sh` / `02_process_run.sh`. Without
it, PISA will still extract CB/UMI but cannot correct barcodes against a
known list -- flag this if you don't have the whitelist, since it changes
results at the margins (barcodes with 1 mismatch to an unknown intended
barcode won't be merged in).

## 5. Run everything

```bash
scripts/03_run_all.sh samples.tsv star_index/ gtf/ results/ 16 /path/to/whitelist.txt
```

This produces, per run, `results/<run_id>/raw_gene_expression/` (MEX: 
`barcodes.tsv.gz`, `features.tsv.gz`, `matrix.mtx.gz`).

## 6. Remove empty droplets

```bash
Rscript scripts/04_droplet_filter.R results/<run_id>/raw_gene_expression results/<run_id>/filtered_gene_expression
```

## 7. Souporcell (ambient RNA / cross-sample contamination check)

```bash
gunzip -k results/<run_id>/filtered_gene_expression/barcodes.tsv.gz
scripts/05_souporcell.sh results/<run_id> ref/<species>.fa \
  results/<run_id>/filtered_gene_expression/barcodes.tsv <K> souporcell_out/<run_id> 16
```
Set `<K>` to however many individuals/genotypes were pooled in that specific
run -- check the run's BioSample metadata for this (souporcell needs it as
an input, it isn't inferred automatically).

## 8. Downstream (not scripted here)

- Load each `filtered_gene_expression/` into Seurat, add metadata with
  `AddMetaData()` if you're integrating with the cross-species datasets
  described in the paper (human, cynomolgus macaque, mouse, rat), and apply
  the same marker-based cell-identity checks described in the methods.
- The rat E19.5 sample GSM8016898_19_5_7-PD is explicitly excluded in the
  original study for low cell yield/quality -- exclude the equivalent low-
  quality samples in your own reprocessing based on the same criteria
  (compare cell yield/quality across replicates at each stage).

## Notes / things to verify before trusting results

- **STAR vs PISA's own alignment tutorial**: PISA's public tutorial uses
  minimap2 for simplicity, but explicitly says STAR is a supported
  alternative (used here, matching the paper's stated method, ref 158).
  Default STAR params are used above; the paper doesn't list non-default
  STAR flags, so if you know their exact STARsolo/STAR invocation, adjust
  `02_process_run.sh` accordingly.
- **`-cr` flag in `PISA corr`**: mirrors PISA's official tutorial command
  and defaults to the standard cell/gene-blocked UMI correction; confirm the
  installed PISA version's `PISA corr --help` still exposes it as documented
  the mismatch threshold if you need to force it to exactly 1bp explicitly.
- **Reference URLs** in `ref/species_refs.tsv` are starting points, not
  verified download links -- confirm exact filenames/versions on Ensembl /
  NCBI before downloading, since assembly patch levels change.
