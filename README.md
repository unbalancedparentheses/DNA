# Genetic Health Analysis Pipeline

Turn raw DNA data into actionable health reports. Works with **23andMe exports** (instant) or **whole-genome sequencing FASTQ files** (full pipeline: QC, alignment, variant calling, analysis).

The pipeline cross-references your genome against ClinVar (340K+ clinical variants), PharmGKB (drug-gene interactions), and a curated database of ~210+ lifestyle/health SNPs to produce a single interactive HTML report with SVG charts, searchable findings, and a print-optimized doctor card.

## Why This Exists

Consumer genetic testing services give you raw data but limited interpretation. Clinical-grade analysis is expensive and often focused on a narrow set of conditions. This pipeline bridges the gap:

- **Comprehensive coverage** -- Analyzes ~210+ curated lifestyle/health variants across 16 categories (drug metabolism, methylation, nutrition, fitness, cardiovascular, sleep, blood type, bitter taste, autoimmune, and more) plus 340K+ ClinVar clinical variants in a single run.
- **Blood type prediction** -- Predicts ABO blood group and Rh factor from proxy SNPs (rs505922, rs8176746, rs590787).
- **Mitochondrial haplogroup** -- Estimates maternal lineage haplogroup from ~25 MT-defining SNPs, covering major lineages (H, J, K, T, U, L, A, B, C, D, and more).
- **Star allele calling** -- CPIC-style pharmacogenomic star allele calling for CYP2C19, CYP2C9, and CYP2D6 with metabolizer phenotype classification (poor/intermediate/normal/rapid/ultrarapid).
- **Data quality metrics** -- Call rate, heterozygosity rate, chromosome coverage, sex inference, and no-call counting from raw genome data.
- **Population frequency annotations** -- Key clinically relevant SNPs include allele frequencies across 5 superpopulations (EUR, AFR, EAS, SAS, AMR) for context.
- **Ancestry-aware** -- Estimates superpopulation proportions from 55 Ancestry-Informative Markers and flags findings that may not apply to your genetic background.
- **Polygenic risk scores** -- Calculates combined risk across multiple SNPs for 5 major conditions (type 2 diabetes, coronary artery disease, hypertension, breast cancer, macular degeneration) using published GWAS effect sizes.
- **Gene-gene interactions** -- Detects 10 known epistatic interactions (e.g., MTHFR + COMT methylation burden) where combined effects differ from individual variants.
- **Drug interactions** -- Merges PharmGKB pharmacogenomic annotations with ClinVar drug-response variants into a single drug interaction card.
- **Interactive exploration** -- The HTML report includes searchable/filterable findings, sortable tables, SVG donut and gauge charts, collapsible sections, and CSV export -- all with zero external dependencies.
- **Works with raw sequencing data** -- The WGS pipeline takes FASTQ files from any sequencer and runs them through the same analysis, so you don't need 23andMe.
- **Fully local** -- All analysis runs on your machine. No data leaves your computer (except optional Ensembl API calls for rsID annotation in the WGS pipeline).
- **Zero runtime dependencies** -- Pure Python with only `pytest` for testing. Bioinformatics tools (minimap2, samtools, bcftools) are provided via Nix for the WGS pipeline.

## Quick Start

### Option A: 23andMe Data (Recommended for most users)

```bash
# 1. Clone and enter the project
git clone <repo-url> && cd DNA

# 2. Enter the Nix development shell (provides Python + tools)
nix develop

# 3. Place your 23andMe raw data export in data/
cp ~/Downloads/genome_Your_Name_*.txt data/genome.txt

# 4. Run the analysis
python -m genetic_health --name "Your Name"

# 5. Open report
open reports/GENETIC_HEALTH_REPORT.html
```

### Option B: Whole-Genome Sequencing (FASTQ)

```bash
# 1. One-time setup: download GRCh37 reference genome (~900MB download, ~10GB disk)
make setup

# 2. Run the full pipeline (FASTQ -> QC -> align -> call -> analyze)
make pipeline FASTQ=/path/to/reads.fastq NAME="Your Name"

# Or with more control:
python -m genetic_health.wgs_pipeline /path/to/reads.fastq \
    --name "Your Name" --threads 12
```

### Option C: Nix Run (No Clone Needed)

```bash
nix run github:<owner>/DNA -- /path/to/reads.fastq --name "Your Name"
```

## Requirements

- **[Nix](https://nixos.org/download/)** -- Provides Python 3.10+, minimap2, samtools, bcftools, and htslib in a reproducible shell. No system-wide installation needed.
- **fastp** (optional, for WGS QC step) -- Not in nixpkgs for Apple Silicon. Install via `brew install fastp` or skip with `--skip-qc`.

### Data Files

Place these in the `data/` directory:

| File | Source | Required For |
|------|--------|--------------|
| `genome.txt` | [23andMe raw data download](https://you.23andme.com/) or WGS pipeline output | All analysis |
| `clinvar_alleles.tsv` | Auto-downloaded via `make update-data` or [ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/) | Disease risk analysis |
| `clinical_annotations.tsv` | [PharmGKB](https://www.pharmgkb.org/downloads) (free account) | Drug-gene interactions |
| `clinical_ann_alleles.tsv` | [PharmGKB](https://www.pharmgkb.org/downloads) (free account) | Drug-gene interactions |

If ClinVar or PharmGKB files are missing, those analysis sections are skipped gracefully -- the pipeline still runs.

## Output

The pipeline generates a single interactive report in `reports/`:

### `GENETIC_HEALTH_REPORT.html`

A self-contained HTML file with 24 sections -- no external dependencies, works offline:

- **ELI5 summary** -- plain-language overview of key findings
- **Data quality dashboard** -- call rate badge, chromosome coverage bar chart, sex inference
- **Research-backed insights** -- multi-gene narratives and gene-specific findings with references
- **APOE haplotype** -- Alzheimer's risk context with visual display
- **Personalized recommendations** -- prioritized actions, drug card, monitoring schedule, specialist referrals
- **SVG dashboard** with impact bar chart, category donut chart, and metabolism gauges
- **Ancestry donut chart** showing superpopulation proportions from ~55 AIMs
- **Blood type & traits** -- ABO/Rh prediction, eye/hair color, earwax type, freckling
- **Mitochondrial haplogroup** -- maternal lineage description
- **PRS gauge visualizations** for 8 conditions (green/blue/orange/red zones)
- **Pharmacogenomic star alleles** -- metabolizer gauges for 6 genes (CYP2C19, CYP2C9, CYP2D6, DPYD, TPMT, UGT1A1)
- **ACMG secondary findings** -- 81 medically actionable genes screened
- **Gene-gene interactions** (epistasis) -- 10 multi-gene interaction models
- **Carrier screening** -- organized by disease system with reproductive context
- **All lifestyle findings** (~210+ SNPs) with population frequency badges, database links, and paper references
- **Drug-gene interactions** -- PharmGKB annotations
- **Print-optimized doctor card**
- **Search box** -- real-time filtering across all findings and table rows
- **Sortable tables** -- click any column header to sort (numeric-aware, toggle asc/desc)
- **CSV export** -- copy findings to clipboard as CSV
- **Collapsible sections** -- expand/collapse each analysis category
- **Dark mode** -- respects system preference
- **Database links** -- every rsID links to dbSNP, ClinVar, SNPedia, and PharmGKB

An optional **PDF** can be generated with `--pdf`.

## Ancestry Estimation

Estimates superpopulation proportions from ~55 Ancestry-Informative Markers (AIMs) using maximum-likelihood scoring with softmax normalization.

**Superpopulations** (matching 1000 Genomes Project):
- EUR (European)
- AFR (African)
- EAS (East Asian)
- SAS (South Asian)
- AMR (Admixed American)

**Key markers include:** SLC24A5, SLC45A2, HERC2 (pigmentation), DARC (malaria resistance), EDAR (hair morphology), ALDH2, ADH1B (alcohol metabolism), MCM6/LCT (lactose tolerance), OCA2, TYR, ABCC11, CYP3A5, and more.

**Confidence levels:** High (40+ markers found), Moderate (20-39), Low (<20).

**Population-specific warnings:** When a finding is primarily relevant to a specific population (e.g., ALDH2 reduced function is ~30-40% in East Asians but rare elsewhere), the reports flag it with context.

## Polygenic Risk Scores

Calculates combined genetic risk across multiple SNPs for 5 conditions using published GWAS log(OR) effect sizes:

| Condition | Reference | SNPs |
|-----------|-----------|------|
| Type 2 Diabetes | Mahajan et al. 2018 | ~25-30 |
| Coronary Artery Disease | Nikpay 2015, Aragam 2022 | ~25-35 |
| Hypertension | Evangelou et al. 2018 | ~20-30 |
| Breast Cancer | Michailidou et al. 2017 | ~20-30 |
| Age-Related Macular Degeneration | Fritsche et al. 2016 | ~15-20 |

**Risk categories:** Low (<20th percentile), Average (20-80th), Elevated (80-95th), High (>95th).

**Ancestry caveat:** PRS models are primarily validated in European populations. When estimated non-European ancestry exceeds 40%, results are flagged as potentially less applicable.

## Gene-Gene Interactions (Epistasis)

Evaluates 10 well-characterized multi-gene interactions where the combined effect differs from individual variants:

- **MTHFR + COMT** -- Methylation-catecholamine dual burden
- **CYP1A2 + ADORA2A** -- Caffeine sensitivity
- **AGT + ACE + ADD1** -- Blood pressure multi-gene risk
- **MTHFR + MTRR** -- Folate cycle compound deficiency
- **COMT + MAO-A** -- Neurotransmitter clearance
- **CYP2D6 + CYP2C19** -- Multi-enzyme drug metabolism
- **Factor V + Prothrombin** -- Clotting cascade compound risk
- **VDR + CYP2R1** -- Vitamin D pathway
- **APOE + CETP** -- Lipid metabolism
- **ACTN3 + COL5A1** -- Athletic injury risk

Each detected interaction includes mechanism explanation, risk level, and specific actionable recommendations.

## Blood Type Prediction

Predicts ABO blood group and Rh factor from proxy SNPs available on 23andMe arrays:

| SNP | Gene | Role |
|-----|------|------|
| rs505922 | ABO | O-type proxy (T allele = O likelihood, C = A or B) |
| rs8176746 | ABO | B antigen (T allele = B antigen present) |
| rs590787 | RHD | Rh factor proxy (TT = Rh-negative) |

```python
from genetic_health.blood_type import predict_blood_type
result = predict_blood_type(genome_by_rsid)
# result = {blood_type: "A+", abo: "A", rh: "+", confidence: "high", details: [...]}
```

**Note:** 23andMe cannot directly genotype the ABO gene's defining indel (rs8176719), so rs505922 is used as a proxy. Results are accurate for most cases but may be ambiguous for A vs AB. Confidence is reported as high/moderate/low depending on how many proxy SNPs are available.

## Mitochondrial Haplogroup

Estimates maternal lineage haplogroup from ~25 mitochondrial-DNA-defining SNPs. The decision tree covers major global haplogroups:

- **African:** L0, L1, L2, L3
- **European:** H, H1, H2, J, J1, T, T2, K, U, U5, V, W, X
- **East Asian / Native American:** A, B, C, D, G
- **South/West Asian:** R, N, M

```python
from genetic_health.mt_haplogroup import estimate_mt_haplogroup
result = estimate_mt_haplogroup(genome_by_rsid)
# result = {haplogroup: "H", description: "Most common European haplogroup (~40%)",
#           confidence: "high", markers_found: 8, markers_tested: 25,
#           lineage: "European maternal", details: [...]}
```

The most specific matching haplogroup wins (e.g., H1 overrides H). Confidence depends on the number of MT SNPs present in the data.

## Pharmacogenomic Star Alleles

CPIC-style star allele calling for three key pharmacogenes:

| Gene | Star Alleles | Defining SNPs |
|------|-------------|---------------|
| CYP2C19 | \*1 (normal), \*2 (rs4244285), \*3 (rs4986893), \*17 (rs12248560) | 3 SNPs |
| CYP2C9 | \*1 (normal), \*2 (rs1799853), \*3 (rs1057910) | 2 SNPs |
| CYP2D6 | \*1 (normal), \*4 (rs3892097), \*10 (rs1065852) | 2 SNPs |

```python
from genetic_health.star_alleles import call_star_alleles
results = call_star_alleles(genome_by_rsid)
# results["CYP2C19"] = {gene: "CYP2C19", diplotype: "*1/*2",
#                        phenotype: "Intermediate Metabolizer", ...}
```

**Metabolizer phenotypes:** Poor, Intermediate, Normal, Rapid, Ultrarapid -- derived from the diplotype's functional impact on enzyme activity.

**Limitation:** 23andMe cannot detect CYP2D6 copy number variants (gene deletions = \*5, duplications = \*1xN). CYP2D6 results always include a caveat about this.

## Data Quality Metrics

Computed from the raw genome file and loaded data:

```python
from genetic_health.quality_metrics import compute_quality_metrics
metrics = compute_quality_metrics(genome_by_rsid, genome_path)
```

| Metric | Description |
|--------|-------------|
| `total_snps` | Total loaded SNPs |
| `no_call_count` | Lines with `--` genotype in raw file |
| `call_rate` | total_snps / (total_snps + no_call_count) |
| `chromosomes` | Dict of chr -> SNP count |
| `has_mt` | Whether MT chromosome SNPs are present |
| `has_y` | Whether Y chromosome SNPs are present (indicates biological male) |
| `het_rate` | Heterozygosity rate across autosomal SNPs |

## WGS Pipeline

For whole-genome sequencing data, the pipeline processes FASTQ files through 5 steps:

```
FASTQ -> [1. QC (fastp)] -> [2. Align (minimap2)] -> [3. Call (bcftools)]
      -> [4. Convert to 23andMe format] -> [5. Full health analysis]
```

**Targeted calling** (default): Only calls variants at ClinVar + SNP database positions -- dramatically faster than genome-wide calling with no loss of relevant information.

**Full calling** (`--full` flag): Genome-wide variant calling for completeness.

The rsID annotation step queries the Ensembl GRCh37 API in parallel with alignment, so it adds no wall-clock time.

```bash
# Basic usage
python -m genetic_health.wgs_pipeline /path/to/reads.fastq

# All options
python -m genetic_health.wgs_pipeline /path/to/reads.fastq \
    --name "Subject Name" \
    --threads 12 \
    --full \               # Genome-wide calling (slower)
    --skip-qc \            # Skip fastp QC
    --skip-analysis \      # Stop after VCF conversion
    --keep-intermediates \ # Keep BAM and temp files
    --output custom.txt    # Custom output path
```

## PDF Export

Export the interactive HTML report to PDF:

```bash
python -m genetic_health --name "Your Name" --pdf
```

Tries Chrome/Chromium headless first, falls back to weasyprint, then offers browser-based print.

## Keeping Data Up to Date

### ClinVar (recommended: quarterly)

```bash
# Automated download + processing
make update-data
# Or directly:
python -m genetic_health.update_data clinvar
```

Downloads `variant_summary.txt.gz` from NCBI FTP (~80MB), filters to GRCh37 assembly, and creates `clinvar_alleles.tsv`. Version metadata is tracked in `data/data_versions.json`.

### PharmGKB (recommended: quarterly)

PharmGKB requires a free account for downloads:

1. Download from https://www.pharmgkb.org/downloads
2. Place `clinical_annotations.tsv` and `clinical_ann_alleles.tsv` in `data/`
3. Validate: `make validate-pharmgkb`

### Check Data Status

```bash
make data-status
# Shows: data versions, file sizes, last update dates
```

## Running for Multiple People

```bash
# Run for each person
python -m genetic_health data/genome_mom.txt --name "Mom"

# Rename output before running the next person
mv reports/GENETIC_HEALTH_REPORT.html reports/GENETIC_HEALTH_REPORT_MOM.html

python -m genetic_health data/genome_dad.txt --name "Dad"
```

## Make Targets

```
make help              Show all available targets
make setup             Download reference genome + create indexes (~10GB disk)
make fastp             Install fastp via Homebrew (Apple Silicon)
make pipeline          Full WGS pipeline: FASTQ -> reports
make analysis          Run analysis on existing data/genome.txt
make test              Run the full test suite
make update-data       Download latest ClinVar data
make validate-pharmgkb Validate PharmGKB files after manual download
make data-status       Show current data versions and file status
make clean             Remove intermediate files (keep BAM/VCF)
make clean-all         Remove all generated files
make shell             Enter Nix development shell
```

## Testing

```bash
make test
# Or: python -m pytest tests/ -v
```

355 tests covering:

- Genome file parsing and loading
- SNP database validation (all ~210+ entries)
- VCF conversion logic
- ClinVar disease analysis (indel filtering, zygosity, significance mapping)
- Ancestry marker validation and scoring
- PRS model validation and calculation
- Epistasis model validation and interaction detection
- Blood type prediction (ABO + Rh inference, missing SNP handling)
- Quality metrics (call rate, heterozygosity, chromosome coverage, sex inference)
- MT haplogroup estimation (decision tree, confidence levels, lineage assignment)
- Star allele calling (diplotype determination, phenotype mapping, compound heterozygotes)
- Recommendation generation from findings
- HTML report generation (inline formatting, tables, SVG generators)
- Enhanced HTML sections (ancestry, PRS, epistasis, blood type, star allele builders)
- Data update processing and metadata
- Path resolution and project structure

All tests use synthetic data -- no genome files or database downloads needed.

## Project Structure

```
DNA/
├── README.md                         # This file
├── CLAUDE.md                         # Development reference (AI assistant context)
├── Makefile                          # Build targets
├── flake.nix                         # Nix flake (reproducible dev environment)
├── pyproject.toml                    # Python project config
│
├── genetic_health/                   # Main Python package
│   ├── __init__.py                   # Public API exports
│   ├── __main__.py                   # python -m genetic_health entry point
│   ├── config.py                     # Centralized paths (BASE_DIR, DATA_DIR, etc.)
│   ├── loading.py                    # Genome and PharmGKB file loaders
│   ├── analysis.py                   # SNP analysis + ClinVar disease scanning
│   ├── snp_database.py              # ~210+ curated lifestyle/health SNP interpretations
│   ├── clinical_context.py          # Detailed clinical interpretations per gene/status
│   ├── ancestry.py                  # 55 AIMs + maximum-likelihood ancestry estimation
│   ├── prs.py                       # Polygenic risk score models (5 conditions)
│   ├── epistasis.py                 # Gene-gene interaction detection (10 models)
│   ├── recommendations.py          # Actionable recommendations from findings
│   ├── blood_type.py               # ABO + Rh blood type prediction from proxy SNPs
│   ├── quality_metrics.py          # Call rate, het rate, chromosome coverage, sex inference
│   ├── mt_haplogroup.py            # Mitochondrial haplogroup (maternal lineage) estimation
│   ├── star_alleles.py             # CPIC-style CYP star allele calling + metabolizer phenotypes
│   ├── pipeline.py                  # Main analysis orchestrator + CLI
│   ├── wgs_pipeline.py             # FASTQ -> alignment -> variant calling -> analysis
│   ├── update_data.py              # ClinVar auto-download + PharmGKB validation
│   └── reports/
│       ├── __init__.py              # Re-exports report generators
│       ├── enhanced_html.py         # Interactive HTML report (SVG, search, sort)
│       └── pdf_export.py           # PDF export (Chrome headless / weasyprint)
│
├── scripts/
│   └── setup_reference.sh           # GRCh37 reference genome download + indexing
│
├── tests/                            # 355 tests (pytest)
│   ├── conftest.py
│   ├── test_genome_loading.py
│   ├── test_snp_analysis.py
│   ├── test_vcf_conversion.py
│   ├── test_disease_analysis.py
│   ├── test_path_resolution.py
│   ├── test_html_reports.py
│   ├── test_ancestry.py
│   ├── test_prs.py
│   ├── test_epistasis.py
│   ├── test_update_data.py
│   ├── test_blood_type.py           # Blood type prediction tests
│   ├── test_quality_metrics.py      # Data quality metrics tests
│   ├── test_mt_haplogroup.py        # MT haplogroup estimation tests
│   ├── test_star_alleles.py         # Star allele calling tests
│   └── test_recommendations.py      # Recommendation generation tests
│
├── data/                             # Input data (gitignored except schema)
│   ├── genome.txt                    # 23andMe raw data or WGS pipeline output
│   ├── clinvar_alleles.tsv           # ClinVar database (~341K variants)
│   ├── clinical_annotations.tsv      # PharmGKB annotations
│   ├── clinical_ann_alleles.tsv      # PharmGKB allele data
│   └── data_versions.json            # Auto-generated version metadata
│
├── reference/                        # Reference genome (gitignored, ~10GB)
│   ├── human_g1k_v37.fasta           # GRCh37 reference
│   ├── human_g1k_v37.fasta.fai       # samtools index
│   └── human_g1k_v37.mmi             # minimap2 index
│
├── wgs_work/                         # WGS intermediate files (gitignored)
│   ├── aligned.bam                   # Sorted BAM
│   ├── variants.vcf.gz              # Called variants
│   └── target_regions.bed            # Targeted calling positions
│
└── reports/                          # Generated output (gitignored)
    ├── GENETIC_HEALTH_REPORT.html
    └── comprehensive_results.json
```

## Interpretation Guide

### Impact Magnitude (0-6)

| Magnitude | Meaning | Example |
|-----------|---------|---------|
| 0 | Informational | Common variant, typical function |
| 1 | Low impact | Minor effect on metabolism |
| 2 | Moderate | Worth noting, may guide choices |
| 3 | High | Actionable, affects health decisions |
| 4-6 | Very high | Requires clinical attention |

### ClinVar Confidence (Gold Stars)

| Stars | Review Status |
|-------|---------------|
| 4 | Practice guideline or expert panel reviewed |
| 3 | Multiple submitters, no conflicts |
| 2 | Multiple submitters with conflicts, or single submitter with criteria |
| 1 | Single submitter |
| 0 | No assertion criteria provided |

### Zygosity

| Zygosity | Inheritance | Implication |
|----------|-------------|-------------|
| Homozygous | Any | Both copies carry variant -- full effect |
| Heterozygous | Dominant | One copy sufficient -- may be affected |
| Heterozygous | Recessive | Carrier only -- reproductive implications, typically unaffected |

## Limitations

1. **Not a clinical diagnosis.** This is an informational tool. Consult a genetic counselor or physician for clinical decisions.
2. **Population bias.** Most GWAS studies and ClinVar submissions are based on European populations. Risk scores and variant classifications may be less accurate for other ancestries.
3. **Incomplete penetrance.** Having a pathogenic variant does not guarantee disease development. Environmental factors, other genetic modifiers, and chance all play roles.
4. **Indels not analyzed.** Only single nucleotide variants (SNPs) are processed -- 23andMe data cannot reliably represent insertions/deletions.
5. **Evolving science.** ClinVar classifications change as research progresses. Update your data files quarterly for the latest information.
6. **WGS coverage matters.** Low-coverage WGS (~1-2x) will miss many positions. Disease analysis (ClinVar position matching) works well, but lifestyle SNP analysis may have fewer hits than 23andMe data.

## Adding Custom SNPs

To add new variants to the lifestyle/health analysis, edit `genetic_health/snp_database.py`:

```python
"rs12345": {
    "gene": "GENE_NAME",
    "category": "Category Name",
    "variants": {
        "AA": {"status": "status_name", "desc": "Description of this genotype", "magnitude": 2},
        "AG": {"status": "carrier", "desc": "Heterozygous carrier", "magnitude": 1},
        "GG": {"status": "reference", "desc": "Typical function", "magnitude": 0},
    },
    "note": "Optional context about this variant"
}
```

## Troubleshooting

**"Genome file not found"** -- Ensure `data/genome.txt` exists or pass the path explicitly: `python -m genetic_health /path/to/genome.txt`

**"ClinVar file not found"** -- Disease risk analysis will be skipped. Run `make update-data` to auto-download, or manually download from ClinVar FTP.

**"PharmGKB files not found"** -- Drug-gene interactions will be skipped. Download from https://www.pharmgkb.org/downloads (free account required).

**"Reference genome not found"** (WGS only) -- Run `make setup` to download and index the GRCh37 reference (~900MB download, ~10GB disk after decompression and indexing).

**"Missing tools: minimap2, samtools..."** -- Enter the Nix shell first: `nix develop`

**fastp not available (Apple Silicon)** -- Install via `brew install fastp`, or the pipeline will skip the QC step automatically.

**Suspicious disease findings** -- The indel filter prevents most false positives. If you see unexpected results for well-known genes (e.g., BRCA), verify the variant is a true SNP and check its ClinVar gold star rating.

## License

This analysis pipeline was developed for personal use. Not for clinical or diagnostic purposes.
