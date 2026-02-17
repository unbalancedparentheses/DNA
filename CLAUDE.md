# Genetic Health Analysis Workspace

This workspace contains a comprehensive genetic health analysis pipeline that processes 23andMe raw data to generate detailed health reports.

## Quick Start

```bash
# Run full analysis with default genome file (data/genome.txt)
cd /path/to/DNA
python -m genetic_health

# Run with a custom genome file
python -m genetic_health /path/to/genome.txt

# Run with subject name included in reports
python -m genetic_health --name "John Doe"
python -m genetic_health /path/to/genome.txt --name "Jane Doe"

# Run tests
python -m pytest tests/ -v
```

## Output

The pipeline generates four reports in the `reports/` directory:

1. **EXHAUSTIVE_GENETIC_REPORT.md** - Lifestyle/health genetics analysis
   - Drug metabolism (CYP enzymes, warfarin sensitivity)
   - Methylation (MTHFR, COMT, MTRR)
   - Nutrition (vitamin D, omega-3, lactose)
   - Fitness (muscle fiber type, exercise response)
   - Cardiovascular (blood pressure genes, clotting)
   - Sleep/circadian rhythm
   - PharmGKB drug-gene interactions

2. **EXHAUSTIVE_DISEASE_RISK_REPORT.md** - Clinical variant analysis
   - Pathogenic variants (affected status)
   - Carrier status for recessive conditions
   - Risk factors
   - Drug response variants
   - Protective variants

3. **ACTIONABLE_HEALTH_PROTOCOL_V3.md** - Comprehensive protocol
   - Critical disease findings summary (pathogenic, carrier status)
   - Ancestry estimation and population context
   - Polygenic risk scores for 5 conditions
   - Supplement and dietary recommendations
   - Lifestyle and exercise guidance
   - Blood pressure management
   - Risk factor monitoring by condition
   - Comprehensive drug interactions (PharmGKB + ClinVar)
   - Testing & monitoring schedule

4. **ENHANCED_HEALTH_REPORT.html** - All-in-one interactive HTML report
   - SVG charts and dashboards
   - Ancestry donut chart with proportions
   - PRS gauge visualizations for 5 conditions
   - Search/filter across all findings and tables
   - Sortable table columns (click headers)
   - CSV export of findings to clipboard
   - Collapsible sections
   - Print-optimized doctor card
   - Database links for every rsID

## Directory Structure

```
DNA/
├── CLAUDE.md
├── .gitignore
├── flake.nix
├── flake.lock
├── Makefile
├── pyproject.toml
├── genetic_health/               # Python package
│   ├── __init__.py               # Public API
│   ├── __main__.py               # python -m genetic_health
│   ├── config.py                 # Centralized paths (BASE_DIR, DATA_DIR, etc.)
│   ├── loading.py                # load_genome(), load_pharmgkb()
│   ├── analysis.py               # analyze_lifestyle_health(), load_clinvar_and_analyze()
│   ├── ancestry.py               # AIMs database + ancestry estimation
│   ├── prs.py                    # Polygenic risk score models + scoring
│   ├── update_data.py            # ClinVar auto-download + PharmGKB validation
│   ├── snp_database.py           # COMPREHENSIVE_SNPS (~200 curated variants)
│   ├── clinical_context.py       # CLINICAL_CONTEXT + PATHWAYS data
│   ├── pipeline.py               # run_full_analysis() orchestrator + CLI
│   ├── wgs_pipeline.py           # WGS FASTQ→reports pipeline
│   └── reports/
│       ├── __init__.py            # Re-exports generators
│       ├── html_converter.py      # Markdown→HTML converter
│       ├── section_builders.py    # Section generators for exhaustive report
│       ├── markdown_reports.py    # 3 Markdown report generators
│       └── enhanced_html.py       # All-in-one interactive HTML report
├── scripts/
│   └── setup_reference.sh        # Reference genome setup (shell script)
├── tests/
│   ├── conftest.py
│   ├── test_genome_loading.py
│   ├── test_snp_analysis.py
│   ├── test_vcf_conversion.py
│   ├── test_disease_analysis.py
│   ├── test_path_resolution.py
│   ├── test_html_reports.py
│   ├── test_ancestry.py           # Ancestry marker + scoring tests
│   ├── test_prs.py                # PRS model + calculation tests
│   └── test_update_data.py        # Data update + metadata tests
├── data/
│   ├── genome.txt                 # 23andMe raw data file
│   ├── clinvar_alleles.tsv        # ClinVar database (~341K variants)
│   ├── clinical_annotations.tsv   # PharmGKB annotations
│   ├── clinical_ann_alleles.tsv   # PharmGKB allele data
│   └── data_versions.json         # Data version metadata (auto-generated)
├── reference/
│   └── human_g1k_v37.fasta.gz    # Reference genome (for WGS pipeline)
└── reports/                       # Generated output (gitignored)
    ├── EXHAUSTIVE_GENETIC_REPORT.md
    ├── EXHAUSTIVE_DISEASE_RISK_REPORT.md
    ├── ACTIONABLE_HEALTH_PROTOCOL_V3.md
    └── ENHANCED_HEALTH_REPORT.html
```

## Data Requirements

### Genome File Format
The pipeline expects 23andMe raw data format (tab-separated):
```
# rsid  chromosome  position  genotype
rs123   1           12345     AG
```

- Lines starting with `#` are ignored
- Genotype `--` indicates no call (ignored)
- Both rsIDs and position-based matching are used

### Required Data Files
The analysis requires these files in `data/`:

1. **genome.txt** - Your 23andMe raw data download
2. **clinvar_alleles.tsv** - ClinVar variant database
   - Download from: https://ftp.ncbi.nlm.nih.gov/pub/clinvar/
3. **clinical_annotations.tsv** - PharmGKB annotations
4. **clinical_ann_alleles.tsv** - PharmGKB allele data
   - Download from: https://www.pharmgkb.org/downloads

## Running Analysis for Family Members

To run analysis for different people:

```bash
# Copy their genome file to data/ with a unique name
cp ~/Downloads/genome_mom.txt data/genome_mom.txt

# Run analysis
python -m genetic_health data/genome_mom.txt --name "Mom"

# Rename outputs to preserve them
mv reports/EXHAUSTIVE_GENETIC_REPORT.md reports/EXHAUSTIVE_GENETIC_REPORT_MOM.md
mv reports/EXHAUSTIVE_DISEASE_RISK_REPORT.md reports/EXHAUSTIVE_DISEASE_RISK_REPORT_MOM.md
mv reports/ACTIONABLE_HEALTH_PROTOCOL_V3.md reports/ACTIONABLE_HEALTH_PROTOCOL_MOM.md
```

## Module Details

### genetic_health.pipeline
**Main entry point** — orchestrates the entire pipeline.

```python
from genetic_health.pipeline import run_full_analysis
results = run_full_analysis(genome_path, subject_name)
```

### genetic_health.snp_database
Contains curated SNP interpretations for ~200 health-relevant variants organized by category:
- Drug Metabolism
- Methylation
- Neurotransmitters
- Caffeine Response
- Cardiovascular
- Nutrition
- Fitness
- Sleep/Circadian
- And more...

Each SNP entry includes:
- Gene name
- Category
- Genotype interpretations
- Status descriptions
- Impact magnitude (0-6 scale)

### genetic_health.analysis — Disease Risk Analysis
The disease risk analysis scans the genome against ClinVar. Important implementation detail:

**Only true SNPs are analyzed** — indels (insertions/deletions) are filtered out because 23andMe data cannot reliably represent them. This prevents false positives.

```python
# Critical filter to avoid false positives
if len(ref_allele) != 1 or len(alt_allele) != 1:
    continue  # Skip indels
```

### genetic_health.ancestry
Estimates ancestry proportions from ~55 Ancestry-Informative Markers (AIMs).

```python
from genetic_health.ancestry import estimate_ancestry
result = estimate_ancestry(genome_by_rsid)
# result = {proportions: {EUR: 0.82, ...}, markers_found, confidence, top_ancestry, details}
```

5 superpopulations (1000 Genomes): EUR, AFR, EAS, SAS, AMR.
Uses maximum-likelihood scoring with softmax normalization.
Includes population-specific warnings via `get_population_warnings(gene, status)`.

### genetic_health.prs
Calculates Polygenic Risk Scores for 5 conditions using published GWAS effect sizes:

| Condition | Reference |
|-----------|-----------|
| Type 2 Diabetes | Mahajan 2018 |
| Coronary Artery Disease | Nikpay 2015, Aragam 2022 |
| Hypertension | Evangelou 2018 |
| Breast Cancer | Michailidou 2017 |
| Age-Related Macular Degeneration | Fritsche 2016 |

```python
from genetic_health.prs import calculate_prs
results = calculate_prs(genome_by_rsid, ancestry_proportions)
# results[condition_id] = {name, percentile, risk_category, ...}
```

Risk categories: low (<20th), average (20-80), elevated (80-95), high (>95th percentile).
Flags non-European ancestry as potentially less applicable.

### genetic_health.update_data
Automated data updates with version tracking.

```bash
python -m genetic_health.update_data clinvar      # Download + process ClinVar
python -m genetic_health.update_data pharmgkb     # Validate PharmGKB files
python -m genetic_health.update_data --status     # Show data versions

# Or via Makefile:
make update-data          # Download ClinVar + print PharmGKB instructions
make validate-pharmgkb    # Validate PharmGKB after manual download
make data-status          # Show current data versions
```

### genetic_health.clinical_context
Contains clinical context database with detailed interpretations, mechanisms, and recommendations for each gene/status combination.

## Interpretation Guide

### Impact Magnitude Scale (0-6)
- **0**: Informational only
- **1**: Low impact - minor effect
- **2**: Moderate impact - worth noting
- **3**: High impact - actionable
- **4-6**: Very high impact - requires attention

### ClinVar Confidence (Gold Stars)
- **4 stars**: Practice guideline / Expert panel
- **3 stars**: Multiple submitters, no conflicts
- **2 stars**: Multiple submitters with conflicts, or single with criteria
- **1 star**: Single submitter
- **0 stars**: No assertion criteria

### Zygosity Interpretation
- **Homozygous**: Both copies of variant - higher effect
- **Heterozygous + Recessive**: Carrier only - reproductive implications
- **Heterozygous + Dominant**: One copy sufficient - may be affected

## Limitations

1. **Not a clinical diagnosis** - For informational purposes only
2. **Population differences** - Associations may vary by ancestry
3. **Incomplete penetrance** - Not everyone with variant develops condition
4. **Evolving science** - Classifications change as research progresses
5. **Indels not analyzed** - Only single nucleotide variants from 23andMe

## Updating Data

### ClinVar (recommended: quarterly)
```bash
# Automated download + processing (preferred)
make update-data
# Or directly:
python -m genetic_health.update_data clinvar
```

This downloads `variant_summary.txt.gz` from NCBI FTP, filters to GRCh37, and creates `clinvar_alleles.tsv`. Version metadata is saved to `data/data_versions.json`.

### PharmGKB (recommended: quarterly)
PharmGKB requires a free account for downloads:
1. Download from https://www.pharmgkb.org/downloads
2. Place `clinical_annotations.tsv` and `clinical_ann_alleles.tsv` in `data/`
3. Validate: `make validate-pharmgkb`

### Check Data Status
```bash
make data-status
# Or: python -m genetic_health.update_data --status
```

## Troubleshooting

### "Genome file not found"
Ensure genome.txt exists in `data/` or provide full path as argument.

### "ClinVar file not found"
Disease risk analysis will be skipped. Download clinvar_alleles.tsv if needed.

### "PharmGKB files not found"
Drug-gene interactions will be skipped. Download from PharmGKB if needed.

### False positives in disease report
The indel filter should prevent this. If you see suspicious findings (especially for BRCA or other well-known genes), verify the variant is a true SNP.

## Adding Custom SNPs

To add new SNPs to the lifestyle/health analysis, edit `genetic_health/snp_database.py`:

```python
"rs12345": {
    "gene": "GENE_NAME",
    "category": "Category Name",
    "variants": {
        "AA": {"status": "status_name", "desc": "Description", "magnitude": 2},
        "AG": {"status": "other_status", "desc": "Description", "magnitude": 1},
        "GG": {"status": "reference", "desc": "Description", "magnitude": 0},
    },
    "note": "Optional additional context"
}
```

## Contact

This analysis pipeline was developed for personal use. Not for clinical or diagnostic purposes.
