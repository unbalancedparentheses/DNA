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

The pipeline generates a unified report in the `reports/` directory:

1. **GENETIC_HEALTH_REPORT.md** — Comprehensive unified Markdown report (18 sections)
   1. Data Quality (call rate, chromosomes, sex inference)
   2. Executive Summary (high-impact findings, critical disease variants, carrier count)
   3. APOE Haplotype (Alzheimer's risk context)
   4. Blood Type & Traits (ABO + Rh, eye color, hair color, earwax type, freckling)
   5. Mitochondrial Haplogroup (maternal lineage)
   6. Ancestry Estimation (5 superpopulations)
   7. Polygenic Risk Scores (8 conditions)
   8. Pharmacogenomic Star Alleles (CYP2C19, CYP2C9, CYP2D6, DPYD, TPMT, UGT1A1)
   9. ACMG Secondary Findings (81 medically actionable genes)
   10. Gene-Gene Interactions (Epistasis)
   11. Personalized Recommendations (priorities, supplements, diet, lifestyle, monitoring)
   12. Complete Lifestyle/Health Findings by Category (~210+ SNPs)
   13. Pathway Analysis
   14. Disease Risk Analysis (pathogenic, likely pathogenic, carriers, risk factors, drug response, protective)
   15. Carrier Screening (organized by disease system, reproductive context)
   16. Drug-Gene Interactions (PharmGKB + ClinVar combined)
   17. Doctor Card (print-friendly summary)
   18. References & Disclaimer

2. **GENETIC_HEALTH_REPORT.html** — Interactive all-in-one HTML report (23 sections)
   - Everything from the Markdown report, plus:
   - SVG charts and dashboards
   - APOE haplotype display with risk context
   - Trait predictions section
   - ACMG findings with actionability notes
   - Carrier screening organized by system
   - Ancestry donut chart with proportions
   - PRS gauge visualizations for 8 conditions
   - Star allele metabolizer gauges (6 genes)
   - Population frequency badges on annotated findings
   - Search/filter across all findings and tables
   - Sortable table columns (click headers)
   - CSV export of findings to clipboard
   - Collapsible sections
   - Print-optimized doctor card
   - Database links for every rsID

3. **GENETIC_HEALTH_REPORT.pdf** (optional, with `--pdf` flag)

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
│   ├── snp_database.py           # COMPREHENSIVE_SNPS (~210+ curated variants)
│   ├── clinical_context.py       # CLINICAL_CONTEXT + PATHWAYS data
│   ├── ancestry.py               # AIMs database + ancestry estimation
│   ├── prs.py                    # Polygenic risk score models (8 conditions)
│   ├── epistasis.py              # Gene-gene interaction detection (10 models)
│   ├── recommendations.py        # Actionable recommendations from findings
│   ├── blood_type.py             # ABO + Rh blood type prediction from proxy SNPs
│   ├── quality_metrics.py        # Call rate, het rate, chromosome coverage, sex inference
│   ├── mt_haplogroup.py          # Mitochondrial haplogroup estimation (~25 MT SNPs)
│   ├── star_alleles.py           # CPIC-style star allele calling (6 genes)
│   ├── apoe.py                   # APOE epsilon haplotype calling (rs429358 + rs7412)
│   ├── acmg.py                   # ACMG SF v3.2 medically actionable gene filtering
│   ├── carrier_screen.py         # Carrier screening organizer (by disease system)
│   ├── traits.py                 # Visible trait predictions (eye/hair color, earwax, freckling)
│   ├── pipeline.py               # run_full_analysis() orchestrator + CLI
│   ├── wgs_pipeline.py           # WGS FASTQ→reports pipeline
│   ├── update_data.py            # ClinVar auto-download + PharmGKB validation
│   └── reports/
│       ├── __init__.py            # Re-exports generators
│       ├── html_converter.py      # Markdown→HTML converter
│       ├── markdown_reports.py    # Unified Markdown report generator (18 sections)
│       ├── enhanced_html.py       # Interactive HTML report (23 sections)
│       └── pdf_export.py          # PDF export (Chrome headless / weasyprint)
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
│   ├── test_epistasis.py          # Epistasis model + interaction tests
│   ├── test_update_data.py        # Data update + metadata tests
│   ├── test_blood_type.py         # Blood type prediction tests
│   ├── test_quality_metrics.py    # Data quality metrics tests
│   ├── test_mt_haplogroup.py      # MT haplogroup estimation tests
│   ├── test_star_alleles.py       # Star allele calling tests (6 genes)
│   ├── test_recommendations.py    # Recommendation generation tests
│   ├── test_apoe.py               # APOE haplotype calling tests
│   ├── test_acmg.py               # ACMG secondary findings tests
│   ├── test_carrier_screen.py     # Carrier screening tests
│   └── test_traits.py             # Trait prediction tests
├── data/
│   ├── genome.txt                 # 23andMe raw data file
│   ├── clinvar_alleles.tsv        # ClinVar database (~341K variants)
│   ├── clinical_annotations.tsv   # PharmGKB annotations
│   ├── clinical_ann_alleles.tsv   # PharmGKB allele data
│   └── data_versions.json         # Data version metadata (auto-generated)
├── reference/
│   └── human_g1k_v37.fasta.gz    # Reference genome (for WGS pipeline)
└── reports/                       # Generated output (gitignored)
    ├── GENETIC_HEALTH_REPORT.md
    ├── GENETIC_HEALTH_REPORT.html
    └── comprehensive_results.json  # Intermediate JSON for HTML generator
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
mv reports/GENETIC_HEALTH_REPORT.md reports/GENETIC_HEALTH_REPORT_MOM.md
mv reports/GENETIC_HEALTH_REPORT.html reports/GENETIC_HEALTH_REPORT_MOM.html
```

## Module Details

### genetic_health.pipeline
**Main entry point** — orchestrates the entire pipeline.

```python
from genetic_health.pipeline import run_full_analysis
results = run_full_analysis(genome_path, subject_name)
```

### genetic_health.snp_database
Contains curated SNP interpretations for ~210+ health-relevant variants organized by 16 categories:
- Drug Metabolism, Methylation, Neurotransmitters, Caffeine Response
- Cardiovascular, Nutrition, Fitness, Sleep/Circadian
- Detoxification, Autoimmune, Iron Metabolism, Immune Function
- Inflammation, Alcohol Metabolism, Blood Type, Taste Perception

Each SNP entry includes:
- Gene name
- Category
- Genotype interpretations
- Status descriptions
- Impact magnitude (0-6 scale)
- Optional `freq` field with population allele frequencies (EUR, AFR, EAS, SAS, AMR)

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
Calculates Polygenic Risk Scores for 8 conditions using published GWAS effect sizes:

| Condition | Reference |
|-----------|-----------|
| Type 2 Diabetes | Mahajan 2018 |
| Coronary Artery Disease | Nikpay 2015, Aragam 2022 |
| Hypertension | Evangelou 2018 |
| Breast Cancer | Michailidou 2017 |
| Age-Related Macular Degeneration | Fritsche 2016 |
| Prostate Cancer | Schumacher 2018 |
| Ischemic Stroke | Malik 2018 |
| Colorectal Cancer | Huyghe 2019 |

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

### genetic_health.epistasis
Evaluates 10 well-characterized multi-gene interactions where combined effects differ from individual variants.

```python
from genetic_health.epistasis import evaluate_epistasis
results = evaluate_epistasis(genome_by_rsid, lifestyle_results)
# results = [{"name": ..., "risk_level": ..., "mechanism": ..., "recommendations": [...]}, ...]
```

### genetic_health.recommendations
Generates actionable supplement, dietary, and lifestyle recommendations from analysis findings.

```python
from genetic_health.recommendations import generate_recommendations
recs = generate_recommendations(lifestyle_results, clinvar_findings, ancestry, prs_results)
```

### genetic_health.blood_type
Predicts ABO blood group + Rh factor from proxy SNPs (rs505922, rs8176746, rs590787).

```python
from genetic_health.blood_type import predict_blood_type
result = predict_blood_type(genome_by_rsid)
# result = {blood_type: "A+", abo: "A", rh: "+", confidence: "high", details: [...]}
```

### genetic_health.quality_metrics
Computes data quality metrics from the loaded genome and raw file.

```python
from genetic_health.quality_metrics import compute_quality_metrics
metrics = compute_quality_metrics(genome_by_rsid, genome_path)
# metrics = {total_snps, no_call_count, call_rate, chromosomes, has_mt, has_y, het_rate, ...}
```

### genetic_health.mt_haplogroup
Estimates maternal lineage haplogroup from ~25 mitochondrial DNA defining SNPs.

```python
from genetic_health.mt_haplogroup import estimate_mt_haplogroup
result = estimate_mt_haplogroup(genome_by_rsid)
# result = {haplogroup: "H", description: "...", confidence: "high", lineage: "European maternal", ...}
```

### genetic_health.star_alleles
CPIC-style star allele calling for 6 pharmacogenes with metabolizer phenotype classification:
- **CYP2C19** — clopidogrel, PPIs, SSRIs
- **CYP2C9** — warfarin, NSAIDs, phenytoin
- **CYP2D6** — codeine, tamoxifen, many antidepressants (caveat: copy number variants undetectable)
- **DPYD** — fluoropyrimidine chemotherapy (5-FU, capecitabine) — deficiency can be fatal
- **TPMT** — thiopurines (azathioprine, 6-MP) — deficiency causes myelosuppression
- **UGT1A1** — irinotecan, atazanavir — *28 causes Gilbert syndrome

```python
from genetic_health.star_alleles import call_star_alleles
results = call_star_alleles(genome_by_rsid)
# results["CYP2C19"] = {gene, diplotype: "*1/*2", phenotype: "Intermediate Metabolizer", ...}
```

### genetic_health.apoe
Combines rs429358 + rs7412 into APOE epsilon haplotypes (e2/e2 through e4/e4) with Alzheimer's risk context.

```python
from genetic_health.apoe import call_apoe_haplotype
result = call_apoe_haplotype(genome_by_rsid)
# result = {apoe_type: "e3/e4", risk_level: "elevated", alzheimer_or: 2.8, confidence: "high", ...}
```

### genetic_health.acmg
Filters ClinVar pathogenic/likely_pathogenic findings against the ACMG SF v3.2 list of 81 medically actionable genes (BRCA1/2, Lynch syndrome genes, cardiac genes, etc.).

```python
from genetic_health.acmg import flag_acmg_findings
result = flag_acmg_findings(disease_findings)
# result = {acmg_findings: [...], genes_screened: 81, genes_with_variants: N, summary: "..."}
```

### genetic_health.carrier_screen
Reorganizes heterozygous recessive ClinVar findings into a structured carrier screening report grouped by disease system with reproductive context.

```python
from genetic_health.carrier_screen import organize_carrier_findings
result = organize_carrier_findings(disease_findings)
# result = {carriers: [...], total_carriers: N, by_system: {...}, couples_relevant: [...]}
```

### genetic_health.traits
Predicts visible traits from well-established SNP associations:
- **Eye color** — rs12913832 (HERC2) + rs1800407 (OCA2)
- **Hair color** — rs1805007, rs1805008 (MC1R)
- **Earwax type** — rs17822931 (ABCC11)
- **Freckling** — MC1R variants

```python
from genetic_health.traits import predict_traits
result = predict_traits(genome_by_rsid)
# result = {eye_color: {prediction, confidence, ...}, hair_color: {...}, ...}
```

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
