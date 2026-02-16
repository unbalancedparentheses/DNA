# Genetic Health Analysis Workspace

This workspace contains a comprehensive genetic health analysis pipeline that processes 23andMe raw data to generate detailed health reports.

## Quick Start

```bash
# Run full analysis with default genome file (data/genome.txt)
cd /path/to/DNA
python scripts/run_full_analysis.py

# Run with a custom genome file
python scripts/run_full_analysis.py /path/to/genome.txt

# Run with subject name included in reports
python scripts/run_full_analysis.py --name "John Doe"
python scripts/run_full_analysis.py /path/to/genome.txt --name "Jane Doe"

# Run tests
python -m pytest tests/ -v
```

## Output

The pipeline generates three reports in the `reports/` directory:

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
   - Daily protocol (morning/midday/evening stacks)
   - Dietary framework with specific targets
   - Exercise protocol with weekly structure
   - Blood pressure management
   - Risk factor monitoring by condition
   - Comprehensive drug interactions (PharmGKB + ClinVar)
   - Testing & monitoring schedule
   - 90-day implementation checklist
   - Quick reference cards

## Directory Structure

```
DNA/
├── CLAUDE.md
├── .gitignore
├── flake.nix
├── flake.lock
├── Makefile
├── data/
│   ├── genome.txt                 # 23andMe raw data file
│   ├── clinvar_alleles.tsv        # ClinVar database (~341K variants)
│   ├── clinical_annotations.tsv   # PharmGKB annotations
│   └── clinical_ann_alleles.tsv   # PharmGKB allele data
├── reference/
│   └── human_g1k_v37.fasta.gz    # Reference genome (for WGS pipeline)
├── scripts/
│   ├── run_full_analysis.py       # MAIN ENTRY POINT - runs everything
│   ├── comprehensive_snp_database.py  # Curated SNP database
│   ├── generate_exhaustive_report.py  # Report generator
│   ├── wgs_pipeline.py           # Whole-genome sequencing pipeline
│   └── setup_reference.sh        # Reference genome setup
├── tests/
│   ├── test_genome_loading.py
│   ├── test_snp_analysis.py
│   ├── test_vcf_conversion.py
│   ├── test_disease_analysis.py
│   └── test_path_resolution.py
└── reports/                       # Generated output (gitignored)
    ├── EXHAUSTIVE_GENETIC_REPORT.md
    ├── EXHAUSTIVE_DISEASE_RISK_REPORT.md
    └── ACTIONABLE_HEALTH_PROTOCOL_V3.md
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
python scripts/run_full_analysis.py data/genome_mom.txt --name "Mom"

# Rename outputs to preserve them
mv reports/EXHAUSTIVE_GENETIC_REPORT.md reports/EXHAUSTIVE_GENETIC_REPORT_MOM.md
mv reports/EXHAUSTIVE_DISEASE_RISK_REPORT.md reports/EXHAUSTIVE_DISEASE_RISK_REPORT_MOM.md
mv reports/ACTIONABLE_HEALTH_PROTOCOL.md reports/ACTIONABLE_HEALTH_PROTOCOL_MOM.md
```

## Script Details

### run_full_analysis.py
**Main entry point** - orchestrates the entire pipeline.

```python
from scripts.run_full_analysis import run_full_analysis
results = run_full_analysis(genome_path, subject_name)
```

### comprehensive_snp_database.py
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

### run_full_analysis.py — Disease Risk Analysis
The disease risk analysis (integrated into `run_full_analysis.py`) scans the genome against ClinVar. Important implementation detail:

**Only true SNPs are analyzed** — indels (insertions/deletions) are filtered out because 23andMe data cannot reliably represent them. This prevents false positives.

```python
# Critical filter to avoid false positives
if len(ref_allele) != 1 or len(alt_allele) != 1:
    continue  # Skip indels
```

### generate_exhaustive_report.py
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
# Download latest ClinVar
wget https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/variant_summary.txt.gz
# Process and convert to required format (see clinvar processing scripts)
```

### PharmGKB (recommended: quarterly)
Download from https://www.pharmgkb.org/downloads after creating free account.

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

To add new SNPs to the lifestyle/health analysis, edit `scripts/comprehensive_snp_database.py`:

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
