# Genetic Health Analysis Pipeline

Processes 23andMe raw data or whole-genome sequencing FASTQ files into detailed health reports covering lifestyle genetics, disease risk, drug interactions, and actionable protocols.

## Requirements

- [Nix](https://nixos.org/download/) (provides Python, minimap2, samtools, bcftools)
- Optional: [fastp](https://github.com/OpenGene/fastp) for QC (`brew install fastp`)

## Quick Start

```bash
# Enter development shell
nix develop

# Run with 23andMe data
cp ~/Downloads/genome.txt data/
python scripts/run_full_analysis.py
python scripts/run_full_analysis.py --name "John Doe"

# Or run full WGS pipeline (FASTQ -> reports)
make setup                          # one-time: download reference genome
make pipeline FASTQ=/path/to/reads.fastq NAME="John Doe"
```

## Output

Three reports in `reports/`:

| Report | Contents |
|--------|----------|
| `EXHAUSTIVE_GENETIC_REPORT.md` | Drug metabolism, methylation, nutrition, fitness, cardiovascular, sleep, PharmGKB interactions |
| `EXHAUSTIVE_DISEASE_RISK_REPORT.md` | Pathogenic variants, carrier status, risk factors, drug response, protective variants |
| `ACTIONABLE_HEALTH_PROTOCOL_V3.md` | Supplement/diet/exercise protocols, monitoring schedule, drug interaction card |

## Make Targets

```
make help       Show all targets
make setup      Download reference genome + create indexes (~10GB disk)
make pipeline   Full WGS pipeline: FASTQ -> QC -> align -> call -> convert -> analyze
make analysis   Run analysis on existing data/genome.txt
make test       Run test suite
make clean      Remove intermediate files
make clean-all  Remove all generated files
```

## Data Files

Place in `data/`:

| File | Source | Required |
|------|--------|----------|
| `genome.txt` | [23andMe raw data download](https://you.23andme.com/) or WGS pipeline output | Yes |
| `clinvar_alleles.tsv` | [ClinVar FTP](https://ftp.ncbi.nlm.nih.gov/pub/clinvar/) | For disease risk |
| `clinical_annotations.tsv` | [PharmGKB](https://www.pharmgkb.org/downloads) | For drug interactions |
| `clinical_ann_alleles.tsv` | [PharmGKB](https://www.pharmgkb.org/downloads) | For drug interactions |

## Tests

```bash
python -m pytest tests/ -v    # 55 tests
```

Covers genome parsing, SNP database validation, VCF conversion, disease analysis logic, and project structure verification. All tests use synthetic data.

## Limitations

- Not a clinical diagnosis — informational purposes only
- Only single nucleotide variants (indels filtered out)
- Associations may vary by ancestry
- Variant classifications evolve with research

## Project Structure

```
DNA/
├── scripts/
│   ├── run_full_analysis.py           # Main entry point
│   ├── comprehensive_snp_database.py  # ~200 curated SNP interpretations
│   ├── generate_exhaustive_report.py  # Report generator with clinical context
│   ├── wgs_pipeline.py               # FASTQ -> 23andMe format
│   └── setup_reference.sh            # Reference genome setup
├── tests/                             # 5 test files, 55 tests
├── data/                              # Input data (gitignored)
├── reference/                         # Reference genome (gitignored)
└── reports/                           # Generated output (gitignored)
```
