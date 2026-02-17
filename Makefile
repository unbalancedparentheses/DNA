SHELL := /bin/bash
NIX := nix --extra-experimental-features 'nix-command flakes'

# Default FASTQ input - override with: make pipeline FASTQ=/path/to/reads.fastq
FASTQ ?= ../22374_1.fastq
NAME ?=
THREADS ?= 12

# Directories
REF_DIR := reference
WORK_DIR := wgs_work
DATA_DIR := data
REPORTS_DIR := reports

# Reference files
REF_FASTA := $(REF_DIR)/human_g1k_v37.fasta
REF_MMI := $(REF_DIR)/human_g1k_v37.mmi
REF_FAI := $(REF_FASTA).fai

# Outputs
GENOME := $(DATA_DIR)/genome.txt
RSID_LOOKUP := $(DATA_DIR)/rsid_positions_grch37.json

# Build name flag if provided
ifdef NAME
  NAME_FLAG := --name "$(NAME)"
else
  NAME_FLAG :=
endif

.PHONY: all setup pipeline analysis test clean clean-all help shell fastp update-data validate-pharmgkb data-status

## Default target: show help
help:
	@echo "Genetic Health Analysis Pipeline"
	@echo "================================"
	@echo ""
	@echo "Setup (one-time):"
	@echo "  make setup          Download reference genome + create indexes"
	@echo "  make fastp          Install fastp via Homebrew (for QC step)"
	@echo ""
	@echo "Run:"
	@echo "  make pipeline       Full pipeline: FASTQ -> reports"
	@echo "  make analysis       Analysis only (requires existing genome.txt)"
	@echo ""
	@echo "Options:"
	@echo "  FASTQ=path          Input FASTQ file (default: ../22374_1.fastq)"
	@echo "  NAME=\"John Doe\"     Subject name for reports"
	@echo "  THREADS=12          CPU threads (default: 12)"
	@echo ""
	@echo "Examples:"
	@echo "  make setup"
	@echo "  make pipeline"
	@echo "  make pipeline FASTQ=../22374_1.fastq.gz NAME=\"Subject\""
	@echo "  make analysis NAME=\"Subject\""
	@echo ""
	@echo "Data:"
	@echo "  make update-data    Download latest ClinVar + print PharmGKB instructions"
	@echo "  make validate-pharmgkb  Validate PharmGKB after manual download"
	@echo "  make data-status    Show current data versions"
	@echo ""
	@echo "Test:"
	@echo "  make test           Run test suite"
	@echo ""
	@echo "Cleanup:"
	@echo "  make clean          Remove intermediate files (keep BAM/VCF)"
	@echo "  make clean-all      Remove all generated files"
	@echo ""
	@echo "Other:"
	@echo "  make shell          Enter Nix development shell"

## Enter Nix development shell
shell:
	$(NIX) develop

## Install fastp via Homebrew (not in nixpkgs for Apple Silicon)
fastp:
	brew install fastp

## One-time setup: download reference genome and build indexes
setup: $(REF_MMI) $(REF_FAI)

$(REF_DIR):
	mkdir -p $(REF_DIR)

$(REF_FASTA): | $(REF_DIR)
	@echo "Downloading GRCh37 reference genome (~900MB)..."
	$(NIX) develop --command bash -c '\
		curl -L --progress-bar \
			-o $(REF_FASTA).gz \
			"https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/human_g1k_v37.fasta.gz" && \
		echo "Decompressing..." && \
		gunzip $(REF_FASTA).gz'

$(REF_FAI): $(REF_FASTA)
	@echo "Creating samtools faidx index..."
	$(NIX) develop --command samtools faidx $(REF_FASTA)

$(REF_MMI): $(REF_FASTA)
	@echo "Creating minimap2 index (~5 min)..."
	$(NIX) develop --command minimap2 -x sr -d $(REF_MMI) $(REF_FASTA)

## Full pipeline: FASTQ -> QC -> align -> call -> convert -> analyze
pipeline: $(REF_MMI)
	@test -f "$(FASTQ)" || { echo "FASTQ not found: $(FASTQ)"; echo "Usage: make pipeline FASTQ=/path/to/reads.fastq"; exit 1; }
	$(NIX) develop --command python -m genetic_health.wgs_pipeline \
		"$(FASTQ)" \
		--threads $(THREADS) \
		$(NAME_FLAG)

## Run analysis only on existing genome.txt
analysis:
	@test -f "$(GENOME)" || { echo "genome.txt not found. Run 'make pipeline' first."; exit 1; }
	$(NIX) develop --command python -m genetic_health \
		$(GENOME) \
		$(NAME_FLAG)

## Run test suite
test:
	$(NIX) develop --command python -m pytest tests/ -v

## Download latest ClinVar + print PharmGKB instructions
update-data:
	$(NIX) develop --command python -m genetic_health.update_data clinvar
	@echo ""
	@echo "PharmGKB requires manual download from https://www.pharmgkb.org/downloads"
	@echo "Run 'make validate-pharmgkb' after downloading."

## Validate PharmGKB after manual download
validate-pharmgkb:
	$(NIX) develop --command python -m genetic_health.update_data pharmgkb

## Show current data versions
data-status:
	$(NIX) develop --command python -m genetic_health.update_data --status

## Remove intermediate files (keep BAM, VCF, and reports)
clean:
	rm -f $(WORK_DIR)/trimmed.fastq.gz
	rm -f $(WORK_DIR)/fastp.*
	@echo "Cleaned intermediate files"

## Remove all generated files (keep reference and data)
clean-all: clean
	rm -rf $(WORK_DIR)
	rm -f $(DATA_DIR)/genome.txt
	rm -f $(DATA_DIR)/genome.txt.bak
	rm -f $(RSID_LOOKUP)
	rm -rf $(REPORTS_DIR)/*.md $(REPORTS_DIR)/*.json
	@echo "Cleaned all generated files"
