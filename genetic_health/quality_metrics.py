"""Compute data quality metrics from 23andMe raw genome data."""


def compute_quality_metrics(genome_by_rsid, genome_path=None):
    """Compute quality metrics from loaded genome data and optionally the raw file.

    Parameters
    ----------
    genome_by_rsid : dict
        Loaded genome dict {rsid: {chromosome, position, genotype}}.
    genome_path : Path or None
        Path to the raw genome file for no-call counting.

    Returns
    -------
    dict with keys: total_snps, no_call_count, call_rate, chromosomes,
        has_mt, has_y, mt_snp_count, autosomal_count, het_rate.
    """
    total_snps = len(genome_by_rsid)

    # Count SNPs per chromosome
    chromosomes = {}
    autosomal_het = 0
    autosomal_total = 0

    for entry in genome_by_rsid.values():
        raw = entry.get("chromosome")
        if raw is None or raw == "":
            continue
        chrom = str(raw).upper().replace("CHR", "")
        if not chrom:
            continue
        chromosomes[chrom] = chromosomes.get(chrom, 0) + 1

        # Heterozygosity rate on autosomes (chr 1-22).
        # Haploid SNPs (MT, X in males) are excluded since they can't be het.
        if chrom.isdigit() and 1 <= int(chrom) <= 22:
            autosomal_total += 1
            gt = entry.get("genotype", "")
            if len(gt) == 2 and gt[0] != gt[1]:
                autosomal_het += 1

    has_mt = "MT" in chromosomes
    has_y = "Y" in chromosomes
    mt_snp_count = chromosomes.get("MT", 0)
    autosomal_count = sum(
        cnt for ch, cnt in chromosomes.items()
        if ch.isdigit() and 1 <= int(ch) <= 22
    )
    het_rate = autosomal_het / autosomal_total if autosomal_total > 0 else 0.0

    # Count no-calls from the raw file
    no_call_count = 0
    if genome_path is not None:
        try:
            with open(genome_path, "r") as fh:
                for line in fh:
                    line = line.strip()
                    if not line or line.startswith("#"):
                        continue
                    parts = line.split("\t")
                    if len(parts) >= 4 and parts[3] == "--":
                        no_call_count += 1
        except (OSError, IOError):
            pass

    call_rate = (
        total_snps / (total_snps + no_call_count)
        if (total_snps + no_call_count) > 0
        else 0.0
    )

    return {
        "total_snps": total_snps,
        "no_call_count": no_call_count,
        "call_rate": call_rate,
        "chromosomes": chromosomes,
        "has_mt": has_mt,
        "has_y": has_y,
        "mt_snp_count": mt_snp_count,
        "autosomal_count": autosomal_count,
        "het_rate": het_rate,
    }
