from os.path import join


UNPHASED_SNPS_DIR = "output/snps/"
UNPHASED_SNPS = join(UNPHASED_SNPS_DIR, "chr{chromosome}.vcf.gz")
PHASED_SNPS_DIR = "output/phased_snps/"
PHASED_SNPS = join(PHASED_SNPS_DIR, "phased.vcf.gz")

normal_baf = join("output", "baf", "{phasing}-snps", "normal.1bed")
tumor_baf = join("output", "baf", "{phasing}-snps", "tumor.1bed")
count_alleles_snps_dir = join("output", "count_alleles", "{phasing}-snps", "temp")

count_reads_dir = join("output", "count_reads", "{phasing}-snps")
total_reads_tsv = join(count_reads_dir, "total.tsv")

combined_counts = join("output", "bb", "{phasing}-snps", "bulk.bb")

clustered_dir = join("output", "bbc", "{phasing}-snps")
clustered_genomic_bins = join(clustered_dir, "bulk.bbc")
clustered_genomic_segments = join(clustered_dir, "bulk.seg")

results_dir = join("output", "results", "{phasing}-snps")
CN_bins = join(results_dir, "best.bbc.ucn")
CN_segs = join(results_dir, "best.seg.ucn")

include: "rules/make-plots.smk"

rule all:
    input:
        expand(BIN_2D_PLOT, phasing=["phased", "unphased"])
        expand(CLUSTERED_BIN_2D_PLOT, phasing=["phased", "unphased"])
        expand(BAF_PLOT, phasing=["phased", "unphased"])
        expand(CLUSTERED_BAF_PLOT, phasing=["phased", "unphased"])
        expand(RDR_PLOT, phasing=["phased", "unphased"])
        expand(CLUSTERED_RDR_PLOT, phasing=["phased", "unphased"])



rule hatchet_compute_cn:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        clustered_genomic_bins,
        clustered_genomic_segments
    params:
        in_dir=clustered_dir,
        out_dir=results_dir,
    output:
        CN_bins,
        CN_segs,
    shell:
        "module use --prepend /data/CDSL_Gurobi_users/modules && "
        "module load gurobi && "
        "export GRB_LICENSE_FILE=/data/CDSL_Gurobi_users/gurobi910/gurobi.lic && "
        "hatchet compute-cn "
        "-i {params.in_dir} "
        "-x {params.out_dir} "
        "--clones 2,8 " # default

rule hatchet_cluster_bins:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        combined_counts
    output:
        clustered_genomic_bins,
        clustered_genomic_segments
    shell:
        "hatchet cluster-bins "
        "{input} "
        "--outsegments {output.clustered_genomic_segments} "
        "--outbins {output.clustered_genomic_bins} "
        "--diploidbaf 0.08 "
        "--seed 0 " # to ensure reproducibility


# combine tumor bin counts, normal bin counts and tumor allele counts to obtain the read-depth ratio
# and the mean B-allele frequency of each bin
rule hatchet_combine_counts:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        total_tsv = total_reads_tsv,
        tumor_baf = tumor_baf
    output:
        combined_counts # tab separated file that include many fields see - http://compbio.cs.brown.edu/hatchet/doc_combine_counts.html
    shell:
        "hatchet combine-counts "
        "--refversion hg19 "
        "--totalcounts {input.total_tsv} "
        "--array output/count_reads "
        "--baffile {input.tumor_baf} "
        "--msr 3000 " # minimum number of SNP-covering reads
        "--mtr 5000 " # minimum number of total reads per bin
        "--outfile {output} "

# count the mapped sequencing reads in bins of fixed and given length, uniformly for a BAM file of a normal sample
# and one or more BAM files of tumor samples
rule hatchet_count_reads:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        normal = "data/normal.bam",
        tumor = ["data/bulk_03clone1_06clone0_01normal.sorted.bam", "data/bulk_08clone1_Noneclone0_02normal.sorted.bam", "data/bulk_Noneclone1_09clone0_01normal.sorted.bam"],
        tumor_baf = tumor_baf
    output:
        out_dir = directory(count_reads_dir),
        total_tsv = total_reads_tsv
    shell:
        "hatchet count-reads "
        "--tumor {input.tumor} "
        "--normal {input.normal} "
        "--samples normal tumor1 tumor2 tumor3 "
        "--refversion hg19 "
        "--baffile {input.tumor_baf} "
        "--chromosomes chr22 "
        "--outdir {output.out_dir} "


def get_snps(wildcards):
    if wildcards["phasing"] == "phased":
        return PHASED_SNPS
    else:
        return expand(PHASED_SNPS, chromosome=22)

rule hatchet_count_alleles:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        normal="data/normal.bam",
        tumor=["data/bulk_03clone1_06clone0_01normal.sorted.bam", "data/bulk_08clone1_Noneclone0_02normal.sorted.bam", "data/bulk_Noneclone1_09clone0_01normal.sorted.bam"],
        snps = get_snps,
    output:
        normal_baf = normal_baf, # contains the number of reads for the major and minor allele
        tumor_baf = tumor_baf,
        snps_dir = directory(count_alleles_snps_dir)
    shell:
        "mkdir {output.snps_dir} && "
        "hatchet count-alleles "
        "--tumors {input.tumor} "
        "--normal {input.normal} "
        "--samples normal tumor1 tumor2 tumor3 "
        "--reference data/hg19.fa "
        "--snps {input.snps} "
        "--chromosomes chr22 "
        "--outputnormal {output.normal_baf} "
        "--outputtumors {output.tumor_baf} "
        "--outputsnps {output.snps_dir} "
        # "--regions " required for WXS (although I don't think this is true, especially since SNPs are passed)
        "--mincov 8 "
        "--maxcov 300 "

rule hatchet_phasing:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        snps=expand(UNPHASED_SNPS, chromosome=22)
    output:
        PHASED_SNPS
    shell:
        "hatchet phase-snps "
        "--snps {input.snps} "
        "--refpaneldir data/reference/panel/1000GP_Phase3 "
        "--refgenome data/hg19.fa "
        "--chrnotation "
        "--outdir {PHASED_SNPS_DIR} "


# used by phase-snps
rule run_hatchet_download_panel:
    # conda:
    #     "envs/HATCHet-env.yaml"
    output:
        directory("data/reference/panel/1000GP_Phase3")
    shell:
        # "hatchet run hatchet-download-panel.ini"
        "hatchet download-panel --refpaneldir data/reference/panel --refpanel 1000GP_Phase3"


rule run_hatchet_genotype_snps:
  # conda:
  #     "envs/HATCHet-env.yaml"
  input:
      bam="data/normal.bam",
      bai="data/normal.bam.bai",
      dict="data/hg19.dict",
      ref="data/hg19.fa"
  output:
      # directory("output/snps/"),
      expand(UNPHASED_SNPS, chromosome=22) #list(range(1, 23)) + ["X", "Y"])
  shell:
      "hatchet genotype-snps "
      "--normal {input.bam} "
      "--reference {input.ref} "
      "--mincov 8 "
      "--maxcov 300 "
      "--chromosomes chr22 "
      "--outputsnps {UNPHASED_SNPS_DIR}"


# rule run_hatchet_init:
#     conda:
#         "envs/HATCHet-env.yaml"
#     input:
#         "data/normal.bam",
#         "data/normal.bam.bai",
#         "data/bulk_03clone1_06clone0_01normal.sorted.bam",
#         "data/bulk_03clone1_06clone0_01normal.sorted.bam.bai",
#         "data/bulk_08clone1_Noneclone0_02normal.sorted.bam",
#         "data/bulk_08clone1_Noneclone0_02normal.sorted.bam.bai",
#         "data/bulk_Noneclone1_09clone0_01normal.sorted.bam",
#         "data/bulk_Noneclone1_09clone0_01normal.sorted.bam.bai",
#         "data/hg19.dict",
#         "data/hg19.fa"
#     shell:
#         "module use --prepend /data/CDSL_Gurobi_users/modules && "
#         "module load gurobi && "
#         "export GRB_LICENSE_FILE=/data/CDSL_Gurobi_users/gurobi910/gurobi.lic && "
#         "hatchet run hatchet.ini"

rule download_bam_files:
    output:
        "data/normal.bam",
        "data/normal.bam.bai",
        "data/bulk_03clone1_06clone0_01normal.sorted.bam",
        "data/bulk_03clone1_06clone0_01normal.sorted.bam.bai",
        "data/bulk_08clone1_Noneclone0_02normal.sorted.bam",
        "data/bulk_08clone1_Noneclone0_02normal.sorted.bam.bai",
        "data/bulk_Noneclone1_09clone0_01normal.sorted.bam",
        "data/bulk_Noneclone1_09clone0_01normal.sorted.bam.bai"
    shell:
        "curl -L 'https://zenodo.org/record/4046906/files/normal.bam?download=1' > data/normal.bam && "
        "curl -L 'https://zenodo.org/record/4046906/files/normal.bam.bai?download=1' > data/normal.bam.bai && "
        "curl -L 'https://zenodo.org/record/4046906/files/bulk_03clone1_06clone0_01normal.sorted.bam?download=1' > data/bulk_03clone1_06clone0_01normal.sorted.bam && "
        "curl -L 'https://zenodo.org/record/4046906/files/bulk_03clone1_06clone0_01normal.sorted.bam.bai?download=1' > data/bulk_03clone1_06clone0_01normal.sorted.bam.bai && "
        "curl -L 'https://zenodo.org/record/4046906/files/bulk_08clone1_Noneclone0_02normal.sorted.bam?download=1' > data/bulk_08clone1_Noneclone0_02normal.sorted.bam && "
        "curl -L 'https://zenodo.org/record/4046906/files/bulk_08clone1_Noneclone0_02normal.sorted.bam.bai?download=1' > data/bulk_08clone1_Noneclone0_02normal.sorted.bam.bai && "
        "curl -L 'https://zenodo.org/record/4046906/files/bulk_Noneclone1_09clone0_01normal.sorted.bam?download=1' > data/bulk_Noneclone1_09clone0_01normal.sorted.bam && "
        "curl -L 'https://zenodo.org/record/4046906/files/bulk_Noneclone1_09clone0_01normal.sorted.bam.bai?download=1' > data/bulk_Noneclone1_09clone0_01normal.sorted.bam.bai"


rule download_reference_genome:
    output:
        "data/hg19.dict",
        "data/hg19.fa"
    shell:
        "module load samtools && "
        "curl -L https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.fa.gz | gzip -d > data/hg19.fa && "
        "samtools faidx data/hg19.fa && "
        "samtools dict data/hg19.fa > data/hg19.dict"
