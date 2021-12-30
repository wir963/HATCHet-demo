

rule all:
    input:
        "data/normal.bam",
        "data/hg19.fa"

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

# download these common snps because this reference genome uses chr1, etc.
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/00-All.vcf.gz
# alternatively, you would download these common snps
# https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/00-All.vcf.gz
rule download_common_snps:
    shell:
        "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/00-All.vcf.gz"

# rule genotype_snps:
#     conda:
#         "envs/HATCHet-env.yaml"
#     shell:
#         "hatchet run hatchet-genotype-snps.ini"

rule genotype_snps:
    conda:
        "envs/HATCHet-env.yaml"
    input:
        "data/normal.bam",
        "data/hg19.fa"
    output:
        directory("new_output/snps/")
    shell:
        "hatchet hatchet-genotype "
        "--normal data/normal.bam "
        "--reference data/hg19.fa "
        "--snps https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/00-All.vcf.gz "# use this version because the ref genome using chr1
        "--mincov 8 "
        "--maxcov 300 "
        "--outputsnps new_output/snps/"

# rule hatchet_count_alleles:
#     conda:
#         "envs/HATCHet-env.yaml"
#     shell:
#         "hatchet run hatchet-count-alleles.ini"


rule hatchet_count_alleles:
    conda:
        "envs/HATCHet-env.yaml"
    output:
        normal="new_output/baf/normal.1bed",
        tumor="new_output/baf/tumor.1bed"
    shell:
        "hatchet count-alleles "
        "--tumors data/bulk_03clone1_06clone0_01normal.sorted.bam data/bulk_08clone1_Noneclone0_02normal.sorted.bam data/bulk_Noneclone1_09clone0_01normal.sorted.bam "
        "--normal data/normal.bam "
        "--reference data/hg19.fa "
        "--snps output/snps/*.vcf.gz "
        "--outputnormal {output.normal} "
        "--outputtumors {output.tumor} "
        # "--outputsnps "
        # "--regions " required for WXS (although I don't think this is true, especially since SNPs are passed)
        "--mincov 8 "
        "--maxcov 300 "

rule hatchet_count_reads:
    conda:
        "envs/HATCHet-env.yaml"
    input:
        ""
    output:
        "output/rdr/normal.1bed",
        "output/rdr/tumor.1bed",
        "output/rdr/total.tsv"
    shell:
        "hatchet count-reads "
        "--normal"

rule hatchet_combine_counts:
    conda:
        "envs/HATCHet-env.yaml"
    output:
        "output/rdr/total.tsv"
    shell:
        "hatchet run hatchet-combine-counts.ini"

rule hatchet_cluster_bins:
    conda:
        "envs/HATCHet-env.yaml"
    shell:
        "hatchet run hatchet-cluster-bins.ini"

rule hatchet_plot_bins:
    conda:
        "envs/HATCHet-env.yaml"
    shell:
        "hatchet run hatchet-plot-bins.ini"

# let's genotype snps using something like platypus
# rule genotype_snps:
#     conda:
#         "envs/HATCHet-env.yaml"
#     shell:
#         "hatchet genotype-snps "
#         "-N {NORMAL} "
#         "-r {REF} "
#         "-j ${J} "
#         "-c ${MINREADS} "
#         "-C ${MAXREADS}"
#         # "-R ${LIST} " use this to limit genotype_snps to a particular region
#         "-o ${SNP}"
