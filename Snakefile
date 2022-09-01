

# rule all:
#     input:
#         #"new_output/snps/",
#         "new_output/baf/normal.1bed"
#         # "data/normal.bam",
#         # "data/hg19.fa"

rule run_hatchet_init:
    conda:
        "envs/HATCHet-env.yaml"
    input:
        "data/normal.bam",
        "data/normal.bam.bai",
        "data/bulk_03clone1_06clone0_01normal.sorted.bam",
        "data/bulk_03clone1_06clone0_01normal.sorted.bam.bai",
        "data/bulk_08clone1_Noneclone0_02normal.sorted.bam",
        "data/bulk_08clone1_Noneclone0_02normal.sorted.bam.bai",
        "data/bulk_Noneclone1_09clone0_01normal.sorted.bam",
        "data/bulk_Noneclone1_09clone0_01normal.sorted.bam.bai",
        "data/hg19.dict",
        "data/hg19.fa"
    shell:
        "module use --prepend /data/CDSL_Gurobi_users/modules && "
        "module load gurobi && "
        "export GRB_LICENSE_FILE=/data/CDSL_Gurobi_users/gurobi910/gurobi.lic && "
        "hatchet run hatchet.ini"

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
# don't need to actually download because you can use a URL in the genotype_snps rule
# rule download_common_snps:
#     shell:
#         "https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/00-All.vcf.gz"

# rule genotype_snps:
#     conda:
#         "envs/HATCHet-env.yaml"
#     shell:
#         "hatchet run hatchet-genotype-snps.ini"

# rule genotype_snps:
#     conda:
#         "envs/HATCHet-env.yaml"
#     input:
#         "data/normal.bam",
#         "data/hg19.fa"
#     output:
#         directory("new_output/snps/")
#     shell:
#         "mkdir {output[0]} && "
#         "hatchet genotype-snps "
#         "--normal data/normal.bam "
#         "--reference data/hg19.fa "
#         "--snps https://ftp.ncbi.nih.gov/snp/organisms/human_9606_b151_GRCh37p13/VCF/GATK/00-All.vcf.gz "# use this version because the ref genome using chr1
#         "--mincov 8 "
#         "--maxcov 300 "
#         "--outputsnps new_output/snps/"

# rule hatchet_count_alleles:
#     conda:
#         "envs/HATCHet-env.yaml"
#     shell:
#         "hatchet run hatchet-count-alleles.ini"

#
# rule hatchet_count_alleles:
#     conda:
#         "envs/HATCHet-env.yaml"
#     input:
#         normal="data/normal.bam",
#         tumor=["data/bulk_03clone1_06clone0_01normal.sorted.bam", "data/bulk_08clone1_Noneclone0_02normal.sorted.bam", "data/bulk_Noneclone1_09clone0_01normal.sorted.bam"],
#     output:
#         normal_baf="new_output/baf/normal.1bed", # contains the number of reads for the major and minor allele
#         tumor_baf="new_output/baf/tumor.1bed",
#         snps_list = directory("new_output/count_alleles/snps")
#     shell:
#         "mkdir {output.snps_list} && "
#         #"touch {output[0]} && "
#         #"touch {output[1]} && "
#         "hatchet count-alleles "
#         "--tumors {input.tumor} "
#         "--normal {input.normal} "
#         "--samples normal tumor1 tumor2 tumor3 "
#         "--reference data/hg19.fa "
#         "--snps new_output/snps/*.vcf.gz "
#         "--outputnormal {output.normal_baf} "
#         "--outputtumors {output.tumor_baf} "
#         "--outputsnps {output.snps_list} "
#         # "--regions " required for WXS (although I don't think this is true, especially since SNPs are passed)
#         "--mincov 8 "
#         "--maxcov 300 "
#         #"|| true"
#
# # count the mapped sequencing reads in bins of fixed and given length, uniformly for a BAM file of a normal sample
# # and one or more BAM files of tumor samples
# rule hatchet_count_reads:
#     conda:
#         "envs/HATCHet-env.yaml"
#     input:
#         normal="data/normal.bam",
#     params: # use as params for now TODO figure out better way
#         tumor="data/bulk_03clone1_06clone0_01normal.sorted.bam data/bulk_08clone1_Noneclone0_02normal.sorted.bam data/bulk_Noneclone1_09clone0_01normal.sorted.bam",
#     output:
#         normal_rdr="output/rdr/normal.1bed",
#         tumor_rdr="output/rdr/tumor.1bed",
#         total_rd="output/rdr/total.tsv"
#     shell:
#         "hatchet count-reads "
#         "--tumors {output.tumor} "
#         "--normal {input.normal} "
#         "--samples normal tumor1 tumor2 tumor3 "
#         "--reference data/hg19.fa "
#         "--size 50kb "
#         "--processes 6 "
#         "--outputnormal {output.normal_rdr} "
#         "--outputtumors {output.tumor_rdr} "
#         "--outputtotal {output.total_rd} "
#
# # combine tumor bin counts, normal bin counts and tumor allele counts to obtain the read-depth ratio
# # and the mean B-allele frequency of each bin
# rule hatchet_combine_counts:
#     conda:
#         "envs/HATCHet-env.yaml"
#     input:
#         normal_rdr="output/rdr/normal.1bed",
#         tumor_rdr="output/rdr/tumor.1bed",
#         total_rd="output/rdr/total.tsv",
#         tumor_baf="output/baf/tumor.1bed"
#     output:
#         "output/bb/bulk.bb" # tab separated file that include many fields see - http://compbio.cs.brown.edu/hatchet/doc_combine_counts.html
#     shell:
#         "hatchet combine-counts "
#         "--normalbins {input.normal_rdr} "
#         "--tumorbins {input.tumor_rdr} "
#         "--tumorbafs {input.tumor_baf} "
#         "--totalcounts {input.total_rd} "
#         "--blocklength 25kb "
#         "--phase None "
#         " > {output} "
#         # "--msr 3000 " # minimum number of SNP-covering reads
#         # "--mtr 5000 " # minimum number of total reads per bin
#         # "--outfile {output} "
#
# rule hatchet_cluster_bins:
#     conda:
#         "envs/HATCHet-env.yaml"
#     input:
#         "output/bb/bulk.bb"
#     output:
#         "output/bbc/bulk.bbc",
#         "output/bbc/bulk.seg"
#     shell:
#         "hatchet cluster-bins "
#         "{input} "
#         "--outsegments {output[1]} "
#         "--outbins {output[0]} "
#         "--diploidbaf 0.08 "
#         "--tolerancebaf 0.04 "
#         "--tolerancerdr 0.15 "
#
#
# rule hatchet_plot_bins:
#     conda:
#         "envs/HATCHet-env.yaml"
#     input:
#         "output/bbc/bulk.bbc"
#     params:
#         "output/plots"
#     output:
#         "output/plots/bb_clustered.png"
#     shell:
#         "hatchet plot-bins "
#         "--rundir {params} "
#         "-tS 0.005 "
#         "{input} "
#
# rule hatchet_compute_cn:
#     conda:
#         "envs/HATCHet-env.yaml"
#     input:
#         "output/bbc/bulk.bbc",
#         "output/bbc/bulk.seg"
#     params:
#         in_dir="output/bbc/bulk",
#         out_dir="output/results"
#     output:
#         "output/results/best.bbc.ucn",
#         "output/results/best.seg.ucn",
#     shell:
        # "module use --prepend /data/CDSL_Gurobi_users/modules && "
        # "module load gurobi && "
        # "export GRB_LICENSE_FILE=/data/CDSL_Gurobi_users/gurobi910/gurobi.lic && "
#         "hatchet compute-cn "
#         "-i {params.in_dir} "
#         "-x {params.out_dir} "
#         "--clones 2,8 " # default
#
# rule hatchet_plot_cn:
#     conda:
#         "envs/HATCHet-env.yaml"
#     input:
#         "output/results/best.bbc.ucn"
#     output:
#         "output/summary/intratumor-clones-allelecn.pdf"
#     shell:
#         "hatchet plot-cn "
#         "--rundir output/summary "
#         "{input} "


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
