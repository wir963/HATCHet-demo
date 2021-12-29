

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

rule genotype_snps:
    conda:
        "envs/HATCHet-env.yaml"
    shell:
        "hatchet run hatchet.ini"
