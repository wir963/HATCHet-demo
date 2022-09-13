

PLOTS_DIR = join("output", "plots", "{phasing}-snps")

BIN_2D_PLOT = join(PLOTS_DIR, "bb.pdf")
CLUSTERED_BIN_2D_PLOT = join(PLOTS_DIR, "bb_clustered.pdf")
BAF_PLOT = join(PLOTS_DIR, "ballelefrequency.pdf")
CLUSTERED_BAF_PLOT = join(PLOTS_DIR, "ballelefrequency_clustered.pdf")
RDR_PLOT = join(PLOTS_DIR, "readdepthratio.pdf")
CLUSTERED_RDR_PLOT = join(PLOTS_DIR, "readdepthratio_clustered.pdf")


rule hatchet_plot_cn:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        "output/results/best.bbc.ucn"
    output:
        "output/summary/intratumor-clones-allelecn.pdf"
    shell:
        "hatchet plot-cn "
        "--rundir output/summary "
        "{input} "


# 2d-scatter plots where x-axis corresponds to the BAF and y-axis to the RDR
rule hatchet_plot_bins_BB:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        clustered_genomic_bins
    params:
        PLOTS_DIR
    output:
        BIN_2D_PLOT
    shell:
        "hatchet plot-bins "
        "-c BB "
        "--rundir {params} "
        "{input} "

# 2d-scatter plots where x-axis corresponds to the BAF and y-axis to the RDR
# points are colored by cluster (as determined by HATCHet) - need to pass clusters in here
rule hatchet_plot_bins_CBB:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        clustered_genomic_bins
    params:
        PLOTS_DIR
    output:
        CLUSTERED_BIN_2D_PLOT
    shell:
        "hatchet plot-bins "
        "-c CBB "
        "--rundir {params} "
        "{input} "

rule hatchet_plot_bins_BAF:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        clustered_genomic_bins
    params:
        PLOTS_DIR
    output:
        BAF_PLOT
    shell:
        "hatchet plot-bins "
        "-c BAF "
        "--rundir {params} "
        "{input} "

rule hatchet_plot_bins_CBAF:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        clustered_genomic_bins
    params:
        PLOTS_DIR
    output:
        CLUSTERED_BAF_PLOT
    shell:
        "hatchet plot-bins "
        "-c CBAF "
        "--rundir {params} "
        "{input} "

rule hatchet_plot_bins_RDR:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        clustered_genomic_bins
    params:
        PLOTS_DIR
    output:
        RDR_PLOT
    shell:
        "hatchet plot-bins "
        "-c RD "
        "--rundir {params} "
        "{input} "

rule hatchet_plot_bins_RDR_clustered:
    # conda:
    #     "envs/HATCHet-env.yaml"
    input:
        clustered_genomic_bins
    params:
        PLOTS_DIR
    output:
        RDR_CLUSTERED_PLOT
    shell:
        "hatchet plot-bins "
        "-c CRD "
        "--rundir {params} "
        "{input} "
