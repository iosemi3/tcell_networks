rule all:
    input:
        "results/top10_markers_per_cluster.csv",
        "results/top10_markers_per_identity.csv",


rule run_step1:
    input:
        input_feature_matrix="/corgi/sebas/tcell_multi/aggregated_tlibs/AGG6789_good/outs/filtered_feature_bc_matrix.h5",
        input_fragments="/corgi/sebas/tcell_multi/aggregated_tlibs/AGG6789_good/outs/atac_fragments.tsv.gz",
    output:
        pmbc_rdata_file="results/pmbc_1103.rds",
    script:
        "step1_create_joint_basic_analysis.R"


rule run_step2:
    input:
        seurat_object="results/pmbc_1103.rds",
    output:
        cluster_markers="results/top10_markers_per_cluster.csv",
        identity_markers="results/top10_markers_per_identity.csv",
    script:
        "step2_get_topDE.R"
