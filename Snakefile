rule run_step1:
  input:
    input_feature_matrix = "/corgi/sebas/tcell_multi/aggregated_tlibs/AGG6789_good/outs/filtered_feature_bc_matrix.h5",
    input_fragments="/corgi/sebas/tcell_multi/aggregated_tlibs/AGG6789_good/outs/atac_fragments.tsv.gz"
  output:
    pmbc_rdata_file = "results/pmbc_1103.rds"
  script:
    "step1_create_joint_basic_analysis.R"
