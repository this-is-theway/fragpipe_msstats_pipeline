rmarkdown::render(
  input       = '/home/jszymanski/Documents/nextflow_pipelines_proteomics/fragpipe_msstats_pipeline/report/run_msstats.Rmd',
  output_file = 'msstats_report_ludovic.html',
  output_dir  = './results',
  params = list(
    input_file           = '/home/jszymanski/Documents/BF_projects/BF_2026_004_Szymanska_MS_data_analysis/results/fragpipe/fragpipe_output/msstats.csv',
    comparisons_file     = '/home/jszymanski/Documents/BF_projects/BF_2026_004_Szymanska_MS_data_analysis/design/comparisons.csv',
    highlight_genes_file = '/home/jszymanski/Documents/BF_projects/BF_2026_004_Szymanska_MS_data_analysis/design/genes_to_highlight.txt',
    fasta_file           = '/home/jszymanski/Documents/nextflow_pipelines_proteomics/fasta_databases/2026-03-18-decoys-contam-uniprotkb_proteome_UP000005640_AND_revi_2026_03_18.fasta.fas',
    scripts_dir          = '/home/jszymanski/Documents/nextflow_pipelines_proteomics/fragpipe_msstats_pipeline/report/scripts/',
    imputation_method    = 'ludovic',
    pvalue_cutoff        = 0.05,
    fc_cutoff            = 1,
    report_title         = 'MS Poject Analysis',
    report_author        = 'Jacek Szymanski'
  )
)


