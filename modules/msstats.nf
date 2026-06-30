/*
========================================================================================
    MSstats Report Process Module
========================================================================================
    Runs an RMarkdown-based MSstats analysis on FragPipe output.
----------------------------------------------------------------------------------------
*/

process MSSTATS_REPORT {
    tag "MSstats Analysis"
    label 'process_medium'
    publishDir "${params.outdir}/msstats",      mode: 'copy', pattern: '*.{html,xlsx}'
    publishDir "${params.outdir}/msstats/logs",  mode: 'copy', pattern: '*.log'

    container 'bioconductor/bioconductor_docker:devel'

    input:
    path msstats_file          // msstats.csv from FragPipe
    path rmd_template          // run_msstats.Rmd
    path scripts               // scripts/ directory (R scripts + comparison_template.Rmd)
    path logo_path
    path comparisons_file      // comparisons.csv (condition_A, condition_B)
    path highlight_genes_file  // genes_to_highlight.txt (optional, can be empty)
    path fasta

    output:
    path "msstats_report.html",      emit: report
    path "*.xlsx",                   emit: excel, optional: true
    path "*.log",                    emit: logs,  optional: true

    script:
    """
    Rscript -e "rmarkdown::render(
        input       = '${rmd_template}',
        output_file = 'msstats_report.html',
        output_dir  = '\$PWD',
        params = list(
            input_file           = '\$PWD/${msstats_file}',
            comparisons_file     = '\$PWD/${comparisons_file}',
            highlight_genes_file = '\$PWD/${highlight_genes_file}',
            scripts_dir          = '\$PWD/${scripts}',
            logo_path            = '\$PWD/${logo_path}',
            fasta_file           = '\$PWD/${fasta}',
            imputation_method    = '${params.imputation_method}',
            pvalue_cutoff        = ${params.fdr_threshold},
            fc_cutoff            = ${params.log2fc_threshold},
            report_title         = '${params.report_title}',
            report_author        = '${params.report_author}'
        )
    )"
    """
}
