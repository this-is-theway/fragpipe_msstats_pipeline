/*
========================================================================================
    FragPipe to R Analysis Pipeline
========================================================================================
    A Nextflow pipeline for proteomics data analysis using FragPipe
    followed by statistical analysis and HTML report generation in R
----------------------------------------------------------------------------------------
*/

nextflow.enable.dsl = 2

/*
========================================================================================
    INCLUDE MODULES
========================================================================================
*/

include { RAW_TO_MZML } from './modules/fragpipe'
include { GENERATE_MANIFEST } from './modules/fragpipe'
include { FRAGPIPE          } from './modules/fragpipe'
include { MSSTATS_REPORT    } from './modules/msstats'



/*
========================================================================================
    MAIN WORKFLOW
========================================================================================
*/

workflow {

    log.info """\

    ===========================================================
               F R A G P I P E  +  M S S T A T S  
    -----------------------------------------------------------
              P R O T E O M I C S  A N A L Y S I S 
    ===========================================================
    input                : ${params.input}
    fasta                : ${params.fasta}
    workflow             : ${params.fragpipe_workflow}
    experimental design  : ${params.experimental_design}
    comparisons          : ${params.comparisons_file}
    outdir               : ${params.outdir}
    threads              : ${params.fragpipe_threads}
    ram (GB)             : ${params.fragpipe_ram}
    ===========================================================
    """
    .stripIndent()

        // -------------------------------------------------------------------
    // Validate required parameters
    // -------------------------------------------------------------------
    if (!params.input) {
        error "Please specify input MS files with --input (glob pattern)"
    }
    if (!params.fasta) {
        error "Please specify FASTA database with --fasta"
    }
    if (!params.fragpipe_workflow) {
        error "Please specify a FragPipe .workflow file with --fragpipe_workflow"
    }
    if (!params.experimental_design) {
        error "Please specify experimental design CSV with --experimental_design"
    }
    if (!params.rmd_template) {
        error "Please specify an RMarkdown template with --rmd_template"
    }
    if (!params.scripts) {
        error "Please specify R scripts directory with --scripts"
    }
    if (!params.comparisons_file) {
        error "Please specify comparisons CSV with --comparisons_file"
    }

    // -------------------------------------------------------------------
    // Create input channels
    // -------------------------------------------------------------------

    // Experimental design, Rmd template, and R scripts as separate channels
    design_ch       = channel.fromPath(params.experimental_design, checkIfExists: true)
    rmd_ch          = channel.fromPath(params.rmd_template, checkIfExists: true)
    scripts_ch      = channel.fromPath(params.scripts, checkIfExists: true)
    comparisons_ch  = channel.fromPath(params.comparisons_file, checkIfExists: true)

    // Highlight genes file (optional — create empty file if not provided)
    highlight_ch = params.highlight_genes_file
        ? channel.fromPath(params.highlight_genes_file, checkIfExists: true)
        : channel.of(file("NO_HIGHLIGHT")).map { file -> 
            def f = file("${workDir}/empty_highlights.txt")
            f.text = ""
            return f
        }
    fasta_ch    = channel.fromPath(params.fasta, checkIfExists: true)
    
    if (params.msstats_csv) {
        // ---- MSstats-only: reuse a pre-computed FragPipe result ----
        log.info "MSstats-only mode: using ${params.msstats_csv} (skipping FragPipe)"
        msstats_ch = channel.fromPath(params.msstats_csv, checkIfExists: true)
    }
    else {
                    
                // input files are collected and converted to mzML if .raw format
                input_ch = channel.fromPath(params.input, type:'any')
                RAW_TO_MZML(input_ch)

                // MS files — collected into a single list so they are all staged together
                ms_files_ch = RAW_TO_MZML.out.mzml
                    .mix(input_ch.filter { !it.name.endsWith('.raw') })
                    .collect()

                workflow_ch = channel.fromPath(params.fragpipe_workflow, checkIfExists: true)

            

                // -------------------------------------------------------------------
                // STEP 0 (conditional): Generate manifest if not supplied
                // -------------------------------------------------------------------
                if (params.fragpipe_manifest) {
                    manifest_ch = channel.fromPath(params.fragpipe_manifest, checkIfExists: true)
                } else {
                    GENERATE_MANIFEST(ms_files_ch, design_ch)
                    manifest_ch = GENERATE_MANIFEST.out.manifest
                }

                // -------------------------------------------------------------------
                // STEP 1: Run FragPipe
                // -------------------------------------------------------------------
                FRAGPIPE(
                    ms_files_ch,
                    fasta_ch,
                    workflow_ch,
                    manifest_ch
                )

                msstats_ch = FRAGPIPE.out.msstats_file
    }

    // -------------------------------------------------------------------
    // STEP 2: MSstats analysis & HTML report
    // -------------------------------------------------------------------
    MSSTATS_REPORT(
        msstats_ch,
        rmd_ch,
        scripts_ch,
        comparisons_ch,
        highlight_ch,
        fasta_ch
    )

    // -------------------------------------------------------------------
    // Completion handler
    // -------------------------------------------------------------------
    workflow.onComplete = {
        println "Pipeline completed at : ${workflow.complete}"
        println "Execution status      : ${workflow.success ? 'Success' : 'Failed'}"
        println "Results directory     : ${params.outdir}"
    }
}
