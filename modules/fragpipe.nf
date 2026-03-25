/*
========================================================================================
    FragPipe Process Module
========================================================================================
    Runs FragPipe headless for mass-spectrometry-based proteomics analysis.
    Includes an optional manifest generator when no pre-made manifest is supplied.
----------------------------------------------------------------------------------------
*/


/*
 * GENERATE_MANIFEST
 * -----------------
 * Builds a .fp-manifest from the staged MS files and an experimental-design CSV.
 *
 * Expected CSV columns (at minimum):
 *   file, experiment, bioreplicate, data_type
 *
 * Where:
 *   file         = original filename (basename, e.g. sample1.mzML)
 *   experiment   = experiment/group label
 *   bioreplicate = integer replicate id
 *   data_type    = DDA or DIA (default DDA)
 */



process RAW_TO_MZML {
    tag "${raw_file.simpleName}"
    label 'process_medium'

    input:
    path raw_file

    output:
    path "*.mzML", emit: mzml

    when:
    raw_file.name.endsWith('.raw')

    script:
    """
    mono ${params.thermorawfileparser} \
        -i=${raw_file} \
        -o=. \
        -f=2
    """
}

process GENERATE_MANIFEST {
    tag "Generate manifest"
    label 'process_low'

    input:
    path ms_files            // staged MS files (collected list)
    path experimental_design // CSV with columns: file, experiment, bioreplicate, data_type

    output:
    path "experiment.fp-manifest", emit: manifest

    script:
    """
    #!/usr/bin/env python3
    import csv, os, sys

    pwd = os.getcwd()
    staged = {f: os.path.join(pwd, f) for f in os.listdir('.')
            if f.endswith(('.mzML', '.raw', '.d', '.mzXML'))}

    with open('${experimental_design}') as fh:
        reader = csv.DictReader(fh)
        rows = list(reader)

    with open('experiment.fp-manifest', 'w') as out:
        for row in rows:
            fname = row['file'].strip()
            if fname in staged:
                path = staged[fname]
            elif fname.endswith('.raw'):
                converted = fname.replace('.raw', '.mzML')
                if converted in staged:
                    path = staged[converted]
                else:
                    print(f"WARNING: {fname} not found (also tried {converted}) — skipping", file=sys.stderr)
                    continue
            else:
                print(f"WARNING: {fname} not staged — skipping", file=sys.stderr)
                continue
            exp    = row.get('experiment', 'default').strip()
            biorep = row.get('bioreplicate', '1').strip()
            dtype  = row.get('data_type', 'DDA').strip()
            out.write(f"{path}\\t{exp}\\t{biorep}\\t{dtype}\\n")

    print("Manifest written with", sum(1 for _ in open('experiment.fp-manifest')), "entries")
    """
}



process FRAGPIPE {
    tag "FragPipe Analysis"
    label 'process_high'
    publishDir "${params.outdir}/fragpipe", mode: 'copy'

    // Uncomment / adjust for your environment:
    // container 'fragpipe:latest'
    // conda 'fragpipe'

    input:
    path ms_files
    path fasta
    path workflow
    path manifest

    output:
    path "fragpipe_output/**",                   emit: results
    path "fragpipe_output{,/**}/msstats.csv",     emit: msstats_file
    path "fragpipe_output/*.log",              emit: logs, optional: true

    script:
    def mem  = params.fragpipe_ram    ?: (task.memory ? task.memory.toGiga() : 0)
    def cpus = params.fragpipe_threads ?: (task.cpus   ?: 1)
    """
    mkdir -p fragpipe_output

    # ---- rewrite manifest so every path points at the Nextflow work dir ----
    awk -F'\\t' -v pwd="\$PWD" 'BEGIN{OFS="\\t"} {
        n = split(\$1, a, "/")
        \$1 = pwd "/" a[n]
        print
    }' ${manifest} > updated_manifest.fp-manifest

     # ---- inject FASTA path into workflow (FragPipe GUI adds this separately) ----
    FASTA_PATH="\$PWD/${fasta}"

    # Insert database.db-path at the top (the workflow file doesn't contain it)
    echo "database.db-path=\${FASTA_PATH}" > updated.workflow
    
    # Add to your updated.workflow generation
    echo "fragpipe-config.bin-thermorawfileparser=${params.thermorawfileparser}" >> updated.workflow
    echo "fragpipe-config.bin-python=${params.python_bin}" >> updated.workflow
    cat ${workflow} >> updated.workflow

    # Verify
    echo "=== FASTA path in workflow ==="
    grep 'database.db-path' updated.workflow
    
    # ---- run FragPipe headless ----
    ${params.fragpipe_bin}/fragpipe \\
        --headless \\
        --workflow updated.workflow \\
        --manifest updated_manifest.fp-manifest \\
        --workdir  fragpipe_output \\
        --ram      ${mem} \\
        --threads  ${cpus}
    """
}
