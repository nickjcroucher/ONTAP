process PYTHON_COVERAGE_OVER_DEFINED_REGIONS {
    label 'cpu_1'
    label 'mem_1'
    label 'time_1'

    conda "conda-forge::pandas=2.2.1"
    container "quay.io/sangerpathogens/pandas:2.2.1"

    // Static publish path to avoid directive-scope issues with qc_stage
    publishDir "${params.outdir}/qc/coverage/coverage_summary", mode: 'copy', overwrite: true

    input:
    tuple val(meta), path(samtools_coverage), path(target_regions_bed)
    val(qc_stage)

    output:
    tuple val(meta), path(coverage_summary), emit: coverage_summary

    script:
    coverage_summary = "*coverage_summary.tsv"
    """
    coverage_over_defined_regions.py \
        -s ${samtools_coverage} \
        -b ${target_regions_bed} \
        -t ${params.coverage_reporting_thresholds},${params.coverage_filtering_threshold} \
        -n ${meta.ID}
    """
}

process PYTHON_PLOT_COVERAGE {
    label 'cpu_1'
    label 'mem_1'
    label 'time_30m'

    conda "plotly::plotly=5.21.0"
    container "quay.io/sangerpathogens/python_graphics:1.0.0"

    // Static publish path to avoid directive-scope issues with qc_stage
    publishDir "${params.outdir}/qc/coverage/coverage_summary", mode: 'copy', overwrite: true

    input:
    path(coverage_summaries)
    val(qc_stage)

    output:
    path("*.html"), emit: coverage_plot

    script:
    """
    plot_coverage.py .
    """
}