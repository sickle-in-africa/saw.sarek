process MultiQC {
    publishDir "${params.outdir}/Reports/MultiQC", mode: params.publish_dir_mode

    input:
        file (multiqcConfig) from ch_multiqc_config
        file (mqc_custom_config) from ch_multiqc_custom_config.collect().ifEmpty([])
        file (versions) from ch_software_versions_yaml.collect()
        file workflow_summary from ch_workflow_summary.collectFile(name: "workflow_summary_mqc.yaml")
        file ('bamQC/*') from bamQCReport.collect().ifEmpty([])
        file ('BCFToolsStats/*') from bcftoolsReport.collect().ifEmpty([])
        file ('FastQC/*') from fastQCReport.collect().ifEmpty([])
        file ('TrimmedFastQC/*') from trimGaloreReport.collect().ifEmpty([])
        file ('MarkDuplicates/*') from duplicates_marked_report.collect().ifEmpty([])
        file ('DuplicatesMarked/*.recal.table') from baseRecalibratorReport.collect().ifEmpty([])
        file ('SamToolsStats/*') from samtoolsStatsReport.collect().ifEmpty([])
        file ('snpEff/*') from snpeffReport.collect().ifEmpty([])
        file ('VCFTools/*') from vcftoolsReport.collect().ifEmpty([])

    output:
        file "*multiqc_report.html" into ch_multiqc_report
        file "*_data"
        file "multiqc_plots"

    when: !('multiqc' in skipQC)

    script:
    rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
    rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
    custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
    """
    multiqc -f ${rtitle} ${rfilename} ${custom_config_file} .
    """
}

ch_multiqc_report.dump(tag:'MultiQC')

// Output Description HTML
process Output_documentation {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode

    input:
        file output_docs from ch_output_docs
        file images from ch_output_docs_images

    output:
        file "results_description.html"

    when: !('documentation' in skipQC)

    script:
    """
    markdown_to_html.py $output_docs -o results_description.html
    """
}