process GetFastqcQualityReport {
    label 'FastQC'
    label 'cpus_2'

    tag "${idPatient}-${idRun}"

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        tuple val(idPatient), val(idSample), val(idRun), path("${idSample}_${idRun}_R1.fastq.gz"), path("${idSample}_${idRun}_R2.fastq.gz")

    output:
        path("*.{html,zip}")

    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}_R1.fastq.gz ${idSample}_${idRun}_R2.fastq.gz
    """
}

process GetUnmappedBamQualityReport {
    label 'FastQC'
    label 'cpus_2'

    tag "${idPatient}-${idRun}"

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        tuple val(idPatient), val(idSample), val(idRun), path("${idSample}_${idRun}.bam")

    output:
        path("*.{html,zip}")

    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}.bam
    """
}

process GetCohortRawReadsQualityReport {
    publishDir "${params.outdir}/Reports/MultiQC", mode: params.publish_dir_mode

    input:
        file (multiqcConfig)
        file (versions)
        file ('FastQC/*')

    script:
        rtitle = custom_runName ? "--title \"$custom_runName\"" : ''
        rfilename = custom_runName ? "--filename " + custom_runName.replaceAll('\\W','_').replaceAll('_+','_') + "_multiqc_report" : ''
        custom_config_file = params.multiqc_config ? "--config $mqc_custom_config" : ''
        """
        multiqc -f ${rtitle} ${rfilename} ${custom_config_file} .
        """
}