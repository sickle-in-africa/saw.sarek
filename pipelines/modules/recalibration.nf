include {
    isChannelActive
} from "${params.modulesDir}/sarek.nf"

process GetBaseRecalibrationReport {
    label 'cpus_1'

    tag "${idPatient}-${idSample}-${intervalBed.baseName}"

    input:
        tuple val(idPatient), val(idSample), file(bam), file(bai), file(intervalBed)
        file(dbsnp)
        file(dbsnpIndex)
        file(fasta)
        file(dict)
        file(fastaFai)
        file(knownIndels)
        file(knownIndelsIndex)

    output:
        tuple val(idPatient), val(idSample), file("${prefix}${idSample}.recal.table")

    when: (params.skip_baserecalibration == false)

    script:
    dbsnpOptions = isChannelActive(dbsnp) ? "--known-sites ${dbsnp}" : ""
    knownOptions = isChannelActive(knownIndels) ? knownIndels.collect{"--known-sites ${it}"}.join(' ') : ""
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    // TODO: --use-original-qualities ???
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        BaseRecalibrator \
        -I ${bam} \
        -O ${prefix}${idSample}.recal.table \
        --tmp-dir . \
        -R ${fasta} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        ${knownOptions} \
        --verbosity INFO
    """
}

process MergeBaseRecalibrationReportsForSample {
    label 'memory_singleCPU_2_task'
    label 'cpus_2'

    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${idSample}.recal.table" && !params.skip_markduplicates) "Preprocessing/${idSample}/DuplicatesMarked/${it}"
            else "Preprocessing/${idSample}/Mapped/${it}"
        }

    input:
        tuple val(idPatient), val(idSample), file(recal)

    output:
        tuple val(idPatient), val(idSample), file("${idSample}.recal.table")
        tuple val(idPatient), val(idSample)

    script:
    input = recal.collect{"-I ${it}"}.join(' ')
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GatherBQSRReports \
        ${input} \
        -O ${idSample}.recal.table \
    """
}


process RecalibrateBasesInReadGroup {
    label 'memory_singleCPU_2_task'
    label 'cpus_2'

    tag "${idPatient}-${idSample}-${intervalBed.baseName}"

    input:
        tuple val(idPatient), val(idSample), file(bam), file(bai), file(recalibrationReport), file(intervalBed)
        file(dict)
        file(fasta)
        file(fastaFai)

    output:
        tuple val(idPatient), val(idSample), file("${prefix}${idSample}.recal.bam")

    script:
    prefix = params.no_intervals ? "" : "${intervalBed.baseName}_"
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        ApplyBQSR \
        -R ${fasta} \
        --input ${bam} \
        --output ${prefix}${idSample}.recal.bam \
        ${intervalsOptions} \
        --bqsr-recal-file ${recalibrationReport}
    """
}

process MergeRecalibratedReadGroupsForSample {
    label 'cpus_8'

    tag "${idPatient}-${idSample}"

    publishDir "${params.outdir}/Preprocessing/${idSample}/Recalibrated", mode: params.publish_dir_mode

    input:
        tuple val(idPatient), val(idSample), file(bam)

    output:
        tuple val(idPatient), val(idSample), file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai")
        tuple val(idPatient), val(idSample), file("${idSample}.recal.bam")
        tuple val(idPatient), val(idSample)

    when: !(params.no_intervals)

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.recal.bam ${bam}
    samtools index ${idSample}.recal.bam
    """
}

process IndexRecalibratedSampleReadGoup {
    label 'cpus_8'

    tag "${idPatient}-${idSample}"

    publishDir "${params.outdir}/Preprocessing/${idSample}/Recalibrated", mode: params.publish_dir_mode

    input:
        tuple val(idPatient), val(idSample), file("${idSample}.recal.bam")

    output:
        tuple val(idPatient), val(idSample), file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai")
        tuple val(idPatient), val(idSample), file("${idSample}.recal.bam")
        tuple val(idPatient), val(idSample)

    script:
    """
    samtools index ${idSample}.recal.bam
    """
}

def writeTsvFilesForRecalibrationReports(
        recalTableTSV,\
        genderMap,\
        statusMap) {

    // Create TSV files to restart from this step
    if (params.skip_markduplicates) {
        recalTableTSV.map { idPatient, idSample ->
            status = statusMap[idPatient, idSample]
            gender = genderMap[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
            recalTable = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.recal.table"
            "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"
        }.collectFile(
            name: 'mapped_no_duplicates_marked.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
        )

        recalTableTSV
            .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV/") {
                idPatient, idSample ->
                status = statusMap[idPatient, idSample]
                gender = genderMap[idPatient]
                bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
                bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
                recalTable = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.recal.table"
                ["mapped_no_duplicates_marked_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"]
        }
    } else {
        recalTableTSV.map { idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
        recalTable = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.recal.table"

            "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"
        }.collectFile(
            name: 'duplicates_marked.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
        )

        recalTableTSV
            .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV/") {
                idPatient, idSample ->
                status = statusMap[idPatient, idSample]
                gender = genderMap[idPatient]
                bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
                bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
                recalTable = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.recal.table"
                ["duplicates_marked_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${recalTable}\n"]
        }
    }
}

def writeTsvFilesForRecalibratedBams(\
        tsv_bam_recalibrated,\
        genderMap,
        statusMap) {

    // Creating a TSV file to restart from this step
    tsv_bam_recalibrated.map { idPatient, idSample ->
        gender = genderMap[idPatient]
        status = statusMap[idPatient, idSample]
        bam = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam.bai"
        "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
    }.collectFile(
        name: 'recalibrated.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
    )

    tsv_bam_recalibrated
        .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") {
            idPatient, idSample ->
            status = statusMap[idPatient, idSample]
            gender = genderMap[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam.bai"
            ["recalibrated_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
    }
}