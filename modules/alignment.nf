include {
    getInputSkipQC
    hasExtension;
} from "${params.modulesDir}/sarek.nf"

process AlignReadsToReferenceSequence {
    label 'cpus_max'

    tag "${idPatient}-${idRun}"

    input:
        tuple val(idPatient), val(idSample), val(idRun), file(inputFile1), file(inputFile2)
        file(bwaIndex)
        file(fasta)
        file(fastaFai) 

    output:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_${idRun}.bam")

    script:
    // -K is an hidden option, used to fix the number of reads processed by bwa mem
    // Chunk size can affect bwa results, if not specified,
    // the number of threads can change which can give not deterministic result.
    // cf https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md
    // and https://github.com/gatk-workflows/gatk4-data-processing/blob/8ffa26ff4580df4ac3a5aa9e272a4ff6bab44ba2/processing-for-variant-discovery-gatk4.b37.wgs.inputs.json#L29
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
    convertToFastq = hasExtension(inputFile1, "bam") ? "gatk --java-options -Xmx${task.memory.toGiga()}g SamToFastq --INPUT=${inputFile1} --FASTQ=/dev/stdout --INTERLEAVE=true --NON_PF=true | \\" : ""
    input = hasExtension(inputFile1, "bam") ? "-p /dev/stdin - 2> >(tee ${inputFile1}.bwa.stderr.log >&2)" : "${inputFile1} ${inputFile2}"
    aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
    sortMemory = "${params.single_cpu_mem}".replaceAll(~/\s/,"").replaceAll(~/"GB"/,"G")
    """
    ${convertToFastq}
    ${aligner} mem -K 100000000 -R \"${readGroup}\" -t ${task.cpus} -M ${fasta} \
    ${input} | \
    samtools sort --threads ${task.cpus} -m ${sortMemory} - > ${idSample}_${idRun}.bam
    """
}

process MergeReadGroupsForSample {
    label 'cpus_8'

    tag "${idPatient}-${idSample}"

    input:
        tuple val(idPatient), val(idSample), val(idRun), file(bam)

    output:
        tuple val(idPatient), val(idSample), file("${idSample}.bam")

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.bam ${bam}
    """
}

process GetIndexOfAlignedSampleReadGroup {
    label 'cpus_8'

    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (save_bam_mapped) "Preprocessing/${idSample}/Mapped/${it}"
            else null
        }

    input:
        tuple val(idPatient), val(idSample), file("${idSample}.bam")

    output:
        tuple val(idPatient), val(idSample), file("${idSample}.bam.bai")
        tuple val(idPatient), val(idSample)

    when: save_bam_mapped

    script:
    """
    samtools index ${idSample}.bam
    """
}

process MarkDuplicatesInSampleReadGroup {
    label 'cpus_16'

    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${idSample}.bam.metrics") "Reports/${idSample}/MarkDuplicates/${it}"
            else "Preprocessing/${idSample}/DuplicatesMarked/${it}"
        }

    input:
        tuple val(idPatient), val(idSample), file("${idSample}.bam")

    output:
        tuple val(idPatient), val(idSample), file("${idSample}.md.bam"), file("${idSample}.md.bam.bai")
        tuple val(idPatient), val(idSample)
        file ("${idSample}.bam.metrics") optional true

    when: !(params.skip_markduplicates)

    script:
    //markdup_java_options = task.memory.toGiga() > 8 ? params.markdup_java_options : "\"-Xms" +  (task.memory.toGiga() / 2).trunc() + "g -Xmx" + (task.memory.toGiga() - 1) + "g\""
    markdup_java_options = params.markdup_java_options
    metrics = 'markduplicates' in getInputSkipQC() ? '' : "-M ${idSample}.bam.metrics"
    if (params.use_gatk_spark)
    """
    gatk --java-options ${markdup_java_options} \
        MarkDuplicatesSpark \
        -I ${idSample}.bam \
        -O ${idSample}.md.bam \
        ${metrics} \
        --tmp-dir . \
        --create-output-bam-index true \
        --spark-master local[${task.cpus}]
    """
    else
    """
    gatk --java-options ${markdup_java_options} \
        MarkDuplicates \
        --INPUT ${idSample}.bam \
        --METRICS_FILE ${idSample}.bam.metrics \
        --TMP_DIR . \
        --ASSUME_SORT_ORDER coordinate \
        --CREATE_INDEX true \
        --OUTPUT ${idSample}.md.bam
    
    mv ${idSample}.md.bai ${idSample}.md.bam.bai
    """
}

def branchIntoSingleOrMultipleGroupChannels(groupsOfReadGroups) {
    result = groupsOfReadGroups
        .branch {
            single: !(it[2].size() > 1)
            multiple: it[2].size() > 1
        }

    groupsWithASingleReadGroup = result.single
    groupsWithMulitpleReadGroups = result.multiple
    
    // the run id distinguishes different read groups for the 
    // same sample id. For groups containing a single read group, 
    // the run id is no longer needed so we remove it. 
    groupsWithASingleReadGroup = groupsWithASingleReadGroup.map {
        idPatient, idSample, idRun, bam ->
        [idPatient, idSample, bam]
    }

    return [groupsWithASingleReadGroup, groupsWithMulitpleReadGroups]   
}

def writeTsvFilesForBams(tsv_bam_indexed, genderMap, statusMap) {

    // Creating a TSV file to restart from this step
    tsv_bam_indexed.map { idPatient, idSample ->
        gender = genderMap[idPatient]
        status = statusMap[idPatient, idSample]
        bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
        "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
    }.collectFile(
        name: 'mapped.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
    )

    tsv_bam_indexed
        .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
            status = statusMap[idPatient, idSample]
            gender = genderMap[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
            ["mapped_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
    }
}

def writeTsvFilesForBamsWithDuplicatesMarked(tsv_bam_duplicates_marked, genderMap, statusMap) {
    tsv_bam_duplicates_marked
        .map { idPatient, idSample ->
            gender = genderMap[idPatient]
            status = statusMap[idPatient, idSample]
            bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
            "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
        }.collectFile(
            name: 'duplicates_marked_no_table.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
        )

    tsv_bam_duplicates_marked
        .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
            status = statusMap[idPatient, idSample]
            gender = genderMap[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
            ["duplicates_marked_no_table_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
    }
}
