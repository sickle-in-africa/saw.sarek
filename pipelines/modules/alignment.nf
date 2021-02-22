include {
    getInputStep;
    getInputTools;
    getInputTsvPath;
    getInputSkipQC;
    hasExtension;
    getInputSampleListAsChannel;
    extractInfos
} from "${params.modulesDir}/sarek.nf"

step = getInputStep()
tools = getInputTools(step)
skipQC = getInputSkipQC()
tsvPath = getInputTsvPath(step)

initializeParamsScope(step, tools)

inputSample = getInputSampleListAsChannel(tsvPath, step)

(genderMap, statusMap, inputSample) = extractInfos(inputSample)

process MapReads {
    label 'cpus_max'

    tag "${idPatient}-${idRun}"

    input:
        tuple val(idPatient), val(idSample), val(idRun), file(inputFile1), file(inputFile2)
        file(bwaIndex)
        file(fasta) 
        file(fastaFai) 

    output:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_${idRun}.bam")
        tuple val(idPatient), val("${idSample}_${idRun}"), file("${idSample}_${idRun}.bam")

    script:
    // -K is an hidden option, used to fix the number of reads processed by bwa mem
    // Chunk size can affect bwa results, if not specified,
    // the number of threads can change which can give not deterministic result.
    // cf https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md
    // and https://github.com/gatk-workflows/gatk4-data-processing/blob/8ffa26ff4580df4ac3a5aa9e272a4ff6bab44ba2/processing-for-variant-discovery-gatk4.b37.wgs.inputs.json#L29
    CN = params.sequencing_center ? "CN:${params.sequencing_center}\\t" : ""
    readGroup = "@RG\\tID:${idRun}\\t${CN}PU:${idRun}\\tSM:${idSample}\\tLB:${idSample}\\tPL:illumina"
    // adjust mismatch penalty for tumor samples
    status = statusMap[idPatient, idSample]
    extra = status == 1 ? "-B 3" : ""
    convertToFastq = hasExtension(inputFile1, "bam") ? "gatk --java-options -Xmx${task.memory.toGiga()}g SamToFastq --INPUT=${inputFile1} --FASTQ=/dev/stdout --INTERLEAVE=true --NON_PF=true | \\" : ""
    input = hasExtension(inputFile1, "bam") ? "-p /dev/stdin - 2> >(tee ${inputFile1}.bwa.stderr.log >&2)" : "${inputFile1} ${inputFile2}"
    aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
    sortMemory = "${params.single_cpu_mem}".replaceAll(~/\s/,"").replaceAll(~/"GB"/,"G")
    """
    ${convertToFastq}
    ${aligner} mem -K 100000000 -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${fasta} \
    ${input} | \
    samtools sort --threads ${task.cpus} -m ${sortMemory} - > ${idSample}_${idRun}.bam
    """
}

process MergeBamMapped {
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

process IndexBamFile {
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
        tuple val(idPatient), val(idSample), file("${idSample}.bam"), file("${idSample}.bam.bai")
        tuple val(idPatient), val(idSample)

    when: save_bam_mapped || !(params.known_indels)

    script:
    """
    samtools index ${idSample}.bam
    """
}

process MarkDuplicates {
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
    metrics = 'markduplicates' in skipQC ? '' : "-M ${idSample}.bam.metrics"
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

def splitMappedBamsIntoSingleAndMultipeLanes(bamMapped) {

    result = bamMapped
        .groupTuple(by:[0, 1])
        .branch {
            single: !(it[2].size() > 1)
            multiple: it[2].size() > 1
        }

    singleBam = result.single
    multipleBam = result.multiple
    
    singleBam = singleBam.map {
        idPatient, idSample, idRun, bam ->
        [idPatient, idSample, bam]
    }

    return [singleBam, multipleBam]

}

def selectPairReadsChannelForMapping(\
    consensus_bam_ch,\
    outputPairReadsTrimGalore,\
    inputPairReadsSplit) {

    if ( params.umi ) {
        preprocessedPairReads = consensus_bam_ch
    }
    else {
        if ( params.trim_fastq ) {
            preprocessedPairReads = outputPairReadsTrimGalore
        } else {
            preprocessedPairReads = inputPairReadsSplit
        }
    }
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

def initializeParamsScope(inputStep, inputToolsList) {
  // Initialize each params in params.genomes, catch the command line first if it was defined
  // params.fasta has to be the first one
  params.fasta = params.genome && !('annotate' in inputStep) ? params.genomes[params.genome].fasta ?: null : null
  // The rest can be sorted
  params.ac_loci = params.genome && 'ascat' in inputToolsList ? params.genomes[params.genome].ac_loci ?: null : null
  params.ac_loci_gc = params.genome && 'ascat' in inputToolsList ? params.genomes[params.genome].ac_loci_gc ?: null : null
  params.bwa = params.genome && params.fasta && 'mapping' in inputStep ? params.genomes[params.genome].bwa ?: null : null
  params.chr_dir = params.genome && 'controlfreec' in inputToolsList ? params.genomes[params.genome].chr_dir ?: null : null
  params.chr_length = params.genome && 'controlfreec' in inputToolsList ? params.genomes[params.genome].chr_length ?: null : null
  params.dbsnp = params.genome && ('mapping' in inputStep || 'preparerecalibration' in inputStep || 'controlfreec' in inputToolsList || 'haplotypecaller' in inputToolsList || 'mutect2' in inputToolsList || params.sentieon) ? params.genomes[params.genome].dbsnp ?: null : null
  params.dbsnp_index = params.genome && params.dbsnp ? params.genomes[params.genome].dbsnp_index ?: null : null
  params.dict = params.genome && params.fasta ? params.genomes[params.genome].dict ?: null : null
  params.fasta_fai = params.genome && params.fasta ? params.genomes[params.genome].fasta_fai ?: null : null
  params.germline_resource = params.genome && 'mutect2' in inputToolsList ? params.genomes[params.genome].germline_resource ?: null : null
  params.germline_resource_index = params.genome && params.germline_resource ? params.genomes[params.genome].germline_resource_index ?: null : null
  params.intervals = params.genome && !('annotate' in inputStep) ? params.genomes[params.genome].intervals ?: null : null
  params.known_indels = params.genome && ('mapping' in inputStep || 'preparerecalibration' in inputStep) ? params.genomes[params.genome].known_indels ?: null : null
  params.known_indels_index = params.genome && params.known_indels ? params.genomes[params.genome].known_indels_index ?: null : null
  params.mappability = params.genome && 'controlfreec' in inputToolsList ? params.genomes[params.genome].mappability ?: null : null
  params.snpeff_db = params.genome && ('snpeff' in inputToolsList || 'merge' in inputToolsList) ? params.genomes[params.genome].snpeff_db ?: null : null
  params.species = params.genome && ('vep' in inputToolsList || 'merge' in inputToolsList) ? params.genomes[params.genome].species ?: null : null
  params.vep_cache_version = params.genome && ('vep' in inputToolsList || 'merge' in inputToolsList) ? params.genomes[params.genome].vep_cache_version ?: null : null
}