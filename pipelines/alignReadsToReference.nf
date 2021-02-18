#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    create_workflow_summary;
    printNfcoreSarekWelcomeGraphic;
    checkHostname
} from "${params.modulesDir}/nfcore.nf"

include {
    printHelpMessageAndExitIfUserAsks;
    getInputStep;
    getInputTools;
    getInputSkipQC;
    getInputListOfAnnotationTools;
    getCustomRunName;
    getSavedBamMapped;
    getInputTsvPath;
    getInputSampleListAsChannel;
    getMultiqcConfigFile;
    getMultiqcCustomConfigFileAsChannel;
    getOutputDocsFile;
    getOutputDocsImagesFile;
    getSummaryMapFromParamsScopeAndArgs;

    initializeDerivedParams;
    
    checkInputReferenceGenomeExists;
    checkInputStepIsValid;
    checkInputToolsExist;
    checkInputSkippedQCToolsExist;
    checkInputListOfAnnotationToolsValid;
    checkInputAscatParametersValid;
    checkInputReadStructureParametersValid;
    checkAwsBatchSettings;
    checkInputTsvPath;

    printSummaryMessage;

    extractInfos;
    hasExtension
} from "${params.modulesDir}/sarek.nf"

include {
    BuildBWAindexes;
    BuildDict;
    BuildFastaFai;
    BuildDbsnpIndex;
    BuildGermlineResourceIndex;
    BuildKnownIndelsIndex;
    BuildPonIndex;
    BuildIntervals
} from "${params.modulesDir}/indices.nf"

include {
    CreateIntervalBeds;
    FastQCFQ;
    FastQCBAM;
    TrimGalore;
    UMIFastqToBAM;
    UMIMapBamFile;
    GroupReadsByUmi;
    CallMolecularConsensusReads
} from "${params.modulesDir}/preprocess.nf"

include {
    MapReads;
    MergeBamMapped;
    IndexBamFile;
    MarkDuplicates
} from "${params.modulesDir}/alignment.nf"

printHelpMessageAndExitIfUserAsks()

step = getInputStep()
tools = getInputTools(step)
skipQC = getInputSkipQC()
annotate_tools = getInputListOfAnnotationTools()
custom_runName = getCustomRunName()
save_bam_mapped = getSavedBamMapped()
tsvPath = getInputTsvPath(step)

initializeParamsScope(step, tools)

//derivedParams = initializeDerivedParams(step, tools)
//summaryMap = getSummaryMapFromParamsScopeAndArgs(step, custom_runName, skipQC, tools)
// getSummaryMapFrom...() needs to be called after initializeParamsScope()

checkHostname()
checkInputReferenceGenomeExists()
checkInputStepIsValid(step)
checkInputToolsExist(tools)
checkInputSkippedQCToolsExist(skipQC)
checkInputListOfAnnotationToolsValid(annotate_tools)
checkInputAscatParametersValid()
checkInputReadStructureParametersValid()
checkAwsBatchSettings()
checkInputTsvPath(tsvPath)

//printNfcoreSarekWelcomeGraphic()
//printSummaryMessage(summaryMap)
//printMutec2Warning(tools)

ch_multiqc_config = getMultiqcConfigFile()
ch_multiqc_custom_config = getMultiqcCustomConfigFileAsChannel()
ch_output_docs = getOutputDocsFile()
ch_output_docs_images = getOutputDocsImagesFile()

inputSample = getInputSampleListAsChannel(tsvPath, step)

(genderMap, statusMap, inputSample) = extractInfos(inputSample)


workflow {

    /*  Get indexes as channels  */

    ch_fasta = params.fasta && !('annotate' in step) ? Channel.value(file(params.fasta)) : "null"
    ch_dbsnp = params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || params.sentieon) ? Channel.value(file(params.dbsnp)) : "null"
    ch_germline_resource = params.germline_resource && 'mutect2' in tools ? Channel.value(file(params.germline_resource)) : "null"
    ch_known_indels = params.known_indels && ('mapping' in step || 'preparerecalibration' in step) ? Channel.value(file(params.known_indels)) : "null"
    ch_pon = params.pon ? Channel.value(file(params.pon)) : "null"
    ch_fai = params.fasta_fai && !('annotate' in step) ? Channel.value(file(params.fasta_fai)) : "null"
    ch_intervals = params.intervals && !params.no_intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"

    // Optional values, not defined within the params.genomes[params.genome] scope
    ch_read_structure1 = params.read_structure1 ? Channel.value(params.read_structure1) : "null"
    ch_read_structure2 = params.read_structure2 ? Channel.value(params.read_structure2) : "null"

    BuildBWAindexes(ch_fasta).set { bwa_built }
    ch_bwa = params.bwa ? Channel.value(file(params.bwa)) : bwa_built

    BuildDict(ch_fasta).set { dictBuilt }
    ch_dict = params.dict ? Channel.value(file(params.dict)) : dictBuilt

    BuildFastaFai(ch_fasta).set { fai_built }
    ch_fai_update = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : fai_built

    BuildDbsnpIndex(ch_dbsnp).set { dbsnp_tbi }
    ch_dbsnp_tbi = params.dbsnp ? params.dbsnp_index ? Channel.value(file(params.dbsnp_index)) : dbsnp_tbi : "null"

    BuildGermlineResourceIndex(ch_germline_resource).set { germline_resource_tbi }
    ch_germline_resource_tbi = params.germline_resource ? params.germline_resource_index ? Channel.value(file(params.germline_resource_index)) : germline_resource_tbi : "null"

    BuildKnownIndelsIndex(ch_known_indels).set { known_indels_tbi }
    ch_known_indels_tbi = params.known_indels ? params.known_indels_index ? Channel.value(file(params.known_indels_index)) : known_indels_tbi.collect() : "null"

    BuildPonIndex(ch_pon).set { pon_tbi }
    ch_pon_tbi = params.pon ? params.pon_index ? Channel.value(file(params.pon_index)) : pon_tbi : "null"

    BuildIntervals(ch_fai_update).set { intervalBuilt }
    ch_intervals_update = params.no_intervals ? "null" : (params.intervals && !('annotate' in step)) ? Channel.value(file(params.intervals)) : intervalBuilt

    /*  preprocess reads  */

    CreateIntervalBeds(ch_intervals_update).flatten().set { bedIntervals }

    bedIntervals = bedIntervals
        .map { intervalFile ->
            def duration = 0.0
            for (line in intervalFile.readLines()) {
                final fields = line.split('\t')
                if (fields.size() >= 5) duration += fields[4].toFloat()
                else {
                    start = fields[1].toInteger()
                    end = fields[2].toInteger()
                    duration += (end - start) / params.nucleotides_per_second
                }
            }
            return [duration, intervalFile]
        }
        .toSortedList({ a, b -> b[0] <=> a[0] })
        .flatten().collate(2)
        .map{duration, intervalFile -> intervalFile}
        .dump(tag:'bedintervals')

    if (params.no_intervals && step != 'annotate') {
        file("${params.outdir}/no_intervals.bed").text = "no_intervals\n"
        bedIntervals = Channel.from(file("${params.outdir}/no_intervals.bed"))
    }

    intBaseRecalibrator = bedIntervals
    intApplyBQSR = bedIntervals
    intHaplotypeCaller = bedIntervals
    intFreebayesSingle = bedIntervals
    intMpileup = bedIntervals

    inputSample.branch {
        bam: hasExtension(it[3], "bam")
        pairReads: !hasExtension(it[3], "bam")
    }
    .set { input }
    inputBam = input.bam
    inputPairReads = input.pairReads

    inputBam.dump(tag: "inputBam")
    inputPairReads.dump(tag: "inputPairReads")

    inputBamFastQC = inputBam

    // Removing inputFile2 which is null in case of uBAM
    inputBamFastQC = inputBamFastQC.map {
        idPatient, idSample, idRun, inputFile1, inputFile2 ->
        [idPatient, idSample, idRun, inputFile1]
    }

    if (params.split_fastq){
        inputPairReads = inputPairReads
            // newly splitfastq are named based on split, so the name is easier to catch
            .splitFastq(by: params.split_fastq, compress:true, file:"split", pe:true)
            .map {idPatient, idSample, idRun, reads1, reads2 ->
                // The split fastq read1 is the 4th element (indexed 3) its name is split_3
                // The split fastq read2's name is split_4
                // It's followed by which split it's acutally based on the mother fastq file
                // Index start at 1
                // Extracting the index to get a new IdRun
                splitIndex = reads1.fileName.toString().minus("split_3.").minus(".gz")
                newIdRun = idRun + "_" + splitIndex
                // Giving the files a new nice name
                newReads1 = file("${idSample}_${newIdRun}_R1.fastq.gz")
                newReads2 = file("${idSample}_${newIdRun}_R2.fastq.gz")
                [idPatient, idSample, newIdRun, reads1, reads2]
            }
    }

    inputPairReads.dump(tag:'INPUT')

    inputPairReadsTrimGalore = inputPairReads
    inputPairReadsFastQC = inputPairReads
    inputPairReadsUMI = inputPairReads

    FastQCFQ(inputPairReadsFastQC).set { fastQCFQReport }

    FastQCBAM(inputBamFastQC).set { fastQCBAMReport }

    fastQCReport = fastQCFQReport.mix(fastQCBAMReport)

    fastQCReport.dump(tag:'FastQC')

    (trimGaloreReport, outputPairReadsTrimGalore) = TrimGalore(inputPairReadsTrimGalore)

    /*  UMIs processing  */

    umi_converted_bams_ch = UMIFastqToBAM(inputPairReadsUMI, ch_read_structure1, ch_read_structure2)

    umi_aligned_bams_ch = UMIMapBamFile(umi_converted_bams_ch, ch_bwa, ch_fasta, ch_fai)

    (umi_histogram_ch, umi_grouped_bams_ch) = GroupReadsByUmi(umi_aligned_bams_ch)

    consensus_bam_ch = CallMolecularConsensusReads(umi_grouped_bams_ch)

    /*  map reads to reference  */

    input_pair_reads_sentieon = Channel.empty()

    if (params.umi) {
        inputPairReads = inputPairReads.dump(tag:'INPUT before mapping')
        if (params.sentieon) input_pair_reads_sentieon = consensus_bam_ch
        else inputPairReads = consensus_bam_ch
    }
    else {
        if (params.trim_fastq) inputPairReads = outputPairReadsTrimGalore
        else inputPairReads = inputPairReads.mix(inputBam)
        inputPairReads = inputPairReads.dump(tag:'INPUT before mapping')

        input_pair_reads_sentieon = inputPairReads
    }

    (bamMapped, bamMappedBamQC) = MapReads(inputPairReads, ch_bwa, ch_fasta, ch_fai)

    bamMapped.dump(tag:'Mapped BAM')

    bamMapped
        .groupTuple(by:[0, 1])
        .branch {
            single: !(it[2].size() > 1)
            multiple: it[2].size() > 1
        }
        .set { result }

    singleBam = result.single.dump(tag: "single")
    multipleBam = result.multiple.dump(tag: "multiple")
    
    singleBam = singleBam.map {
        idPatient, idSample, idRun, bam ->
        [idPatient, idSample, bam]
    }

    singleBam.dump(tag:'Single BAM')

    bam_mapped_merged = MergeBamMapped(multipleBam)

    //bam_mapped_merged.dump(tag:'Merged BAM')

    bam_mapped_merged = bam_mapped_merged.mix(singleBam)

    bam_sentieon_mapped_merged = bam_mapped_merged

    bam_mapped_merged.dump(tag:'BAMs for MD')

    bam_mapped_merged_to_index = bam_mapped_merged

    (bam_mapped_merged_indexed, tsv_bam_indexed) = IndexBamFile(bam_mapped_merged_to_index)

    tsv_bam_indexed_sample = tsv_bam_indexed

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

    tsv_bam_indexed_sample
        .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
            status = statusMap[idPatient, idSample]
            gender = genderMap[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/Mapped/${idSample}.bam.bai"
            ["mapped_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
    }

    /*  mark duplicates  */

    (bam_duplicates_marked, tsv_bam_duplicates_marked, duplicates_marked_report) = MarkDuplicates(bam_mapped_merged)

    tsv_bam_duplicates_marked_sample = tsv_bam_duplicates_marked

    // Creating a TSV file to restart from this step
    tsv_bam_duplicates_marked.map { idPatient, idSample ->
        gender = genderMap[idPatient]
        status = statusMap[idPatient, idSample]
        bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
        "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
    }.collectFile(
        name: 'duplicates_marked_no_table.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
    )

    tsv_bam_duplicates_marked_sample
        .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
            status = statusMap[idPatient, idSample]
            gender = genderMap[idPatient]
            bam = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam"
            bai = "${params.outdir}/Preprocessing/${idSample}/DuplicatesMarked/${idSample}.md.bam.bai"
            ["duplicates_marked_no_table_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
    }


    bam_duplicates_marked.dump(tag:'MD BAM')
   
    duplicates_marked_report.dump(tag:'MD Report')

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



