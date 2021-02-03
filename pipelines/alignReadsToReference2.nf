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

    extractInfos
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
    CreateIntervalBeds
} from "${params.modulesDir}/preprocess.nf"

printHelpMessageAndExitIfUserAsks()

step = getInputStep()
tools = getInputTools(step)
skipQC = getInputSkipQC()
annotate_tools = getInputListOfAnnotationTools()
custom_runName = getCustomRunName()
save_bam_mapped = getSavedBamMapped()
tsvPath = getInputTsvPath()

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

    CreateIntervalBeds(ch_intervals_update).set { bedIntervals }

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



