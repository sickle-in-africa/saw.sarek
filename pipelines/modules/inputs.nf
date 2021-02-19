include {
    create_workflow_summary;
    printNfcoreSarekWelcomeGraphic;
    checkHostname
} from "${params.modulesDir}/nfcore.nf"

include {
    getInputStep;
    getInputTools;
    getInputTsvPath;
    getInputSkipQC;
    getInputListOfAnnotationTools;
    getInputSampleListAsChannel;

    hasExtension;
    extractInfos;

    checkInputReferenceGenomeExists;
    checkInputStepIsValid;
    checkInputToolsExist;
    checkInputSkippedQCToolsExist;
    checkInputListOfAnnotationToolsValid;
    checkInputAscatParametersValid;
    checkInputReadStructureParametersValid;
    checkAwsBatchSettings;
    checkInputTsvPath
} from "${params.modulesDir}/sarek.nf"


def initializeInputChannelsForMapping(tools) {

    step = 'mapping'
    tools = getInputTools(step)
    skipQC = getInputSkipQC()
    annotate_tools = getInputListOfAnnotationTools()
    tsvPath = getInputTsvPath(step)

    initializeParamsScope(step, tools)

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

    inputSample = getInputSampleListAsChannel(tsvPath, step)

    (genderMap, statusMap, inputSample) = extractInfos(inputSample)

    _step_ = Channel.value(step)
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

    _bwa_ = params.bwa ? Channel.value(file(params.bwa)) : Channel.value(file('NULL'))
    //_bwa_ = Channel.value(file('NULL'))

    return [
        genderMap,\
        statusMap,\
        inputSample,\
        _step_,\
        ch_fasta,\
        ch_dbsnp,\
        ch_germline_resource,\
        ch_known_indels,\
        ch_pon,\
        ch_fai,\
        ch_intervals,
        ch_read_structure1,
        ch_read_structure2,
        _bwa_]

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