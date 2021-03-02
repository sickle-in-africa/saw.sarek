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
    getInactiveChannel;

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
    checkInputTsvPath;
    checkTrimmingAndUmiFlags;
    checkBaseRecalibrationIsDesired;

    isChannelActive
} from "${params.modulesDir}/sarek.nf"

process GetSoftwareVersions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: {it.indexOf(".csv") > 0 ? it : null}

    output:
        file 'software_versions_mqc.yaml'
        //file "software_versions.csv"

    when: !('versions' in skipQC)

    script:
    aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
    """
    alleleCounter --version &> v_allelecount.txt 2>&1 || true
    bcftools --version &> v_bcftools.txt 2>&1 || true
    ${aligner} version &> v_bwa.txt 2>&1 || true
    cnvkit.py version &> v_cnvkit.txt 2>&1 || true
    configManta.py --version &> v_manta.txt 2>&1 || true
    configureStrelkaGermlineWorkflow.py --version &> v_strelka.txt 2>&1 || true
    echo "${workflow.manifest.version}" &> v_pipeline.txt 2>&1 || true
    echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true
    snpEff -version &> v_snpeff.txt 2>&1 || true
    fastqc --version &> v_fastqc.txt 2>&1 || true
    freebayes --version &> v_freebayes.txt 2>&1 || true
    freec &> v_controlfreec.txt 2>&1 || true
    gatk ApplyBQSR --help &> v_gatk.txt 2>&1 || true
    msisensor &> v_msisensor.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt 2>&1 || true
    qualimap --version &> v_qualimap.txt 2>&1 || true
    R --version &> v_r.txt 2>&1 || true
    R -e "library(ASCAT); help(package='ASCAT')" &> v_ascat.txt 2>&1 || true
    samtools --version &> v_samtools.txt 2>&1 || true
    tiddit &> v_tiddit.txt 2>&1 || true
    trim_galore -v &> v_trim_galore.txt 2>&1 || true
    vcftools --version &> v_vcftools.txt 2>&1 || true
    vep --help &> v_vep.txt 2>&1 || true

    ${params.sarekDir}/bin/scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

def initializeInputChannelsForMapping() {

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

    (__genderMap__, __statusMap__, inputSample) = extractInfos(inputSample)

    _fasta_ = params.fasta ? Channel.value(file(params.fasta)) : getInactiveChannel()
    _bwa_ = params.bwa ? Channel.value(file(params.bwa)) : getInactiveChannel()
    _fastaFai_ = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : getInactiveChannel()
    _intervalsList_ = params.intervals ? Channel.value(file(params.intervals)) : getInactiveChannel()

    return [\
        inputSample,\
        _fasta_,\
        _bwa_,\
        _fastaFai_,\
        __genderMap__,\
        __statusMap__]

}

def initializeInputChannelsForCalling() {

    step = 'variantcalling'
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

    checkTrimmingAndUmiFlags(params.trim_fastq, params.umi)

    inputSample = getInputSampleListAsChannel(tsvPath, step)

    (__genderMap__, __statusMap__, inputSample) = extractInfos(inputSample)

    _fasta_ = params.fasta ? Channel.value(file(params.fasta)) : getInactiveChannel()
    _dict_ = params.dict ? Channel.value(file(params.dict)) : getInactiveChannel()
    _fastaFai_ = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : getInactiveChannel()
    _dbsnp_ = params.dbsnp ? Channel.value(file(params.dbsnp)) : getInactiveChannel()
    _knownIndels_ = params.known_indels ? Channel.value(file(params.known_indels)) : getInactiveChannel()
    _intervalsList_ = params.intervals ? Channel.value(file(params.intervals)) : getInactiveChannel()

    _targetBed_ = params.target_bed ? Channel.value(file(params.target_bed)) : getInactiveChannel()

    return [\
        inputSample,\
        _fasta_,\
        _dict_,\
        _fastaFai_,\
        _dbsnp_,\
        _intervalsList_,\
        _targetBed_,\
        __genderMap__,\
        __statusMap__]

}

def initializeInputChannelsForRecalibration() {

    step = 'preparerecalibration'
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

    checkBaseRecalibrationIsDesired()

    // Initialize channels with files based on params
    _fasta_ = params.fasta ? Channel.value(file(params.fasta)) : getInactiveChannel()
    _dict_ = params.dict ? Channel.value(file(params.dict)) : getInactiveChannel()
    _fastaFai_ = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : getInactiveChannel()
    _dbsnp_ = params.dbsnp ? Channel.value(file(params.dbsnp)) : getInactiveChannel()
    _knownIndels_ = params.known_indels ? Channel.value(file(params.known_indels)) : getInactiveChannel()
    _intervalsList_ = params.intervals ? Channel.value(file(params.intervals)) : getInactiveChannel()
    _known_indels_ = params.known_indels ? Channel.value(file(params.known_indels)) : getInactiveChannel()

    inputSample = getInputSampleListAsChannel(tsvPath, step)

    (__genderMap__, __statusMap__, inputSample) = extractInfos(inputSample)

    return [\
        inputSample,\
        _fasta_,\
        _dict_,\
        _fastaFai_,\
        _dbsnp_,\
        _intervalsList_,\
        _known_indels_,\
        __genderMap__,\
        __statusMap__]

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