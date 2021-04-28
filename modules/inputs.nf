include {
    create_workflow_summary;
    printNfcoreSarekWelcomeGraphic;
    checkHostname
} from "${params.modulesDir}/nfcore.nf"

include {
    getInputTools;
    getInputTsvPath;
    getInputSkipQC;
    getInputListOfAnnotationTools;
    getInputSampleListAsChannel;
    getInactiveChannel;
    getInactiveValueChannel;

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
    checkBaseRecalibrationIsPossible;

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

    return [\
        inputSample,\
        _fasta_,\
        _bwa_,\
        _fastaFai_,\
        __genderMap__,\
        __statusMap__]

}

def initializeInputChannelsForRawReadQualityReporting() {

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

    ch_multiqc_config = Channel.value(file("${params.sarekDir}/assets/multiqc_config.yaml"))


    return [
        inputSample,
        _fasta_,
        _bwa_,
        _fastaFai_,
        __genderMap__,
        __statusMap__,
        ch_multiqc_config]

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

    _fasta_ = params.fasta ? Channel.value(file(params.fasta)) : getInactiveChannel('fasta')
    _dict_ = params.dict ? Channel.value(file(params.dict)) : getInactiveChannel('dict')
    _fastaFai_ = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : getInactiveChannel('fastaFai')
    _dbsnp_ = params.dbsnp ? Channel.value(file(params.dbsnp)) : getInactiveChannel('dbsnp')
    dbsnp_index = params.dbsnp_index ? Channel.value(file(params.dbsnp_index)) : getInactiveChannel('dbsnpIndex')
    _knownIndels_ = params.known_indels ? Channel.value(file(params.known_indels)) : getInactiveChannel('knownIndels')
    referenceIntervalList = params.intervals ? Channel.value(file(params.intervals)) : Channel.empty()

    _targetBed_ = params.target_bed ? Channel.value(file(params.target_bed)) : Channel.empty()

    return [
        inputSample,
        _fasta_,
        _dict_,
        _fastaFai_,
        _dbsnp_,
        dbsnp_index,
        referenceIntervalList,
        _targetBed_,
        __genderMap__,
        __statusMap__]

}

def initializeInputChannelsForPreparation() {

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
    _fasta_ = params.fasta ? Channel.value(file(params.fasta)) : getInactiveChannel('fasta')
    _dict_ = params.dict ? Channel.value(file(params.dict)) : getInactiveChannel('dict')
    _fastaFai_ = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : getInactiveChannel('fastaFai')
    _dbsnp_ = params.dbsnp ? Channel.value(file(params.dbsnp)) : getInactiveChannel('dbsnp')
    dbsnp_index = params.dbsnp_index ? Channel.value(file(params.dbsnp_index)) : getInactiveChannel('dbsnpIndex')
    referenceIntervalList = params.intervals ? Channel.value(file(params.intervals)) : Channel.empty()
    known_indels = params.known_indels ? Channel.value(file(params.known_indels)) : getInactiveChannel('knownIndels')
    known_indels_index = params.known_indels_index ? Channel.value(file(params.known_indels_index)) : getInactiveChannel('knownIndelsIndex')

    inputSample = getInputSampleListAsChannel(tsvPath, step)

    checkBaseRecalibrationIsPossible(params.dbsnp, params.known_indels)

    (__genderMap__, __statusMap__, inputSample) = extractInfos(inputSample)

    return [\
        inputSample,\
        _fasta_,\
        _dict_,\
        _fastaFai_,\
        _dbsnp_,\
        dbsnp_index,\
        referenceIntervalList,\
        known_indels,\
        known_indels_index,\
        __genderMap__,\
        __statusMap__]

}

def initializeInputChannelsForRecalibration() {

    step = 'recalibrate'
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

    // Initialize channels with files based on params
    _fasta_ = params.fasta ? Channel.value(file(params.fasta)) : getInactiveChannel('fasta')
    _dict_ = params.dict ? Channel.value(file(params.dict)) : getInactiveChannel('dict')
    _fastaFai_ = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : getInactiveChannel('fastaFai')
    _dbsnp_ = params.dbsnp ? Channel.value(file(params.dbsnp)) : getInactiveChannel('dbsnp')
    dbsnp_index = params.dbsnp_index ? Channel.value(file(params.dbsnp_index)) : getInactiveChannel('dbsnpIndex')
    referenceIntervalList = params.intervals ? Channel.value(file(params.intervals)) : Channel.empty()
    known_indels = params.known_indels ? Channel.value(file(params.known_indels)) : getInactiveChannel('knownIndels')
    known_indels_index = params.known_indels_index ? Channel.value(file(params.known_indels_index)) : getInactiveChannel('knownIndelsIndex')

    inputSample = getInputSampleListAsChannel(tsvPath, step)

    (__genderMap__, __statusMap__, inputSample) = extractInfos(inputSample)

    return [\
        inputSample,\
        _fasta_,\
        _dict_,\
        _fastaFai_,\
        _dbsnp_,\
        dbsnp_index,\
        referenceIntervalList,\
        known_indels,\
        known_indels_index,\
        __genderMap__,\
        __statusMap__]

}


def initializeInputChannelsForVQSR() {

    step = 'recalibratevariants'
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

    // Initialize channels with files based on params
    _fasta_ = params.fasta ? Channel.value(file(params.fasta)) : getInactiveChannel()
    _dict_ = params.dict ? Channel.value(file(params.dict)) : getInactiveChannel()
    _fastaFai_ = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : getInactiveChannel()
    _dbsnp_ = params.dbsnp ? Channel.value(file(params.dbsnp)) : getInactiveChannel()
    dbsnp_index = params.dbsnp_index ? Channel.value(file(params.dbsnp_index)) : getInactiveChannel()
    referenceIntervalList = params.intervals ? Channel.value(file(params.intervals)) : Channel.empty()
    known_indels = params.known_indels ? Channel.value(file(params.known_indels)) : getInactiveChannel()
    known_indels_index = params.known_indels_index ? Channel.value(file(params.known_indels_index)) : getInactiveChannel()

    inputSample = getInputSampleListAsChannel(tsvPath, step)

    (__genderMap__, __statusMap__, inputSample) = extractInfos(inputSample)

    return [\
        inputSample,\
        _fasta_,\
        _dict_,\
        _fastaFai_,\
        _dbsnp_,\
        dbsnp_index,\
        referenceIntervalList,\
        known_indels,\
        known_indels_index,\
        __genderMap__,\
        __statusMap__]

}

def initializeInputChannelsForAnnotation() {

    variantSets = Channel.empty().mix(
        Channel.fromPath("${params.outdir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
          .flatten().map{vcf -> ['HaplotypeCaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
          .flatten().map{vcf -> ['Manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/Mutect2/*.vcf.gz")
          .flatten().map{vcf -> ['Mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonDNAseq/*.vcf.gz")
          .flatten().map{vcf -> ['SentieonDNAseq', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonDNAscope/*.vcf.gz")
          .flatten().map{vcf -> ['SentieonDNAscope', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonTNscope/*.vcf.gz")
          .flatten().map{vcf -> ['SentieonTNscope', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/Strelka/*{somatic,variant}*.vcf.gz")
          .flatten().map{vcf -> ['Strelka', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/TIDDIT/*.vcf.gz")
          .flatten().map{vcf -> ['TIDDIT', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
        )

    params.snpeff_db = params.genomes[params.genome].snpeff_db
    snpeff_cache = params.genomes[params.genome].snpeff_cache

    vep_cache = params.genomes[params.genome].vep_cache
    params.vep_cache_version = params.genomes[params.genome].vep_cache_version

    ch_snpeff_config = Channel.value(file("${params.sarekDir}/conf/snpEff.config"))
    ch_snpeff_cache = snpeff_cache ? Channel.value(file(snpeff_cache)) : getInactiveChannel('snpeffCache')
    ch_snpeff_db = params.snpeff_db ? Channel.value(params.snpeff_db) : getInactiveValueChannel()
    ch_vep_cache = vep_cache ? Channel.value(file(vep_cache)) : getInactiveChannel('vepCache')
    ch_vep_cache_version = params.vep_cache_version ? Channel.value(params.vep_cache_version) : getInactiveValueChannel()

    ch_cadd_cache = params.cadd_cache ? Channel.value(file(params.cadd_cache)) : getInactiveChannel('caddCache')
    ch_cadd_indels = params.cadd_indels ? Channel.value(file(params.cadd_indels)) : getInactiveChannel('caddIndels')
    ch_cadd_indels_tbi = params.cadd_indels_tbi ? Channel.value(file(params.cadd_indels_tbi)) : getInactiveChannel('caddIndelsTbi')
    ch_cadd_wg_snvs = params.cadd_wg_snvs ? Channel.value(file(params.cadd_wg_snvs)) : getInactiveChannel('caddWgSnvs')
    ch_cadd_wg_snvs_tbi = params.cadd_wg_snvs_tbi ? Channel.value(file(params.cadd_wg_snvs_tbi)) : getInactiveChannel('caddWgSnvsTbi')

    return [
        variantSets,
        ch_snpeff_config,
        ch_snpeff_cache,
        ch_snpeff_db,
        ch_vep_cache,
        ch_vep_cache_version,
        ch_cadd_cache,
        ch_cadd_indels,
        ch_cadd_indels_tbi,
        ch_cadd_wg_snvs,
        ch_cadd_wg_snvs_tbi]
}

def initializeInputChannelsForVariantRecalibration() {
    
    variantSets = Channel.empty().mix(
        Channel.fromPath("${params.outdir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
            .flatten()
            .map{
                vcf -> 
                def variantCaller = 'HaplotypeCaller'
                def idSample = vcf.minus(vcf.fileName)[-2].toString()
                [variantCaller, idSample, vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/Strelka/*variants.vcf.gz")
            .flatten()
            .map{
                vcf -> 
                def variantCaller = 'Strelka'
                def idSample = vcf.minus(vcf.fileName)[-2].toString()
                [variantCaller, idSample, vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/FreeBayes/*.vcf.gz")
            .flatten()
            .map{
                vcf ->
                def variantCaller = 'FreeBayes'
                def idSample = vcf.minus(vcf.fileName)[-2].toString()
                [variantCaller, idSample, vcf]})

    _fasta_ = Channel.value(file(params.genomes[params.genome].fasta))
    _dict_ = Channel.value(file(params.genomes[params.genome].dict))
    _fastaFai_ = Channel.value(file(params.genomes[params.genome].fasta_fai))
    _dbsnp_ = Channel.value(file(params.genomes[params.genome].dbsnp))
    dbsnp_index = Channel.value(file(params.genomes[params.genome].dbsnp_index))
    hapmap = Channel.value(file(params.genomes[params.genome].hapmap))
    hapmap_index = Channel.value(file(params.genomes[params.genome].hapmap_index))
    onekg_snps = Channel.value(file(params.genomes[params.genome].onekg_snps))
    onekg_snps_index = Channel.value(file(params.genomes[params.genome].onekg_snps_index))
    onekg_indels = Channel.value(file(params.genomes[params.genome].onekg_indels))
    onekg_indels_index = Channel.value(file(params.genomes[params.genome].onekg_indels_index))
    onekg_omni = Channel.value(file(params.genomes[params.genome].onekg_omni))
    onekg_omni_index = Channel.value(file(params.genomes[params.genome].onekg_omni_index))
    axiom_exome_plus = Channel.value(file(params.genomes[params.genome].axiom_exome_plus))
    axiom_exome_plus_index = Channel.value(file(params.genomes[params.genome].axiom_exome_plus_index))


    return [
        variantSets,
        _fasta_,
        _dict_,
        _fastaFai_,
        _dbsnp_,
        dbsnp_index,
        hapmap,
        hapmap_index,
        onekg_snps,
        onekg_snps_index,
        onekg_indels,
        onekg_indels_index,
        onekg_omni,
        onekg_omni_index,
        axiom_exome_plus,
        axiom_exome_plus_index]

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
