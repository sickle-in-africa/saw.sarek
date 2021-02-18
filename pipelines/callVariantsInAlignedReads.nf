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
    printMutec2Warning;

    extractInfos;
    hasExtension
} from "${params.modulesDir}/sarek.nf"

include {
    get_software_versions;
    HaplotypeCaller;
    GenotypeGVCFs;
    StrelkaSingle;
    MantaSingle;
    TIDDIT;
    FreebayesSingle;
    ConcatVCF
} from "${params.modulesDir}/calling.nf"

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

step = 'variantcalling'
tools = getInputTools(step)
skipQC = getInputSkipQC()
annotate_tools = getInputListOfAnnotationTools()
custom_runName = getCustomRunName()
save_bam_mapped = getSavedBamMapped()
tsvPath = getInputTsvPath(step)

initializeParamsObject(step, tools)

//summaryMap = getSummaryMapFromParamsObjectAndArgs(step, custom_runName, skipQC, tools)

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

    // Initialize channels with files based on params
    ch_ac_loci = params.ac_loci && 'ascat' in tools ? Channel.value(file(params.ac_loci)) : "null"
    ch_ac_loci_gc = params.ac_loci_gc && 'ascat' in tools ? Channel.value(file(params.ac_loci_gc)) : "null"
    ch_chr_dir = params.chr_dir && 'controlfreec' in tools ? Channel.value(file(params.chr_dir)) : "null"
    ch_chr_length = params.chr_length && 'controlfreec' in tools ? Channel.value(file(params.chr_length)) : "null"
    ch_dbsnp = params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || params.sentieon) ? Channel.value(file(params.dbsnp)) : "null"
    ch_fasta = params.fasta && !('annotate' in step) ? Channel.value(file(params.fasta)) : "null"
    ch_fai = params.fasta_fai && !('annotate' in step) ? Channel.value(file(params.fasta_fai)) : "null"
    ch_germline_resource = params.germline_resource && 'mutect2' in tools ? Channel.value(file(params.germline_resource)) : "null"
    ch_intervals = params.intervals && !params.no_intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : "null"
    ch_known_indels = params.known_indels && ('mapping' in step || 'preparerecalibration' in step) ? Channel.value(file(params.known_indels)) : "null"
    ch_mappability = params.mappability && 'controlfreec' in tools ? Channel.value(file(params.mappability)) : "null"

    // Initialize channels with values based on params
    ch_snpeff_cache = params.snpeff_cache ? Channel.value(file(params.snpeff_cache)) : "null"
    ch_snpeff_db = params.snpeff_db ? Channel.value(params.snpeff_db) : "null"
    ch_vep_cache_version = params.vep_cache_version ? Channel.value(params.vep_cache_version) : "null"
    ch_vep_cache = params.vep_cache ? Channel.value(file(params.vep_cache)) : "null"

    // Optional files, not defined within the params.genomes[params.genome] scope
    ch_cadd_indels = params.cadd_indels ? Channel.value(file(params.cadd_indels)) : "null"
    ch_cadd_indels_tbi = params.cadd_indels_tbi ? Channel.value(file(params.cadd_indels_tbi)) : "null"
    ch_cadd_wg_snvs = params.cadd_wg_snvs ? Channel.value(file(params.cadd_wg_snvs)) : "null"
    ch_cadd_wg_snvs_tbi = params.cadd_wg_snvs_tbi ? Channel.value(file(params.cadd_wg_snvs_tbi)) : "null"
    ch_pon = params.pon ? Channel.value(file(params.pon)) : "null"
    ch_target_bed = params.target_bed ? Channel.value(file(params.target_bed)) : "null"

    // Optional values, not defined within the params.genomes[params.genome] scope
    ch_read_structure1 = params.read_structure1 ? Channel.value(params.read_structure1) : "null"
    ch_read_structure2 = params.read_structure2 ? Channel.value(params.read_structure2) : "null"

    ch_software_versions_yaml = get_software_versions()

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


    // When starting with variant calling, Channel bam_recalibrated is inputSample
    if (step == 'variantcalling') bam_recalibrated = inputSample

    bamMantaSingle = bam_recalibrated
    bamStrelkaSingle = bam_recalibrated
    bamTIDDIT = bam_recalibrated
    bamFreebayesSingleNoIntervals = bam_recalibrated
    bamHaplotypeCallerNoIntervals = bam_recalibrated
    bamRecalAll = bam_recalibrated

    // To speed Variant Callers up we are chopping the reference into smaller pieces
    // Do variant calling by this intervals, and re-merge the VCFs

    bamHaplotypeCaller = bamHaplotypeCallerNoIntervals.combine(intHaplotypeCaller)
    bamFreebayesSingle = bamFreebayesSingleNoIntervals.combine(intFreebayesSingle)

    (gvcfHaplotypeCaller,\
     gvcfGenotypeGVCFs)\
         = HaplotypeCaller(\
            bamHaplotypeCaller,\
            ch_dict,\
            ch_fasta,\
            ch_fai_update)

    gvcfHaplotypeCaller.groupTuple(by:[0, 1, 2]).set { gvcfHaplotypeCallerGrouped }

    gvcfHaplotypeCaller.dump(tag: 'GVCF HaplotypeCaller')

    (vcfGenotypeGVCFs)\
        = GenotypeGVCFs(\
            gvcfGenotypeGVCFs,\
            ch_dict,\
            ch_fasta,\
            ch_fai_update)

    vcfGenotypeGVCFs.groupTuple(by:[0, 1, 2]).set { vcfGenotypeGVCFsGrouped }

    (vcfStrelkaSingle)\
        = StrelkaSingle(\
            bamStrelkaSingle,\
            ch_fasta,\
            ch_fai_update,\
            ch_target_bed)
    
    vcfStrelkaSingle.dump(tag:'Strelka - Single Mode')

    (vcfMantaSingle)\
        = MantaSingle(\
            bamMantaSingle,\
            ch_fasta,\
            ch_fai_update,\
            ch_target_bed)

    vcfMantaSingle.dump(tag:'Single Manta')

    (vcfTIDDIT,\
     tidditOut)\
        = TIDDIT(\
            bamTIDDIT,\
            ch_fasta,\
            ch_fai_update)

    vcfTIDDIT.dump(tag:'TIDDIT')

    (vcfFreebayesSingle)\
        = FreebayesSingle(\
            bamFreebayesSingle,\
            ch_fasta,\
            ch_software_versions_yaml)

    vcfFreebayesSingle.groupTuple(by: [0,1,2]).set{ vcfFreebayesSingleGrouped }

    vcfConcatenateVCFs = Channel
        .empty()
        .mix(vcfFreebayesSingle, vcfGenotypeGVCFs, gvcfHaplotypeCaller)

    vcfConcatenateVCFs.dump(tag:'VCF to merge')

    (vcfConcatenated)\
        = ConcatVCF(\
            vcfConcatenateVCFs,\
            ch_fai_update,\
            ch_target_bed)

    vcfConcatenated.dump(tag:'VCF')
}



/*













// When starting with variant calling, Channel bam_recalibrated is inputSample
if (step == 'variantcalling') bam_recalibrated = inputSample

bam_recalibrated = bam_recalibrated.dump(tag:'BAM for Variant Calling')

// Here we have a recalibrated bam set
// The TSV file is formatted like: "idPatient status idSample bamFile baiFile"
// Manta will be run in Germline mode, or in Tumor mode depending on status
// HaplotypeCaller, TIDDIT and Strelka will be run for Normal and Tumor samples

(bamMantaSingle, bamStrelkaSingle, bamTIDDIT, bamFreebayesSingleNoIntervals, bamHaplotypeCallerNoIntervals, bamRecalAll) = bam_recalibrated.into(6)


// To speed Variant Callers up we are chopping the reference into smaller pieces
// Do variant calling by this intervals, and re-merge the VCFs

bamHaplotypeCaller = bamHaplotypeCallerNoIntervals.combine(intHaplotypeCaller)
bamFreebayesSingle = bamFreebayesSingleNoIntervals.combine(intFreebayesSingle)


// STEP GATK HAPLOTYPECALLER.1

process HaplotypeCaller {
    label 'memory_singleCPU_task_sq'
    label 'cpus_2'

    tag "${idSample}-${intervalBed.baseName}"

    input:
        set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamHaplotypeCaller
        //file(dbsnp) from ch_dbsnp
        //file(dbsnpIndex) from ch_dbsnp_tbi
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set val("HaplotypeCallerGVCF"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.g.vcf") into gvcfHaplotypeCaller
        set idPatient, idSample, file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf") into gvcfGenotypeGVCFs

    //when: 'haplotypecaller' in tools

    script:
    //javaOptions = "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
    javaOptions = ""
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    dbsnpOptions = "" //params.dbsnp ? "--D ${dbsnp}" : ""
    """
    gatk --java-options "${javaOptions}" \
        HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        -O ${intervalBed.baseName}_${idSample}.g.vcf \
        -ERC GVCF
    """
}

gvcfHaplotypeCaller = gvcfHaplotypeCaller.groupTuple(by:[0, 1, 2])

if (!params.generate_gvcf) gvcfHaplotypeCaller.close()
else gvcfHaplotypeCaller = gvcfHaplotypeCaller.dump(tag:'GVCF HaplotypeCaller')

// STEP GATK HAPLOTYPECALLER.2

process GenotypeGVCFs {
    tag "${idSample}-${intervalBed.baseName}"

    input:
        set idPatient, idSample, file(intervalBed), file(gvcf) from gvcfGenotypeGVCFs
        //file(dbsnp) from ch_dbsnp
        //file(dbsnpIndex) from ch_dbsnp_tbi
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
    set val("HaplotypeCaller"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.vcf") into vcfGenotypeGVCFs

    //when: 'haplotypecaller' in tools

    script:
    // Using -L is important for speed and we have to index the interval files also
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    dbsnpOptions = params.dbsnp ? "--D ${dbsnp}" : ""
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        IndexFeatureFile \
        -I ${gvcf}

    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GenotypeGVCFs \
        -R ${fasta} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        -V ${gvcf} \
        -O ${intervalBed.baseName}_${idSample}.vcf
    """
}

vcfGenotypeGVCFs = vcfGenotypeGVCFs.groupTuple(by:[0, 1, 2])



// STEP STRELKA.1 - SINGLE MODE

process StrelkaSingle {
    label 'cpus_max'
    label 'memory_max'

    tag "${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/Strelka", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file(bam), file(bai) from bamStrelkaSingle
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(targetBED) from ch_target_bed

    output:
        set val("Strelka"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfStrelkaSingle

    //when: 'strelka' in tools

    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    """
    ${beforeScript}
    configureStrelkaGermlineWorkflow.py \
        --bam ${bam} \
        --referenceFasta ${fasta} \
        ${options} \
        --runDir Strelka

    python Strelka/runWorkflow.py -m local -j ${task.cpus}

    mv Strelka/results/variants/genome.*.vcf.gz \
        Strelka_${idSample}_genome.vcf.gz
    mv Strelka/results/variants/genome.*.vcf.gz.tbi \
        Strelka_${idSample}_genome.vcf.gz.tbi
    mv Strelka/results/variants/variants.vcf.gz \
        Strelka_${idSample}_variants.vcf.gz
    mv Strelka/results/variants/variants.vcf.gz.tbi \
        Strelka_${idSample}_variants.vcf.gz.tbi
    """
}

vcfStrelkaSingle = vcfStrelkaSingle.dump(tag:'Strelka - Single Mode')

// STEP MANTA.1 - SINGLE MODE

process MantaSingle {
    label 'cpus_max'
    label 'memory_max'

    tag "${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/Manta", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file(bam), file(bai) from bamMantaSingle
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai
        file(targetBED) from ch_target_bed

    output:
        set val("Manta"), idPatient, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into vcfMantaSingle

    //when: 'manta' in tools

    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    status = statusMap[idPatient, idSample]
    inputbam = status == 0 ? "--bam" : "--tumorBam"
    vcftype = status == 0 ? "diploid" : "tumor"
    """
    ${beforeScript}
    configManta.py \
        ${inputbam} ${bam} \
        --reference ${fasta} \
        ${options} \
        --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}

    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${idSample}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${idSample}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${idSample}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/${vcftype}SV.vcf.gz \
        Manta_${idSample}.${vcftype}SV.vcf.gz
    mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
        Manta_${idSample}.${vcftype}SV.vcf.gz.tbi
    """
}

vcfMantaSingle = vcfMantaSingle.dump(tag:'Single Manta')

// STEP TIDDIT

process TIDDIT {
    tag "${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "TIDDIT_${idSample}.vcf") "VariantCalling/${idSample}/TIDDIT/${it}"
            else "Reports/${idSample}/TIDDIT/${it}"
        }

    input:
        set idPatient, idSample, file(bam), file(bai) from bamTIDDIT
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set val("TIDDIT"), idPatient, idSample, file("*.vcf.gz"), file("*.tbi") into vcfTIDDIT
        set file("TIDDIT_${idSample}.old.vcf"), file("TIDDIT_${idSample}.ploidy.tab"), file("TIDDIT_${idSample}.signals.tab"), file("TIDDIT_${idSample}.wig"), file("TIDDIT_${idSample}.gc.wig") into tidditOut

    //when: 'tiddit' in tools

    script:
    """
    tiddit --sv -o TIDDIT_${idSample} --bam ${bam} --ref ${fasta}

    mv TIDDIT_${idSample}.vcf TIDDIT_${idSample}.old.vcf

    grep -E "#|PASS" TIDDIT_${idSample}.old.vcf > TIDDIT_${idSample}.vcf

    bgzip --threads ${task.cpus} -c TIDDIT_${idSample}.vcf > TIDDIT_${idSample}.vcf.gz

    tabix TIDDIT_${idSample}.vcf.gz
    """
}

vcfTIDDIT = vcfTIDDIT.dump(tag:'TIDDIT')

// STEP FREEBAYES SINGLE MODE

process FreebayesSingle {
    tag "${idSample}-${intervalBed.baseName}"

    label 'cpus_1'
    
    input:
        set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamFreebayesSingle
        file(fasta) from ch_fasta
        file(fastaFai) from ch_software_versions_yaml
    
    output:
        set val("FreeBayes"), idPatient, idSample, file("${intervalBed.baseName}_${idSample}.vcf") into vcfFreebayesSingle
    
    //when: 'freebayes' in tools

    script:
    intervalsOptions = params.no_intervals ? "" : "-t ${intervalBed}"
    """
    freebayes \
        -f ${fasta} \
        --min-alternate-fraction 0.1 \
        --min-mapping-quality 1 \
        ${intervalsOptions} \
        ${bam} > ${intervalBed.baseName}_${idSample}.vcf
    """
}

vcfFreebayesSingle = vcfFreebayesSingle.groupTuple(by: [0,1,2])





















vcfConcatenateVCFs = Channel.empty().mix(vcfFreebayesSingle, vcfGenotypeGVCFs, gvcfHaplotypeCaller)

//vcfConcatenateVCFs = Channel.empty().mix(vcfGenotypeGVCFs, gvcfHaplotypeCaller)

vcfConcatenateVCFs = vcfConcatenateVCFs.dump(tag:'VCF to merge')

process ConcatVCF {
    label 'concat_vcf'
    label 'cpus_8'

    tag "${variantCaller}-${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/${"$variantCaller"}", mode: params.publish_dir_mode

    input:
        set variantCaller, idPatient, idSample, file(vcf) from vcfConcatenateVCFs
        file(fastaFai) from ch_fai
        file(targetBED) from ch_target_bed

    output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
        set variantCaller, idPatient, idSample, file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi") into vcfConcatenated

    //when: ('haplotypecaller' in tools || 'mutect2' in tools || 'freebayes' in tools)

    script:
    if (variantCaller == 'HaplotypeCallerGVCF')
        outputFile = "HaplotypeCaller_${idSample}.g.vcf"
    else
        outputFile = "${variantCaller}_${idSample}.vcf"
    options = params.target_bed ? "-t ${targetBED}" : ""
    intervalsOptions = params.no_intervals ? "-n" : ""
    """
    ${params.sarekDir}/bin/concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options} ${intervalsOptions}
    """
}

vcfConcatenated = vcfConcatenated.dump(tag:'VCF')


*/




def initializeParamsObject(inputStep, inputToolsList) {
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