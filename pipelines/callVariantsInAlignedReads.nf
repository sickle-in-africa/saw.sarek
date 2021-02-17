#!/usr/bin/env nextflow

printHelpMessageAndExitIfUserAsks()

step = 'variantcalling'
tools = getInputTools()
skipQC = getInputSkipQC()
annotate_tools = getInputListOfAnnotationTools()
custom_runName = getCustomRunName()
save_bam_mapped = getSavedBamMapped()
tsvPath = getInputTsvPath()

initializeParamsObject(step, tools)

summaryMap = getSummaryMapFromParamsObjectAndArgs(step, custom_runName, skipQC, tools)

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

printNfcoreSarekWelcomeGraphic()
printSummaryMessage(summaryMap)
printMutec2Warning(tools)

ch_multiqc_config = getMultiqcConfigFile()
ch_multiqc_custom_config = getMultiqcCustomConfigFileAsChannel()
ch_output_docs = getOutputDocsFile()
ch_output_docs_images = getOutputDocsImagesFile()

inputSample = getInputSampleListAsChannel(tsvPath)

(genderMap, statusMap, inputSample) = extractInfos(inputSample)


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


// Parse software version numbers

process get_software_versions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: {it.indexOf(".csv") > 0 ? it : null}

    output:
        file 'software_versions_mqc.yaml' into ch_software_versions_yaml
        file "software_versions.csv"

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

ch_software_versions_yaml = ch_software_versions_yaml.dump(tag:'SOFTWARE VERSIONS')






























/*
================================================================================
                                BUILDING INDEXES
================================================================================
*/

// And then initialize channels based on params or indexes that were just built

process BuildBWAindexes {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/BWAIndex/${it}" : null }

    input:
        file(fasta) from ch_fasta

    output:
        file("${fasta}.*") into bwa_built

    when: !(params.bwa) && params.fasta && 'mapping' in step

    script:
    aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
    """
    ${aligner} index ${fasta}
    """
}

ch_bwa = params.bwa ? Channel.value(file(params.bwa)) : bwa_built

process BuildDict {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(fasta) from ch_fasta

    output:
        file("${fasta.baseName}.dict") into dictBuilt

    when: !(params.dict) && params.fasta && !('annotate' in step) && !('controlfreec' in step)

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT ${fasta.baseName}.dict
    """
}

ch_dict = params.dict ? Channel.value(file(params.dict)) : dictBuilt

process BuildFastaFai {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(fasta) from ch_fasta

    output:
        file("${fasta}.fai") into fai_built

    when: !(params.fasta_fai) && params.fasta && !('annotate' in step)

    script:
    """
    samtools faidx ${fasta}
    """
}

ch_fai = params.fasta_fai ? Channel.value(file(params.fasta_fai)) : fai_built

process BuildDbsnpIndex {
    tag "${dbsnp}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(dbsnp) from ch_dbsnp

    output:
        file("${dbsnp}.tbi") into dbsnp_tbi

    when: !(params.dbsnp_index) && params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || 'tnscope' in tools)

    script:
    """
    tabix -p vcf ${dbsnp}
    """
}

ch_dbsnp_tbi = params.dbsnp ? params.dbsnp_index ? Channel.value(file(params.dbsnp_index)) : dbsnp_tbi : "null"

process BuildGermlineResourceIndex {
    tag "${germlineResource}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(germlineResource) from ch_germline_resource

    output:
        file("${germlineResource}.tbi") into germline_resource_tbi

    when: !(params.germline_resource_index) && params.germline_resource && 'mutect2' in tools

    script:
    """
    tabix -p vcf ${germlineResource}
    """
}

ch_germline_resource_tbi = params.germline_resource ? params.germline_resource_index ? Channel.value(file(params.germline_resource_index)) : germline_resource_tbi : "null"

process BuildKnownIndelsIndex {
    tag "${knownIndels}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        each file(knownIndels) from ch_known_indels

    output:
        file("${knownIndels}.tbi") into known_indels_tbi

    when: !(params.known_indels_index) && params.known_indels && ('mapping' in step || 'preparerecalibration' in step)

    script:
    """
    tabix -p vcf ${knownIndels}
    """
}

ch_known_indels_tbi = params.known_indels ? params.known_indels_index ? Channel.value(file(params.known_indels_index)) : known_indels_tbi.collect() : "null"

process BuildPonIndex {
    tag "${pon}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(pon) from ch_pon

    output:
        file("${pon}.tbi") into pon_tbi

    when: !(params.pon_index) && params.pon && ('tnscope' in tools || 'mutect2' in tools)

    script:
    """
    tabix -p vcf ${pon}
    """
}

ch_pon_tbi = params.pon ? params.pon_index ? Channel.value(file(params.pon_index)) : pon_tbi : "null"

process BuildIntervals {
    tag "${fastaFai}"

    publishDir params.outdir, mode: params.publish_dir_mode,
    saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(fastaFai) from ch_fai

    output:
        file("${fastaFai.baseName}.bed") into intervalBuilt

    when: !(params.intervals) && !('annotate' in step) && !('controlfreec' in step) 

    script:
    """
    awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fastaFai} > ${fastaFai.baseName}.bed
    """
}

ch_intervals = params.no_intervals ? "null" : params.intervals && !('annotate' in step) ? Channel.value(file(params.intervals)) : intervalBuilt




























/*
================================================================================
                                  PREPROCESSING
================================================================================
*/

// STEP 0: CREATING INTERVALS FOR PARALLELIZATION (PREPROCESSING AND VARIANT CALLING)

process CreateIntervalBeds {
    tag "${intervals}"

    input:
        file(intervals) from ch_intervals

    output:
        file '*.bed' into bedIntervals mode flatten

    when: (!params.no_intervals) && step != 'annotate'

    script:
    // If the interval file is BED format, the fifth column is interpreted to
    // contain runtime estimates, which is then used to combine short-running jobs
    if (hasExtension(intervals, "bed"))
        """
        awk -vFS="\t" '{
          t = \$5  # runtime estimate
          if (t == "") {
            # no runtime estimate in this row, assume default value
            t = (\$3 - \$2) / ${params.nucleotides_per_second}
          }
          if (name == "" || (chunk > 600 && (chunk + t) > longest * 1.05)) {
            # start a new chunk
            name = sprintf("%s_%d-%d.bed", \$1, \$2+1, \$3)
            chunk = 0
            longest = 0
          }
          if (t > longest)
            longest = t
          chunk += t
          print \$0 > name
        }' ${intervals}
        """
    else if (hasExtension(intervals, "interval_list"))
        """
        grep -v '^@' ${intervals} | awk -vFS="\t" '{
          name = sprintf("%s_%d-%d", \$1, \$2, \$3);
          printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }'
        """
    else
        """
        awk -vFS="[:-]" '{
          name = sprintf("%s_%d-%d", \$1, \$2, \$3);
          printf("%s\\t%d\\t%d\\n", \$1, \$2-1, \$3) > name ".bed"
        }' ${intervals}
        """
}

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
        [duration, intervalFile]
        }.toSortedList({ a, b -> b[0] <=> a[0] })
    .flatten().collate(2)
    .map{duration, intervalFile -> intervalFile}

bedIntervals = bedIntervals.dump(tag:'bedintervals')

if (params.no_intervals && step != 'annotate') {
    file("${params.outdir}/no_intervals.bed").text = "no_intervals\n"
    bedIntervals = Channel.from(file("${params.outdir}/no_intervals.bed"))
}

(intBaseRecalibrator, intApplyBQSR, intHaplotypeCaller, intFreebayesSingle, intMpileup, bedIntervals) = bedIntervals.into(6)


















/*
================================================================================
                            GERMLINE VARIANT CALLING
================================================================================
*/




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












/*
================================================================================
                                nf-core functions
================================================================================
*/

def create_workflow_summary(summary) {
    def yaml_file = workDir.resolve('workflow_summary_mqc.yaml')
    yaml_file.text  = """
    id: 'nf-core-sarek-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/sarek Workflow Summary'
    section_href: 'https://github.com/nf-core/sarek'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
${summary.collect { k, v -> "            <dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }.join("\n")}
        </dl>
    """.stripIndent()

   return yaml_file
}

def nfcoreHeader() {
    // Log colors ANSI codes
    c_black  = params.monochrome_logs ? '' : "\033[0;30m";
    c_blue   = params.monochrome_logs ? '' : "\033[0;34m";
    c_dim    = params.monochrome_logs ? '' : "\033[2m";
    c_green  = params.monochrome_logs ? '' : "\033[0;32m";
    c_purple = params.monochrome_logs ? '' : "\033[0;35m";
    c_reset  = params.monochrome_logs ? '' : "\033[0m";
    c_white  = params.monochrome_logs ? '' : "\033[0;37m";
    c_yellow = params.monochrome_logs ? '' : "\033[0;33m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                            ${c_green},--.${c_black}/${c_green},-.${c_reset}
    ${c_blue}        ___     __   __   __   ___     ${c_green}/,-._.--~\'${c_reset}
    ${c_blue}  |\\ | |__  __ /  ` /  \\ |__) |__         ${c_yellow}}  {${c_reset}
    ${c_blue}  | \\| |       \\__, \\__/ |  \\ |___     ${c_green}\\`-._,-`-,${c_reset}
                                            ${c_green}`._,._,\'${c_reset}
        ${c_white}____${c_reset}
      ${c_white}.´ _  `.${c_reset}
     ${c_white}/  ${c_green}|\\${c_reset}`-_ \\${c_reset}     ${c_blue} __        __   ___     ${c_reset}
    ${c_white}|   ${c_green}| \\${c_reset}  `-|${c_reset}    ${c_blue}|__`  /\\  |__) |__  |__/${c_reset}
     ${c_white}\\ ${c_green}|   \\${c_reset}  /${c_reset}     ${c_blue}.__| /¯¯\\ |  \\ |___ |  \\${c_reset}
      ${c_white}`${c_green}|${c_reset}____${c_green}\\${c_reset}´${c_reset}

    ${c_purple}  nf-core/sarek v${workflow.manifest.version}${c_reset}
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

def checkHostname() {
    def c_reset       = params.monochrome_logs ? '' : "\033[0m"
    def c_white       = params.monochrome_logs ? '' : "\033[0;37m"
    def c_red         = params.monochrome_logs ? '' : "\033[1;91m"
    def c_yellow_bold = params.monochrome_logs ? '' : "\033[1;93m"
    if (params.hostnames) {
        def hostname = "hostname".execute().text.trim()
        params.hostnames.each { prof, hnames ->
            hnames.each { hname ->
                if (hostname.contains(hname) && !workflow.profile.contains(prof)) {
                    log.error "====================================================\n" +
                            "  ${c_red}WARNING!${c_reset} You are running with `-profile $workflow.profile`\n" +
                            "  but your machine hostname is ${c_white}'$hostname'${c_reset}\n" +
                            "  ${c_yellow_bold}It's highly recommended that you use `-profile $prof${c_reset}`\n" +
                            "============================================================"
                }
            }
        }
    }
}

/*
================================================================================
                                 sarek functions
================================================================================
*/

// Check if a row has the expected number of item
def checkNumberOfItem(row, number) {
    if (row.size() != number) exit 1, "Malformed row in TSV file: ${row}, see --help for more information"
    return true
}

// Check parameter existence
def checkParameterExistence(it, list) {
    if (!list.contains(it)) {
        log.warn "Unknown parameter: ${it}"
        return false
    }
    return true
}

// Compare each parameter with a list of parameters
def checkParameterList(list, realList) {
    return list.every{ checkParameterExistence(it, realList) }
}

// Define list of available tools to annotate
def defineAnnoList() {
    return [
        'haplotypecaller',
        'manta',
        'mutect2',
        'strelka',
        'tiddit'
    ]
}

// Define list of skipable QC tools
def defineSkipQClist() {
    return [
        'bamqc',
        'baserecalibrator',
        'bcftools',
        'documentation',
        'fastqc',
        'markduplicates',
        'multiqc',
        'samtools',
        'sentieon',
        'vcftools',
        'versions'
    ]
}

// Define list of available step
def defineStepList() {
    return [
        'annotate',
        'controlfreec',
        'mapping',
        'preparerecalibration',
        'recalibrate',
        'variantcalling'
    ]
}

// Define list of available tools
def defineToolList() {
    return [
        'ascat',
        'cnvkit',
        'controlfreec',
        'dnascope',
        'dnaseq',
        'freebayes',
        'haplotypecaller',
        'manta',
        'merge',
        'mpileup',
        'mutect2',
        'snpeff',
        'strelka',
        'tiddit',
        'tnscope',
        'vep',
        'msisensor'
    ]
}

// Channeling the TSV file containing BAM.
// Format is: "subject gender status sample bam bai"
def extractBam(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 6)
            def idPatient = row[0]
            def gender    = row[1]
            def status    = returnStatus(row[2].toInteger())
            def idSample  = row[3]
            def bamFile   = returnFile(row[4])
            def baiFile   = returnFile(row[5])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"

            return [idPatient, gender, status, idSample, bamFile, baiFile]
        }
}

// Create a channel of germline FASTQs from a directory pattern: "my_samples/*/"
// All FASTQ files in subdirectories are collected and emitted;
// they must have _R1_ and _R2_ in their names.
def extractFastqFromDir(pattern) {
    def fastq = Channel.create()
    // a temporary channel does all the work
    Channel
        .fromPath(pattern, type: 'dir')
        .ifEmpty { error "No directories found matching pattern '${pattern}'" }
        .subscribe onNext: { sampleDir ->
            // the last name of the sampleDir is assumed to be a unique sample id
            sampleId = sampleDir.getFileName().toString()

            for (path1 in file("${sampleDir}/**_R1_*.fastq.gz")) {
                assert path1.getName().contains('_R1_')
                path2 = file(path1.toString().replace('_R1_', '_R2_'))
                if (!path2.exists()) error "Path '${path2}' not found"
                (flowcell, lane) = flowcellLaneFromFastq(path1)
                String random = org.apache.commons.lang.RandomStringUtils.random(8, true, true) // random string to avoid duplicate names
                patient = sampleId
                gender = 'ZZ'  // unused
                status = 0  // normal (not tumor)
                rgId = "${flowcell}.${sampleId}.${lane}.${random}"
                result = [patient, gender, status, sampleId, rgId, path1, path2]
                fastq.bind(result)
            }
    }, onComplete: { fastq.close() }
    fastq
}

// Extract gender and status from Channel
def extractInfos(channel) {
    def genderMap = [:]
    def statusMap = [:]
    channel = channel.map{ it ->
        def idPatient = it[0]
        def gender = it[1]
        def status = it[2]
        def idSample = it[3]
        genderMap[idPatient] = gender
        statusMap[idPatient, idSample] = status
        [idPatient] + it[3..-1]
    }
    [genderMap, statusMap, channel]
}

// Channeling the TSV file containing FASTQ or BAM
// Format is: "subject gender status sample lane fastq1 fastq2"
// or: "subject gender status sample lane bam"
def extractFastq(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = returnStatus(row[2].toInteger())
            def idSample   = row[3]
            def idRun      = row[4]
            def file1      = returnFile(row[5])
            def file2      = "null"
            if (hasExtension(file1, "fastq.gz") || hasExtension(file1, "fq.gz") || hasExtension(file1, "fastq") || hasExtension(file1, "fq")) {
                checkNumberOfItem(row, 7)
                file2 = returnFile(row[6])
                if (!hasExtension(file2, "fastq.gz") && !hasExtension(file2, "fq.gz")  && !hasExtension(file2, "fastq") && !hasExtension(file2, "fq")) exit 1, "File: ${file2} has the wrong extension. See --help for more information"
                if (hasExtension(file1, "fastq") || hasExtension(file1, "fq") || hasExtension(file2, "fastq") || hasExtension(file2, "fq")) {
                    exit 1, "We do recommend to use gziped fastq file to help you reduce your data footprint."
                }
            }
            else if (hasExtension(file1, "bam")) checkNumberOfItem(row, 6)
            else "No recognisable extention for input file: ${file1}"

            [idPatient, gender, status, idSample, idRun, file1, file2]
        }
}

// Channeling the TSV file containing mpileup
// Format is: "subject gender status sample pileup"
def extractPileup(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 5)
            def idPatient = row[0]
            def gender    = row[1]
            def status    = returnStatus(row[2].toInteger())
            def idSample  = row[3]
            def mpileup   = returnFile(row[4])

            if (!hasExtension(mpileup, "pileup")) exit 1, "File: ${mpileup} has the wrong extension. See --help for more information"

            return [idPatient, gender, status, idSample, mpileup]
        }
}

// Channeling the TSV file containing Recalibration Tables.
// Format is: "subject gender status sample bam bai recalTable"
def extractRecal(tsvFile) {
    Channel.from(tsvFile)
        .splitCsv(sep: '\t')
        .map { row ->
            checkNumberOfItem(row, 7)
            def idPatient  = row[0]
            def gender     = row[1]
            def status     = returnStatus(row[2].toInteger())
            def idSample   = row[3]
            def bamFile    = returnFile(row[4])
            def baiFile    = returnFile(row[5])
            def recalTable = returnFile(row[6])

            if (!hasExtension(bamFile, "bam")) exit 1, "File: ${bamFile} has the wrong extension. See --help for more information"
            if (!hasExtension(baiFile, "bai")) exit 1, "File: ${baiFile} has the wrong extension. See --help for more information"
            if (!hasExtension(recalTable, "recal.table")) exit 1, "File: ${recalTable} has the wrong extension. See --help for more information"

            [idPatient, gender, status, idSample, bamFile, baiFile, recalTable]
    }
}

// Parse first line of a FASTQ file, return the flowcell id and lane number.
def flowcellLaneFromFastq(path) {
    // expected format:
    // xx:yy:FLOWCELLID:LANE:... (seven fields)
    // or
    // FLOWCELLID:LANE:xx:... (five fields)
    InputStream fileStream = new FileInputStream(path.toFile())
    InputStream gzipStream = new java.util.zip.GZIPInputStream(fileStream)
    Reader decoder = new InputStreamReader(gzipStream, 'ASCII')
    BufferedReader buffered = new BufferedReader(decoder)
    def line = buffered.readLine()
    assert line.startsWith('@')
    line = line.substring(1)
    def fields = line.split(' ')[0].split(':')
    String fcid
    int lane
    if (fields.size() == 7) {
        // CASAVA 1.8+ format
        fcid = fields[2]
        lane = fields[3].toInteger()
    } else if (fields.size() == 5) {
        fcid = fields[0]
        lane = fields[1].toInteger()
    }
    [fcid, lane]
}

// Check file extension
def hasExtension(it, extension) {
    it.toString().toLowerCase().endsWith(extension.toLowerCase())
}

// Return file if it exists
def returnFile(it) {
    if (!file(it).exists()) exit 1, "Missing file in TSV file: ${it}, see --help for more information"
    return file(it)
}

// Remove .ann .gz and .vcf extension from a VCF file
def reduceVCF(file) {
    return file.fileName.toString().minus(".ann").minus(".vcf").minus(".gz")
}

// Return status [0,1]
// 0 == Normal, 1 == Tumor
def returnStatus(it) {
    if (!(it in [0, 1])) exit 1, "Status is not recognized in TSV file: ${it}, see --help for more information"
    return it
}


def getInputStep() {
    return params.step ? params.step.toLowerCase().replaceAll('-', '').replaceAll('_', '') : ''
}
def getInputTools() {
    if (step == 'controlfreec') {
        return ['controlfreec']
    } else {
      return params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []
    } 
}
def getInputSkipQC() {
  return params.skip_qc ? params.skip_qc == 'all' ? skipQClist : params.skip_qc.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []
}
def getInputListOfAnnotationTools() {
  return params.annotate_tools ? params.annotate_tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '')} : []
}
def getCustomRunName() {
  // Has the run name been specified by the user?
  // This has the bonus effect of catching both -name and --name
  if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) return workflow.runName
  else return params.name
}
def getMultiqcConfigFile() {
  return file("${params.sarekDir}/assets/multiqc_config.yaml", checkIfExists: true)
}
def getMultiqcCustomConfigFileAsChannel() {
  return params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
}
def getOutputDocsFile() {
  return file("${params.sarekDir}/docs/output.md", checkIfExists: true)
}
def getOutputDocsImagesFile(){
  return file("${params.sarekDir}/docs/images/", checkIfExists: true)
}
def getSavedBamMapped() {
  params.skip_markduplicates ? true : params.save_bam_mapped ? true : false
}
def getInputTsvPath() {

  if ( isInputTsvFileSpecified() ) {
    return getTsvPathFromUserInput()
  }
  else {
    return getTsvPathFromOutputOfPreviousStep()
    // only for steps:
    //    + preparerecalibration
    //    + recalibrate
    //    + variantcalling
    //    + controlfreec
  }
}
def getTsvPathFromUserInput() {
  return params.input
}
def getTsvPathFromOutputOfPreviousStep() {
  if (!params.skip_markduplicates) {
      switch (step) {
          case 'mapping': break
          case 'preparerecalibration': return "${params.outdir}/Preprocessing/TSV/duplicates_marked_no_table.tsv"; break
          case 'recalibrate': return "${params.outdir}/Preprocessing/TSV/duplicates_marked.tsv"; break
          case 'variantcalling':
            if ( params.skip_recalibration == false ) {
                return "${params.outdir}/Preprocessing/TSV/recalibrated.tsv"
            } else {
                return "${params.outdir}/Preprocessing/TSV/duplicates_marked_no_table.tsv"
            }
            break
          case 'controlfreec': return "${params.outdir}/VariantCalling/TSV/control-freec_mpileup.tsv"; break
          case 'annotate': break
          default: exit 1, "Unknown step ${step}"
      }
  } else if (params.skip_markduplicates) {
      switch (step) {
          case 'mapping': break
          case 'preparerecalibration': return "${params.outdir}/Preprocessing/TSV/mapped.tsv"; break
          case 'recalibrate': return "${params.outdir}/Preprocessing/TSV/mapped_no_duplicates_marked.tsv"; break
          case 'variantcalling': return "${params.outdir}/Preprocessing/TSV/recalibrated.tsv"; break
          case 'controlfreec': return"${params.outdir}/VariantCalling/TSV/control-freec_mpileup.tsv"; break
          case 'annotate': break
          default: exit 1, "Unknown step ${step}"
      }
  } 
}
def getInputSampleListAsChannel(inputTsvPath) {
  tsvFile = file(inputTsvPath)
  switch (step) {
      case 'mapping': return extractFastq(tsvFile)
      case 'preparerecalibration': return extractBam(tsvFile)
      case 'recalibrate': return extractRecal(tsvFile)
      case 'variantcalling': return extractBam(tsvFile)
      case 'controlfreec': return extractPileup(tsvFile)
      case 'annotate': break
      default: exit 1, "Unknown step ${step}"
  }
}

def checkInputReferenceGenomeExists() {
    if (params.genomes && !params.genomes.containsKey(params.genome) && !params.igenomes_ignore) {
        exit 1, "The provided genome '${params.genome}' is not available in the iGenomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
    } else if (params.genomes && !params.genomes.containsKey(params.genome) && params.igenomes_ignore) {
        exit 1, "The provided genome '${params.genome}' is not available in the genomes file. Currently the available genomes are ${params.genomes.keySet().join(", ")}"
    }
}
def checkInputStepIsValid(inputStep) {
    stepList = defineStepList()
    if (inputStep.contains(',')) exit 1, 'You can choose only one step, see --help for more information'
    if (!checkParameterExistence(inputStep, stepList)) exit 1, "Unknown step ${inputStep}, see --help for more information"
}
def checkInputToolsExist(inputTools) {
  toolList = defineToolList()
  if (!checkParameterList(inputTools, toolList)) exit 1, 'Unknown tool(s), see --help for more information'
}
def checkInputSkippedQCToolsExist(inputSkipQC) {
  skipQClist = defineSkipQClist()
  if (!checkParameterList(inputSkipQC, skipQClist)) exit 1, 'Unknown QC tool(s), see --help for more information'
}
def checkInputListOfAnnotationToolsValid(inputAnnotationTools) {
  annoList = defineAnnoList()
  if (!checkParameterList(annotate_tools,annoList)) exit 1, 'Unknown tool(s) to annotate, see --help for more information'
}
def checkInputAscatParametersValid() {
  if ((params.ascat_ploidy && !params.ascat_purity) || (!params.ascat_ploidy && params.ascat_purity)) exit 1, 'Please specify both --ascat_purity and --ascat_ploidy, or none of them'
}
def checkInputReadStructureParametersValid() {
  if (params.umi && !(params.read_structure1 && params.read_structure2)) exit 1, 'Please specify both --read_structure1 and --read_structure2, when using --umi'
}
def checkAwsBatchSettings() {
  if (workflow.profile.contains('awsbatch')) {
      // AWSBatch sanity checking
      if (!params.awsqueue || !params.awsregion) exit 1, "Specify correct --awsqueue and --awsregion parameters on AWSBatch!"
      // Check outdir paths to be S3 buckets if running on AWSBatch
      // related: https://github.com/nextflow-io/nextflow/issues/813
      if (!params.outdir.startsWith('s3:')) exit 1, "Outdir not on S3 - specify S3 Bucket to run on AWSBatch!"
      // Prevent trace files to be stored on S3 since S3 does not support rolling files.
      if (params.tracedir.startsWith('s3:')) exit 1, "Specify a local tracedir or run without trace! S3 cannot be used for tracefiles."
  }
}
def checkInputTsvPath(inputTsvPath) {
  if ( !hasExtension(inputTsvPath, "tsv") ) exit 1, "your input sample tsv file has the wrong extension. Please check and try again."
  if ( hasExtension(inputTsvPath, "vcf") || hasExtension(inputTsvPath, "vcf.gz") ) exit 1, "you have specified a vcf input instead of tsv. That was perhaps meant for the annotation step."
}

def isInputTsvFileSpecified() {
  if (params.input) return true
  else return false
}

def printHelpMessageAndExitIfUserAsks() {
  if (params.help == true) exit 0, helpMessage()
}

def helpMessage() {
    log.info nfcoreHeader()
    log.info"""

    Usage:

    The typical command for running the pipeline is as follows:

    nextflow run nf-core/sarek --input sample.tsv -profile docker

    Mandatory arguments:
      --input                      [file] Path to input TSV file on mapping, prepare_recalibration, recalibrate, variant_calling and Control-FREEC steps
                                          Multiple TSV files can be specified surrounded with quotes
                                          Works also with the path to a directory on mapping step with a single germline sample only
                                          Alternatively, path to VCF input file on annotate step
                                          Multiple VCF files can be specified surrounded with quotes
      -profile                      [str] Configuration profile to use. Can use multiple (comma separated)
                                          Available: conda, docker, singularity, test, awsbatch, <institute> and more
      --step                       [list] Specify starting step (only one)
                                          Available: mapping, prepare_recalibration, recalibrate, variant_calling, annotate, ControlFREEC
                                          Default: ${params.step}
      --genome                      [str] Name of iGenomes reference
                                          Default: ${params.genome}

    Main options:
      --help                       [bool] You're reading it
      --no_intervals               [bool] Disable usage of intervals
                                          Intervals are part of the genome chopped up, used to speed up preprocessing and variant calling
      --nucleotides_per_second      [int] To estimate interval size
                                          Default: ${params.nucleotides_per_second}
      --sentieon                   [bool] If sentieon is available, will enable it for Preprocessing, and Variant Calling
                                          Adds the following options for --tools: DNAseq, DNAscope and TNscope
      --skip_qc                     [str] Specify which QC tools to skip when running Sarek (multiple separated with commas)
                                          Available: all, bamQC, BaseRecalibrator, BCFtools, Documentation
                                          FastQC, MultiQC, samtools, vcftools, versions
                                          Default: None
      --target_bed                 [file] Target BED file for whole exome or targeted sequencing
                                          Default: None
      --tools                       [str] Specify tools to use for variant calling (multiple separated with commas):
                                          Available: ASCAT, CNVkit, ControlFREEC, FreeBayes, HaplotypeCaller
                                          Manta, mpileup, MSIsensor, Mutect2, Strelka, TIDDIT
                                          and/or for annotation:
                                          snpEff, VEP, merge
                                          Default: None

    Modify fastqs (trim/split):
      --trim_fastq                 [bool] Run Trim Galore
      --clip_r1                     [int] Instructs Trim Galore to remove bp from the 5' end of read 1 (or single-end reads)
      --clip_r2                     [int] Instructs Trim Galore to remove bp from the 5' end of read 2 (paired-end reads only)
      --three_prime_clip_r1         [int] Instructs Trim Galore to remove bp from the 3' end of read 1 AFTER adapter/quality trimming has been performed
      --three_prime_clip_r2         [int] Instructs Trim Galore to remove bp from the 3' end of read 2 AFTER adapter/quality trimming has been performed
      --trim_nextseq                [int] Instructs Trim Galore to apply the --nextseq=X option, to trim based on quality after removing poly-G tails
      --save_trimmed               [bool] Save trimmed FastQ file intermediates
      --split_fastq                 [int] Specify how many reads should be contained in the split fastq file
                                          Default: no split

    Preprocessing:
      --markdup_java_options        [str] Establish values for markDuplicates memory consumption
                                          Default: ${params.markdup_java_options}
      --use_gatk_spark             [bool] Enable usage of GATK Spark implementation of their tools in local mode
      --save_bam_mapped            [bool] Save Mapped BAMs
      --skip_markduplicates        [bool] Skip MarkDuplicates

    Variant Calling:
      --ascat_ploidy                [int] Use this parameter to overwrite default behavior from ASCAT regarding ploidy
                                          Requires that --ascat_purity is set
      --ascat_purity                [int] Use this parameter to overwrite default behavior from ASCAT regarding purity
                                          Requires that --ascat_ploidy is set
      --cf_coeff                    [str] Control-FREEC coefficientOfVariation
                                          Default: ${params.cf_coeff}
      --cf_ploidy                   [str] Control-FREEC ploidy
                                          Default: ${params.cf_ploidy}
      --cf_window                   [int] Control-FREEC window size
                                          Default: Disabled
      --generate_gvcf              [bool] Enable g.vcf output from GATK HaplotypeCaller
      --no_strelka_bp              [bool] Will not use Manta candidateSmallIndels for Strelka (not recommended by Best Practices)
      --pon                        [file] Panel-of-normals VCF (bgzipped) for GATK Mutect2 / Sentieon TNscope
                                          See: https://software.broadinstitute.org/gatk/documentation/tooldocs/current/org_broadinstitute_hellbender_tools_walkers_mutect_CreateSomaticPanelOfNormals.php
      --pon_index                  [file] Index of pon panel-of-normals VCF
                                          If none provided, will be generated automatically from the PON
      --ignore_soft_clipped_bases  [bool] Do not analyze soft clipped bases in the reads for GATK Mutect2
                                          Default: Do not use
      --umi                        [bool] If provided, UMIs steps will be run to extract and annotate the reads with UMI and create consensus reads
      --read_structure1          [string] When processing UMIs, a read structure should always be provided for each of the fastq files. If the read does not contain any UMI, the structure will be +T (i.e. only template of any length). 
                                          See: https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures
      --read_structure2          [string] When processing UMIs, a read structure should always be provided for each of the fastq files. If the read does not contain any UMI, the structure will be +T (i.e. only template of any length). 
                                          See: https://github.com/fulcrumgenomics/fgbio/wiki/Read-Structures

    Annotation:
      --annotate_tools              [str] Specify from which tools Sarek should look for VCF files to annotate, only for step Annotate
                                          Available: HaplotypeCaller, Manta, Mutect2, Strelka, TIDDIT
                                          Default: None
      --annotation_cache           [bool] Enable the use of cache for annotation, to be used with --snpeff_cache and/or --vep_cache
      --snpeff_cache               [file] Specity the path to snpEff cache, to be used with --annotation_cache
      --vep_cache                  [file] Specity the path to VEP cache, to be used with --annotation_cache
      --cadd_cache                 [bool] Enable CADD cache
      --cadd_indels                [file] Path to CADD InDels file
      --cadd_indels_tbi            [file] Path to CADD InDels index
      --cadd_wg_snvs               [file] Path to CADD SNVs file
      --cadd_wg_snvs_tbi           [file] Path to CADD SNVs index
      --genesplicer                [file] Enable genesplicer within VEP

    References options:
      --igenomes_base              [file] Specify base path to AWS iGenomes
                                          Default: ${params.igenomes_base}
      --igenomes_ignore            [bool] Do not use AWS iGenomes. Will load genomes.config instead of igenomes.config
      --genomes_base               [file] Specify base path to reference genome
      --save_reference             [bool] Save built references
      
    References:                           If not specified in the configuration file or you wish to overwrite any of the references.
      --ac_loci                    [file] Loci file for ASCAT
      --ac_loci_gc                 [file] Loci GC file for ASCAT
      --bwa                        [file] BWA indexes
                                          If none provided, will be generated automatically from the fasta reference
      --chr_dir                    [file] Chromosomes folder
      --chr_length                 [file] Chromosomes length file
      --dbsnp                      [file] Dbsnp file
      --dbsnp_index                [file] Dbsnp index
                                          If none provided, will be generated automatically if a dbsnp file is provided
      --dict                       [file] Fasta dictionary file
                                          If none provided, will be generated automatically from the fasta reference
      --fasta                      [file] Fasta reference
      --fasta_fai                  [file] Fasta reference index
                                          If none provided, will be generated automatically from the fasta reference
      --germline_resource          [file] Germline Resource File for GATK Mutect2
      --germline_resource_index    [file] Germline Resource Index for GATK Mutect2
                                          if none provided, will be generated automatically if a germlineResource file is provided
      --intervals                  [file] Intervals
                                          If none provided, will be generated automatically from the fasta reference
                                          Use --no_intervals to disable automatic generation
      --known_indels               [file] Known indels file
      --known_indels_index         [file] Known indels index
                                          If none provided, will be generated automatically if a knownIndels file is provided
      --mappability                [file] Mappability file for Control-FREEC
      --snpeff_db                   [str] snpEff Database version
      --species                     [str] Species for VEP
      --vep_cache_version           [int] VEP cache version

    Other options:
      --outdir                     [file] The output directory where the results will be saved
      --publish_dir_mode           [list] Mode for publishing results in the output directory (only one)
                                          Available: symlink, rellink, link, copy, copyNoFollow, move
                                          Default: copy
      --sequencing_center           [str] Name of sequencing center to be displayed in BAM file
      --multiqc_config             [file] Specify a custom config file for MultiQC
      --monochrome_logs            [bool] Logs will be without colors
      --email                       [str] Set this parameter to your e-mail address to get a summary e-mail with details of the run sent to you when the workflow exits
      --email_on_fail               [str] Same as --email, except only send mail if the workflow is not successful
      --plaintext_email            [bool] Enable plaintext email
      --max_multiqc_email_size      [str] Threshold size for MultiQC report to be attached in notification email. If file generated by pipeline exceeds the threshold, it will not be attached
                                          Default: 25MB
      -name                         [str] Name for the pipeline run. If not specified, Nextflow will automatically generate a random mnemonic

    AWSBatch options:
      --awsqueue                    [str] The AWSBatch JobQueue that needs to be set when running on AWSBatch
      --awsregion                   [str] The AWS Region for your AWS Batch job to run on
      --awscli                      [str] Path to the AWS CLI tool
    """.stripIndent()
}

def printNfcoreSarekWelcomeGraphic() {
  log.info nfcoreHeader()
}

def getSummaryMapFromParamsObjectAndArgs(step, custom_runName, skipQC, tools) {
  // also requires params and workflow maps to be initialized
  
  def summary = [:]

  if (workflow.revision) summary['Pipeline Release'] = workflow.revision

  summary['Run Name'] = custom_runName ?: workflow.runName
  summary['Max Resources'] = "${params.max_memory} memory, ${params.max_cpus} cpus, ${params.max_time} time per job"

  if (workflow.containerEngine) summary['Container']  = "${workflow.containerEngine} - ${workflow.container}"

  summary['Input'] = params.input
  summary['Step']  = step
  summary['Genome'] = params.genome

  if (params.no_intervals && step != 'annotate')  summary['Intervals'] = 'Do not use'

  summary['Nucleotides/s'] = params.nucleotides_per_second

  if (params.sentieon) summary['Sention'] = "Using Sentieon for Preprocessing and/or Variant Calling"
  if (params.skip_qc) summary['QC tools skipped']  = skipQC.join(', ')
  if (params.target_bed) summary['Target BED']  = params.target_bed
  if (params.tools) summary['Tools'] = tools.join(', ')

  if (params.trim_fastq || params.split_fastq) summary['Modify fastqs (trim/split)'] = ""

  if (params.trim_fastq) {
      summary['Fastq trim']         = "Fastq trim selected"
      summary['Trim R1']            = "${params.clip_r1} bp"
      summary['Trim R2']            = "${params.clip_r2} bp"
      summary["Trim 3' R1"]         = "${params.three_prime_clip_r1} bp"
      summary["Trim 3' R2"]         = "${params.three_prime_clip_r2} bp"
      summary['NextSeq Trim']       = "${params.trim_nextseq} bp"
      summary['Saved Trimmed Fastq'] = params.save_trimmed ? 'Yes' : 'No'
  }

  if (params.split_fastq)          summary['Reads in fastq']                   = params.split_fastq

  summary['MarkDuplicates'] = "Options"
  summary['Java options'] = params.markdup_java_options
  summary['GATK Spark']   = params.use_gatk_spark ? 'Yes' : 'No'

  summary['Save BAMs mapped']   = params.save_bam_mapped ? 'Yes' : 'No'
  summary['Skip MarkDuplicates']   = params.skip_markduplicates ? 'Yes' : 'No'

  if ('ascat' in tools) {
      summary['ASCAT'] = "Options"
      if (params.ascat_purity) summary['purity'] = params.ascat_purity
      if (params.ascat_ploidy) summary['ploidy'] = params.ascat_ploidy
  }

  if ('controlfreec' in tools) {
      summary['Control-FREEC'] = "Options"
      if (params.cf_window)    summary['window']                 = params.cf_window
      if (params.cf_coeff)     summary['coefficientOfVariation'] = params.cf_coeff
      if (params.cf_ploidy)    summary['ploidy']                 = params.cf_ploidy
  }

  if ('haplotypecaller' in tools)             summary['GVCF']       = params.generate_gvcf ? 'Yes' : 'No'
  if ('strelka' in tools && 'manta' in tools) summary['Strelka BP'] = params.no_strelka_bp ? 'No' : 'Yes'
  if (params.pon && ('mutect2' in tools || (params.sentieon && 'tnscope' in tools))) summary['Panel of normals'] = params.pon

  if (params.annotate_tools) summary['Tools to annotate'] = annotate_tools.join(', ')

  if (params.annotation_cache) {
      summary['Annotation cache'] = "Enabled"
      if (params.snpeff_cache) summary['snpEff cache'] = params.snpeff_cache
      if (params.vep_cache)    summary['VEP cache']    = params.vep_cache
  }

  if (params.cadd_cache) {
      summary['CADD cache'] = "Enabled"
      if (params.cadd_indels)  summary['CADD indels']  = params.cadd_indels
      if (params.cadd_wg_snvs) summary['CADD wg snvs'] = params.cadd_wg_snvs
  }

  if (params.genesplicer) summary['genesplicer'] = "Enabled"

  if (params.igenomes_base && !params.igenomes_ignore) summary['AWS iGenomes base'] = params.igenomes_base
  if (params.igenomes_ignore)                          summary['AWS iGenomes']      = "Do not use"
  if (params.genomes_base && !params.igenomes_ignore)  summary['Genomes base']      = params.genomes_base

  summary['Save Reference']    = params.save_reference ? 'Yes' : 'No'

  if (params.ac_loci)                 summary['Loci']                    = params.ac_loci
  if (params.ac_loci_gc)              summary['Loci GC']                 = params.ac_loci_gc
  if (params.bwa)                     summary['BWA indexes']             = params.bwa
  if (params.chr_dir)                 summary['Chromosomes']             = params.chr_dir
  if (params.chr_length)              summary['Chromosomes length']      = params.chr_length
  if (params.dbsnp)                   summary['dbsnp']                   = params.dbsnp
  if (params.dbsnp_index)             summary['dbsnpIndex']              = params.dbsnp_index
  if (params.dict)                    summary['dict']                    = params.dict
  if (params.fasta)                   summary['fasta reference']         = params.fasta
  if (params.fasta_fai)               summary['fasta index']             = params.fasta_fai
  if (params.germline_resource)       summary['germline resource']       = params.germline_resource
  if (params.germline_resource_index) summary['germline resource index'] = params.germline_resource_index
  if (params.intervals)               summary['intervals']               = params.intervals
  if (params.known_indels)            summary['known indels']            = params.known_indels
  if (params.known_indels_index)      summary['known indels index']      = params.known_indels_index
  if (params.mappability)             summary['Mappability']             = params.mappability
  if (params.snpeff_cache)            summary['snpEff cache']            = params.snpeff_cache
  if (params.snpeff_db)               summary['snpEff DB']               = params.snpeff_db
  if (params.species)                 summary['species']                 = params.species
  if (params.vep_cache)               summary['VEP cache']               = params.vep_cache
  if (params.vep_cache_version)       summary['VEP cache version']       = params.vep_cache_version

  summary['Output dir']        = params.outdir
  summary['Publish dir mode']  = params.publish_dir_mode
  if (params.sequencing_center) summary['Sequenced by'] = params.sequencing_center

  summary['Launch dir']  = workflow.launchDir
  summary['Working dir'] = workflow.workDir
  summary['Script dir']  = workflow.projectDir
  summary['User']        = workflow.userName

  if (params.multiqc_config) summary['MultiQC config'] = params.multiqc_config

  summary['Config Profile'] = workflow.profile

  if (params.config_profile_description) summary['Config Description'] = params.config_profile_description
  if (params.config_profile_contact)     summary['Config Contact']     = params.config_profile_contact
  if (params.config_profile_url)         summary['Config URL']         = params.config_profile_url

  summary['Config Files']   = workflow.configFiles.join(', ')


  if (workflow.profile.contains('awsbatch')) {
      summary['AWS Region'] = params.awsregion
      summary['AWS Queue']  = params.awsqueue
      summary['AWS CLI']    = params.awscli
  }

  if (params.email || params.email_on_fail) {
      summary['E-mail Address']    = params.email
      summary['E-mail on failure'] = params.email_on_fail
      summary['MultiQC maxsize']   = params.max_multiqc_email_size
  }

  return summary
}

def printSummaryMessage(summaryMap) {
  log.info summaryMap.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")
  log.info "-\033[2m--------------------------------------------------\033[0m-"
}

def printMutec2Warning(inputToolsList) {
  if ('mutect2' in inputToolsList && !(params.pon)) log.warn "[nf-core/sarek] Mutect2 was requested, but as no panel of normals were given, results will not be optimal"

}

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