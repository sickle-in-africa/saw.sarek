#!/usr/bin/env nextflow


/*
================================================================================
                         SET UP CONFIGURATION VARIABLES
================================================================================
*/

step = 'preparerecalibration'

stepList = defineStepList()



toolList = defineToolList()
tools = params.tools ? params.tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []


skipQClist = defineSkipQClist()
skipQC = params.skip_qc ? params.skip_qc == 'all' ? skipQClist : params.skip_qc.split(',').collect{it.trim().toLowerCase().replaceAll('-', '').replaceAll('_', '')} : []
if (!checkParameterList(skipQC, skipQClist)) exit 1, 'Unknown QC tool(s), see --help for more information'

annoList = defineAnnoList()
annotate_tools = params.annotate_tools ? params.annotate_tools.split(',').collect{it.trim().toLowerCase().replaceAll('-', '')} : []
if (!checkParameterList(annotate_tools,annoList)) exit 1, 'Unknown tool(s) to annotate, see --help for more information'

// Check parameters
if ((params.ascat_ploidy && !params.ascat_purity) || (!params.ascat_ploidy && params.ascat_purity)) exit 1, 'Please specify both --ascat_purity and --ascat_ploidy, or none of them'
if (params.umi && !(params.read_structure1 && params.read_structure2)) exit 1, 'Please specify both --read_structure1 and --read_structure2, when using --umi'

// Has the run name been specified by the user?
// This has the bonus effect of catching both -name and --name
custom_runName = params.name
if (!(workflow.runName ==~ /[a-z]+_[a-z]+/)) custom_runName = workflow.runName


// MultiQC
// Stage config files
ch_multiqc_config = file("${params.sarekDir}/assets/multiqc_config.yaml", checkIfExists: true)
ch_multiqc_custom_config = params.multiqc_config ? Channel.fromPath(params.multiqc_config, checkIfExists: true) : Channel.empty()
ch_output_docs = file("${params.sarekDir}/docs/output.md", checkIfExists: true)
ch_output_docs_images = file("${params.sarekDir}/docs/images/", checkIfExists: true)

// Handle input
tsvPath = null
if (params.input && (hasExtension(params.input, "tsv") || hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz"))) tsvPath = params.input
if (params.input && (hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz"))) step = "annotate"

save_bam_mapped = params.skip_markduplicates ? true : params.save_bam_mapped ? true : false

// If no input file specified, trying to get TSV files corresponding to step in the TSV directory
// only for steps preparerecalibration, recalibrate, variantcalling and controlfreec
if (!params.input && params.sentieon) {
    switch (step) {
        case 'mapping': break
        case 'recalibrate': tsvPath = "${params.outdir}/Preprocessing/TSV/sentieon_deduped.tsv"; break
        case 'variantcalling': tsvPath = "${params.outdir}/Preprocessing/TSV/sentieon_recalibrated.tsv"; break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (!params.input && !params.sentieon && !params.skip_markduplicates) {
    switch (step) {
        case 'mapping': break
        case 'preparerecalibration': tsvPath = "${params.outdir}/Preprocessing/TSV/duplicates_marked_no_table.tsv"; break
        case 'recalibrate': tsvPath = "${params.outdir}/Preprocessing/TSV/duplicates_marked.tsv"; break
        case 'variantcalling': 
          if (!params.genomes[params.genome].dbsnp && !params.genomes[params.genome].known_indels) {
            tsvPath = "${params.outdir}/Preprocessing/TSV/duplicates_marked_no_table.tsv"
          }
          else {
            tsvPath = "${params.outdir}/Preprocessing/TSV/recalibrated.tsv"
          }
          break
        case 'controlfreec': tsvPath = "${params.outdir}/VariantCalling/TSV/control-freec_mpileup.tsv"; break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (!params.input && !params.sentieon && params.skip_markduplicates) {
    switch (step) {
        case 'mapping': break
        case 'preparerecalibration': tsvPath = "${params.outdir}/Preprocessing/TSV/mapped.tsv"; break
        case 'recalibrate': tsvPath = "${params.outdir}/Preprocessing/TSV/mapped_no_duplicates_marked.tsv"; break
        case 'variantcalling': tsvPath = "${params.outdir}/Preprocessing/TSV/recalibrated.tsv"; break
        case 'controlfreec': tsvPath = "${params.outdir}/VariantCalling/TSV/control-freec_mpileup.tsv"; break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
}

inputSample = Channel.empty()
if (tsvPath) {
    tsvFile = file(tsvPath)
    switch (step) {
        case 'mapping': inputSample = extractFastq(tsvFile); break
        case 'preparerecalibration': inputSample = extractBam(tsvFile); break
        case 'recalibrate': inputSample = extractRecal(tsvFile); break
        case 'variantcalling': inputSample = extractBam(tsvFile); break
        case 'controlfreec': inputSample = extractPileup(tsvFile); break
        case 'annotate': break
        default: exit 1, "Unknown step ${step}"
    }
} else if (params.input && !hasExtension(params.input, "tsv")) {
    log.info "No TSV file"
    if (step != 'mapping') exit 1, 'No step other than "mapping" supports a directory as an input'
    log.info "Reading ${params.input} directory"
    log.warn "[nf-core/sarek] in ${params.input} directory, all fastqs are assuming to be from the same sample, which is assumed to be a germline one"
    inputSample = extractFastqFromDir(params.input)
    (inputSample, fastqTMP) = inputSample.into(2)
    fastqTMP.toList().subscribe onNext: {
        if (it.size() == 0) exit 1, "No FASTQ files found in --input directory '${params.input}'"
    }
    tsvFile = params.input  // used in the reports
} else if (tsvPath && step == 'annotate') {
    log.info "Annotating ${tsvPath}"
} else if (step == 'annotate') {
    log.info "Trying automatic annotation on files in the VariantCalling/ directory"
} else exit 1, 'No sample were defined, see --help'

(genderMap, statusMap, inputSample) = extractInfos(inputSample)

/*
================================================================================
                               CHECKING REFERENCES
================================================================================
*/

// Initialize each params in params.genomes, catch the command line first if it was defined
// params.fasta has to be the first one
params.fasta = params.genome && !('annotate' in step) ? params.genomes[params.genome].fasta ?: null : null
// The rest can be sorted
params.ac_loci = params.genome && 'ascat' in tools ? params.genomes[params.genome].ac_loci ?: null : null
params.ac_loci_gc = params.genome && 'ascat' in tools ? params.genomes[params.genome].ac_loci_gc ?: null : null
params.bwa = params.genome && params.fasta && 'mapping' in step ? params.genomes[params.genome].bwa ?: null : null
params.chr_dir = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chr_dir ?: null : null
params.chr_length = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].chr_length ?: null : null
params.dbsnp = params.genome && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || params.sentieon) ? params.genomes[params.genome].dbsnp ?: null : null
params.dbsnp_index = params.genome && params.dbsnp ? params.genomes[params.genome].dbsnp_index ?: null : null
params.dict = params.genome && params.fasta ? params.genomes[params.genome].dict ?: null : null
params.fasta_fai = params.genome && params.fasta ? params.genomes[params.genome].fasta_fai ?: null : null
params.germline_resource = params.genome && 'mutect2' in tools ? params.genomes[params.genome].germline_resource ?: null : null
params.germline_resource_index = params.genome && params.germline_resource ? params.genomes[params.genome].germline_resource_index ?: null : null
params.intervals = params.genome && !('annotate' in step) ? params.genomes[params.genome].intervals ?: null : null
params.known_indels = params.genome && ('mapping' in step || 'preparerecalibration' in step) ? params.genomes[params.genome].known_indels ?: null : null
params.known_indels_index = params.genome && params.known_indels ? params.genomes[params.genome].known_indels_index ?: null : null
params.mappability = params.genome && 'controlfreec' in tools ? params.genomes[params.genome].mappability ?: null : null
params.snpeff_db = params.genome && ('snpeff' in tools || 'merge' in tools) ? params.genomes[params.genome].snpeff_db ?: null : null
params.species = params.genome && ('vep' in tools || 'merge' in tools) ? params.genomes[params.genome].species ?: null : null
params.vep_cache_version = params.genome && ('vep' in tools || 'merge' in tools) ? params.genomes[params.genome].vep_cache_version ?: null : null

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


Channel.from(summary.collect{ [it.key, it.value] })
    .map { k,v -> "<dt>$k</dt><dd><samp>${v ?: '<span style=\"color:#999999;\">N/A</a>'}</samp></dd>" }
    .reduce { a, b -> return [a, b].join("\n            ") }
    .map { x -> """
    id: 'sarek-summary'
    description: " - this information is collected when the pipeline is started."
    section_name: 'nf-core/sarek Workflow Summary'
    section_href: 'https://github.com/nf-core/sarek'
    plot_type: 'html'
    data: |
        <dl class=\"dl-horizontal\">
            $x
        </dl>
    """.stripIndent() }
    .set { ch_workflow_summary }

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

// PREPARING CHANNELS FOR PREPROCESSING AND QC

inputBam = Channel.create()
inputPairReads = Channel.create()

if (step in ['preparerecalibration', 'recalibrate', 'variantcalling', 'controlfreec', 'annotate']) {
    inputBam.close()
    inputPairReads.close()
} else inputSample.choice(inputPairReads, inputBam) {hasExtension(it[3], "bam") ? 1 : 0}

(inputBam, inputBamFastQC) = inputBam.into(2)

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
            [idPatient, idSample, newIdRun, reads1, reads2]}
}

inputPairReads = inputPairReads.dump(tag:'INPUT')

(inputPairReads, inputPairReadsTrimGalore, inputPairReadsFastQC, inputPairReadsUMI) = inputPairReads.into(4)

if (params.umi) inputPairReads.close()
else inputPairReadsUMI.close()

if (params.trim_fastq) inputPairReads.close()
else inputPairReadsTrimGalore.close()

// STEP 0.5: QC ON READS

// TODO: Use only one process for FastQC for FASTQ files and uBAM files
// FASTQ and uBAM files are renamed based on the sample name

process FastQCFQ {
    label 'FastQC'
    label 'cpus_2'

    tag "${idPatient}-${idRun}"

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz") from inputPairReadsFastQC

    output:
        file("*.{html,zip}") into fastQCFQReport

    when: !('fastqc' in skipQC)

    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}_R1.fastq.gz ${idSample}_${idRun}_R2.fastq.gz
    """
}

process FastQCBAM {
    label 'FastQC'
    label 'cpus_2'

    tag "${idPatient}-${idRun}"

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") from inputBamFastQC

    output:
        file("*.{html,zip}") into fastQCBAMReport

    when: !('fastqc' in skipQC)

    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}.bam
    """
}

fastQCReport = fastQCFQReport.mix(fastQCBAMReport)

fastQCReport = fastQCReport.dump(tag:'FastQC')

process TrimGalore {
    label 'TrimGalore'

    tag "${idPatient}-${idRun}"

    publishDir "${params.outdir}/Reports/${idSample}/TrimGalore/${idSample}_${idRun}", mode: params.publish_dir_mode,
      saveAs: {filename ->
        if (filename.indexOf("_fastqc") > 0) "FastQC/$filename"
        else if (filename.indexOf("trimming_report.txt") > 0) "logs/$filename"
        else if (params.save_trimmed) filename
        else null
      }

    input:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz") from inputPairReadsTrimGalore

    output:
        file("*.{html,zip,txt}") into trimGaloreReport
        set idPatient, idSample, idRun, file("${idSample}_${idRun}_R1_val_1.fq.gz"), file("${idSample}_${idRun}_R2_val_2.fq.gz") into outputPairReadsTrimGalore

    when: params.trim_fastq

    script:
    // Calculate number of --cores for TrimGalore based on value of task.cpus
    // See: https://github.com/FelixKrueger/TrimGalore/blob/master/Changelog.md#version-060-release-on-1-mar-2019
    // See: https://github.com/nf-core/atacseq/pull/65
    def cores = 1
    if (task.cpus) {
      cores = (task.cpus as int) - 4
      if (cores < 1) cores = 1
      if (cores > 4) cores = 4
      }
    c_r1 = params.clip_r1 > 0 ? "--clip_r1 ${params.clip_r1}" : ''
    c_r2 = params.clip_r2 > 0 ? "--clip_r2 ${params.clip_r2}" : ''
    tpc_r1 = params.three_prime_clip_r1 > 0 ? "--three_prime_clip_r1 ${params.three_prime_clip_r1}" : ''
    tpc_r2 = params.three_prime_clip_r2 > 0 ? "--three_prime_clip_r2 ${params.three_prime_clip_r2}" : ''
    nextseq = params.trim_nextseq > 0 ? "--nextseq ${params.trim_nextseq}" : ''
    """
    trim_galore \
         --cores ${cores} \
        --paired \
        --fastqc \
        --gzip \
        ${c_r1} ${c_r2} \
        ${tpc_r1} ${tpc_r2} \
        ${nextseq} \
        ${idSample}_${idRun}_R1.fastq.gz ${idSample}_${idRun}_R2.fastq.gz

    mv *val_1_fastqc.html "${idSample}_${idRun}_R1.trimmed_fastqc.html"
    mv *val_2_fastqc.html "${idSample}_${idRun}_R2.trimmed_fastqc.html"
    mv *val_1_fastqc.zip "${idSample}_${idRun}_R1.trimmed_fastqc.zip"
    mv *val_2_fastqc.zip "${idSample}_${idRun}_R2.trimmed_fastqc.zip"
    """
}

/*
================================================================================
                            UMIs PROCESSING
================================================================================
*/

// UMI - STEP 1 - ANNOTATE
// the process needs to convert fastq to unmapped bam
// and while doing the conversion, tag the bam field RX with the UMI sequence

process UMIFastqToBAM {
    publishDir "${params.outdir}/Reports/${idSample}/UMI/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz") from inputPairReadsUMI
        val read_structure1 from ch_read_structure1
        val read_structure2 from ch_read_structure2

    output:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_umi_converted.bam") into umi_converted_bams_ch

    when: params.umi

    // tmp folder for fgbio might be solved more elengantly?

    script:
    """
    mkdir tmp

    fgbio --tmp-dir=${PWD}/tmp \
    FastqToBam \
    -i "${idSample}_${idRun}_R1.fastq.gz" "${idSample}_${idRun}_R2.fastq.gz" \
    -o "${idSample}_umi_converted.bam" \
    --read-structures ${read_structure1} ${read_structure2} \
    --sample ${idSample} \
    --library ${idSample}
    """
}

// UMI - STEP 2 - MAP THE BAM FILE
// this is necessary because the UMI groups are created based on
// mapping position + same UMI tag

process UMIMapBamFile {
    input:
        set idPatient, idSample, idRun, file(convertedBam) from umi_converted_bams_ch
        file(bwaIndex) from ch_bwa
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_umi_unsorted.bam") into umi_aligned_bams_ch

    when: params.umi

    script:
    aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
    """
    samtools bam2fq -T RX ${convertedBam} | \
    ${aligner} mem -p -t ${task.cpus} -C -M -R \"@RG\\tID:${idSample}\\tSM:${idSample}\\tPL:Illumina\" \
    ${fasta} - | \
    samtools view -bS - > ${idSample}_umi_unsorted.bam
    """
}

// UMI - STEP 3 - GROUP READS BY UMIs
// We have chose the Adjacency method, following the nice paper and blog explanation integrated in both
// UMItools and FGBIO
// https://cgatoxford.wordpress.com/2015/08/14/unique-molecular-identifiers-the-problem-the-solution-and-the-proof/
// alternatively we can define this as input for the user to choose from

process GroupReadsByUmi {
    publishDir "${params.outdir}/Reports/${idSample}/UMI/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, idRun, file(alignedBam) from umi_aligned_bams_ch

    output:
        file("${idSample}_umi_histogram.txt") into umi_histogram_ch
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_umi-grouped.bam") into umi_grouped_bams_ch

    when: params.umi

    script:
    """
    mkdir tmp

    samtools view -h ${alignedBam} | \
    samblaster -M --addMateTags | \
    samtools view -Sb - >${idSample}_unsorted_tagged.bam

    fgbio --tmp-dir=${PWD}/tmp \
    GroupReadsByUmi \
    -s Adjacency \
    -i ${idSample}_unsorted_tagged.bam \
    -o ${idSample}_umi-grouped.bam \
    -f ${idSample}_umi_histogram.txt
    """
}

// UMI - STEP 4 - CALL MOLECULAR CONSENSUS
// Now that the reads are organised by UMI groups a molecular consensus will be created
// the resulting bam file will be again unmapped and therefore can be fed into the
// existing workflow from the step mapping

process CallMolecularConsensusReads {
    publishDir "${params.outdir}/Reports/${idSample}/UMI/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, idRun, file(groupedBamFile) from umi_grouped_bams_ch

    output:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_umi-consensus.bam"), val("null") into consensus_bam_ch

    when: params.umi

    script:
    """
    mkdir tmp

    fgbio --tmp-dir=${PWD}/tmp \
    CallMolecularConsensusReads \
    -i $groupedBamFile \
    -o ${idSample}_umi-consensus.bam \
    -M 1 -S Coordinate
    """
}

// ################# END OF UMI READS PRE-PROCESSING
// from this moment on the generated uBam files can feed into the existing tools

// STEP 1: MAPPING READS TO REFERENCE GENOME WITH BWA MEM

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

    (inputPairReads, input_pair_reads_sentieon) = inputPairReads.into(2)
    if (params.sentieon) inputPairReads.close()
    else input_pair_reads_sentieon.close()
}

process MapReads {
    label 'cpus_max'

    tag "${idPatient}-${idRun}"

    input:
        set idPatient, idSample, idRun, file(inputFile1), file(inputFile2) from inputPairReads
        file(bwaIndex) from ch_bwa
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") into bamMapped
        set idPatient, val("${idSample}_${idRun}"), file("${idSample}_${idRun}.bam") into bamMappedBamQC

    when: !(params.sentieon)

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

bamMapped = bamMapped.dump(tag:'Mapped BAM')
// Sort BAM whether they are standalone or should be merged

singleBam = Channel.create()
multipleBam = Channel.create()
bamMapped.groupTuple(by:[0, 1])
    .choice(singleBam, multipleBam) {it[2].size() > 1 ? 1 : 0}
singleBam = singleBam.map {
    idPatient, idSample, idRun, bam ->
    [idPatient, idSample, bam]
}
singleBam = singleBam.dump(tag:'Single BAM')

// STEP 1': MAPPING READS TO REFERENCE GENOME WITH SENTIEON BWA MEM

process Sentieon_MapReads {
    label 'cpus_max'
    label 'memory_max'
    label 'sentieon'

    tag "${idPatient}-${idRun}"

    input:
        set idPatient, idSample, idRun, file(inputFile1), file(inputFile2) from input_pair_reads_sentieon
        file(bwaIndex) from ch_bwa
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set idPatient, idSample, idRun, file("${idSample}_${idRun}.bam") into bam_sentieon_mapped

    when: params.sentieon

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
    """
    sentieon bwa mem -K 100000000 -R \"${readGroup}\" ${extra} -t ${task.cpus} -M ${fasta} \
    ${inputFile1} ${inputFile2} | \
    sentieon util sort -r ${fasta} -o ${idSample}_${idRun}.bam -t ${task.cpus} --sam2bam -i -
    """
}

bam_sentieon_mapped = bam_sentieon_mapped.dump(tag:'Sentieon Mapped BAM')
// Sort BAM whether they are standalone or should be merged

singleBamSentieon = Channel.create()
multipleBamSentieon = Channel.create()
bam_sentieon_mapped.groupTuple(by:[0, 1])
    .choice(singleBamSentieon, multipleBamSentieon) {it[2].size() > 1 ? 1 : 0}
singleBamSentieon = singleBamSentieon.map {
    idPatient, idSample, idRun, bam ->
    [idPatient, idSample, bam]
}
singleBamSentieon = singleBamSentieon.dump(tag:'Single BAM')

// STEP 1.5: MERGING BAM FROM MULTIPLE LANES

multipleBam = multipleBam.mix(multipleBamSentieon)

process MergeBamMapped {
    label 'cpus_8'

    tag "${idPatient}-${idSample}"

    input:
        set idPatient, idSample, idRun, file(bam) from multipleBam

    output:
        set idPatient, idSample, file("${idSample}.bam") into bam_mapped_merged

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.bam ${bam}
    """
}

bam_mapped_merged = bam_mapped_merged.dump(tag:'Merged BAM')

bam_mapped_merged = bam_mapped_merged.mix(singleBam,singleBamSentieon)

(bam_mapped_merged, bam_sentieon_mapped_merged) = bam_mapped_merged.into(2)

if (!params.sentieon) bam_sentieon_mapped_merged.close()
else bam_mapped_merged.close()

bam_mapped_merged = bam_mapped_merged.dump(tag:'BAMs for MD')
bam_sentieon_mapped_merged = bam_sentieon_mapped_merged.dump(tag:'Sentieon BAMs to Index')

process IndexBamMergedForSentieon {
    label 'cpus_8'

    tag "${idPatient}-${idSample}"

    input:
        set idPatient, idSample, file("${idSample}.bam") from bam_sentieon_mapped_merged

    output:
        set idPatient, idSample, file("${idSample}.bam"), file("${idSample}.bam.bai") into bam_sentieon_mapped_merged_indexed

    script:
    """
    samtools index ${idSample}.bam
    """
}

(bam_mapped_merged, bam_mapped_merged_to_index) = bam_mapped_merged.into(2)

process IndexBamFile {
    label 'cpus_8'

    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (save_bam_mapped) "Preprocessing/${idSample}/Mapped/${it}"
            else null
        }

    input:
        set idPatient, idSample, file("${idSample}.bam") from bam_mapped_merged_to_index

    output:
        set idPatient, idSample, file("${idSample}.bam"), file("${idSample}.bam.bai") into bam_mapped_merged_indexed
        set idPatient, idSample into tsv_bam_indexed

    when: save_bam_mapped || !(params.known_indels)

    script:
    """
    samtools index ${idSample}.bam
    """
}

if (!save_bam_mapped) tsv_bam_indexed.close()

(tsv_bam_indexed, tsv_bam_indexed_sample) = tsv_bam_indexed.into(2)

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
// STEP 2: MARKING DUPLICATES

process MarkDuplicates {
    label 'cpus_16'

    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${idSample}.bam.metrics") "Reports/${idSample}/MarkDuplicates/${it}"
            else "Preprocessing/${idSample}/DuplicatesMarked/${it}"
        }

    input:
        set idPatient, idSample, file("${idSample}.bam") from bam_mapped_merged

    output:
        set idPatient, idSample, file("${idSample}.md.bam"), file("${idSample}.md.bam.bai") into bam_duplicates_marked
        set idPatient, idSample into tsv_bam_duplicates_marked
        file ("${idSample}.bam.metrics") optional true into duplicates_marked_report

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

(tsv_bam_duplicates_marked, tsv_bam_duplicates_marked_sample) = tsv_bam_duplicates_marked.into(2)

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

if ('markduplicates' in skipQC) duplicates_marked_report.close()

if (step == 'preparerecalibration') bam_duplicates_marked = inputSample

bam_duplicates_marked = bam_duplicates_marked.dump(tag:'MD BAM')
duplicates_marked_report = duplicates_marked_report.dump(tag:'MD Report')

if (params.skip_markduplicates) bam_duplicates_marked = bam_mapped_merged_indexed

(bamMD, bamMDToJoin, bam_duplicates_marked) = bam_duplicates_marked.into(3)

bamBaseRecalibrator = bamMD.combine(intBaseRecalibrator)

bamBaseRecalibrator = bamBaseRecalibrator.dump(tag:'BAM FOR BASERECALIBRATOR')

// STEP 2': SENTIEON DEDUP

process Sentieon_Dedup {
    label 'cpus_max'
    label 'memory_max'
    label 'sentieon'

    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${idSample}_*.txt" && 'sentieon' in skipQC) null
            else if (it == "${idSample}_*.txt") "Reports/${idSample}/Sentieon/${it}"
            else "Preprocessing/${idSample}/DedupedSentieon/${it}"
        }

    input:
        set idPatient, idSample, file(bam), file(bai) from bam_sentieon_mapped_merged_indexed
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set idPatient, idSample, file("${idSample}.deduped.bam"), file("${idSample}.deduped.bam.bai") into bam_sentieon_dedup

    when: params.sentieon

    script:
    """
    sentieon driver \
        -t ${task.cpus} \
        -i ${bam} \
        -r ${fasta} \
        --algo GCBias --summary ${idSample}_gc_summary.txt ${idSample}_gc_metric.txt \
        --algo MeanQualityByCycle ${idSample}_mq_metric.txt \
        --algo QualDistribution ${idSample}_qd_metric.txt \
        --algo InsertSizeMetricAlgo ${idSample}_is_metric.txt  \
        --algo AlignmentStat ${idSample}_aln_metric.txt

    sentieon driver \
        -t ${task.cpus} \
        -i ${bam} \
        --algo LocusCollector \
        --fun score_info ${idSample}_score.gz

    sentieon driver \
        -t ${task.cpus} \
        -i ${bam} \
        --algo Dedup \
        --rmdup \
        --score_info ${idSample}_score.gz  \
        --metrics ${idSample}_dedup_metric.txt ${idSample}.deduped.bam
    """
}




































// STEP 3: CREATING RECALIBRATION TABLES

process BaseRecalibrator {
    label 'cpus_1'

    tag "${idPatient}-${idSample}-${intervalBed.baseName}"

    input:
        set idPatient, idSample, file(bam), file(bai), file(intervalBed) from bamBaseRecalibrator
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnp_tbi
        file(fasta) from ch_fasta
        file(dict) from ch_dict
        file(fastaFai) from ch_fai
        file(knownIndels) from ch_known_indels
        file(knownIndelsIndex) from ch_known_indels_tbi

    output:
        set idPatient, idSample, file("${prefix}${idSample}.recal.table") into tableGatherBQSRReports
        set idPatient, idSample into recalTableTSVnoInt

    when: params.known_indels

    script:
    dbsnpOptions = params.dbsnp ? "--known-sites ${dbsnp}" : ""
    knownOptions = params.known_indels ? knownIndels.collect{"--known-sites ${it}"}.join(' ') : ""
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

if (!params.no_intervals) tableGatherBQSRReports = tableGatherBQSRReports.groupTuple(by:[0, 1])

tableGatherBQSRReports = tableGatherBQSRReports.dump(tag:'BQSR REPORTS')

if (params.no_intervals) {
    (tableGatherBQSRReports, tableGatherBQSRReportsNoInt) = tableGatherBQSRReports.into(2)
    recalTable = tableGatherBQSRReportsNoInt
} else recalTableTSVnoInt.close()

// STEP 3.5: MERGING RECALIBRATION TABLES

process GatherBQSRReports {
    label 'memory_singleCPU_2_task'
    label 'cpus_2'

    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${idSample}.recal.table" && !params.skip_markduplicates) "Preprocessing/${idSample}/DuplicatesMarked/${it}"
            else "Preprocessing/${idSample}/Mapped/${it}"
        }

    input:
        set idPatient, idSample, file(recal) from tableGatherBQSRReports

    output:
        set idPatient, idSample, file("${idSample}.recal.table") into recalTable
        file("${idSample}.recal.table") into baseRecalibratorReport
        set idPatient, idSample into recalTableTSV

    when: !(params.no_intervals)

    script:
    input = recal.collect{"-I ${it}"}.join(' ')
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GatherBQSRReports \
        ${input} \
        -O ${idSample}.recal.table \
    """
}

if ('baserecalibrator' in skipQC) baseRecalibratorReport.close()

recalTable = recalTable.dump(tag:'RECAL TABLE')

(recalTableTSV, recalTableSampleTSV) = recalTableTSV.mix(recalTableTSVnoInt).into(2)

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

    recalTableSampleTSV
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

    recalTableSampleTSV
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

bamApplyBQSR = bamMDToJoin.join(recalTable, by:[0,1])




if (step == 'recalibrate') bamApplyBQSR = inputSample

bamApplyBQSR = bamApplyBQSR.dump(tag:'BAM + BAI + RECAL TABLE')

bamApplyBQSR = bamApplyBQSR.combine(intApplyBQSR)

bamApplyBQSR = bamApplyBQSR.dump(tag:'BAM + BAI + RECAL TABLE + INT')


































// STEP 4: RECALIBRATING

process ApplyBQSR {
    label 'memory_singleCPU_2_task'
    label 'cpus_2'

    tag "${idPatient}-${idSample}-${intervalBed.baseName}"

    input:
        set idPatient, idSample, file(bam), file(bai), file(recalibrationReport), file(intervalBed) from bamApplyBQSR
        file(dict) from ch_dict
        file(fasta) from ch_fasta
        file(fastaFai) from ch_fai

    output:
        set idPatient, idSample, file("${prefix}${idSample}.recal.bam") into bam_recalibrated_to_merge

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

(bam_recalibrated_to_merge, bam_recalibrated_to_index) = bam_recalibrated_to_merge.groupTuple(by:[0, 1]).into(2)

// STEP 4': SENTIEON BQSR

bam_sentieon_dedup = bam_sentieon_dedup.dump(tag:'deduped.bam')

process Sentieon_BQSR {
    label 'cpus_max'
    label 'memory_max'
    label 'sentieon'

    tag "${idPatient}-${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "${idSample}_recal_result.csv" && 'sentieon' in skipQC) "Reports/${idSample}/Sentieon/${it}"
            else "Preprocessing/${idSample}/RecalSentieon/${it}"
        }

    input:
        set idPatient, idSample, file(bam), file(bai) from bam_sentieon_dedup
        file(dbsnp) from ch_dbsnp
        file(dbsnpIndex) from ch_dbsnp_tbi
        file(fasta) from ch_fasta
        file(dict) from ch_dict
        file(fastaFai) from ch_fai
        file(knownIndels) from ch_known_indels
        file(knownIndelsIndex) from ch_known_indels_tbi

    output:
        set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai") into bam_sentieon_recal
        set idPatient, idSample, file(bam), file(bai), file("${idSample}.recal.table") into bam_sentieon_deduped_table
        set idPatient, idSample into tsv_sentieon

    when: params.sentieon

    script:
    known = knownIndels.collect{"--known-sites ${it}"}.join(' ')
    """
    sentieon driver  \
        -t ${task.cpus} \
        -r ${fasta} \
        -i ${idSample}.deduped.bam \
        --algo QualCal \
        -k ${dbsnp} \
        ${idSample}.recal.table

    sentieon driver \
        -t ${task.cpus} \
        -r ${fasta} \
        -i ${idSample}.deduped.bam \
        -q ${idSample}.recal.table \
        --algo QualCal \
        -k ${dbsnp} \
        ${idSample}.table.post \
        --algo ReadWriter ${idSample}.recal.bam

    sentieon driver \
        -t ${task.cpus} \
        --algo QualCal \
        --plot \
        --before ${idSample}.recal.table \
        --after ${idSample}.table.post \
        ${idSample}_recal_result.csv
    """
}

(tsv_sentieon_deduped, tsv_sentieon_deduped_sample, tsv_sentieon_recal, tsv_sentieon_recal_sample) = tsv_sentieon.into(4)

// Creating a TSV file to restart from this step
tsv_sentieon_deduped.map { idPatient, idSample ->
    gender = genderMap[idPatient]
    status = statusMap[idPatient, idSample]
    bam = "${params.outdir}/Preprocessing/${idSample}/DedupedSentieon/${idSample}.deduped.bam"
    bai = "${params.outdir}/Preprocessing/${idSample}/DedupedSentieon/${idSample}.deduped.bam.bai"
    table = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.table"
    "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${table}\n"
}.collectFile(
    name: 'sentieon_deduped.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
)

tsv_sentieon_deduped_sample
    .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        bam = "${params.outdir}/Preprocessing/${idSample}/DedupedSentieon/${idSample}.deduped.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/DedupedSentieon/${idSample}.deduped.bam.bai"
        table = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.table"
        ["sentieon_deduped_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\t${table}\n"]
}

// Creating a TSV file to restart from this step
tsv_sentieon_recal.map { idPatient, idSample ->
    gender = genderMap[idPatient]
    status = statusMap[idPatient, idSample]
    bam = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.bam"
    bai = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.bam.bai"
    "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"
}.collectFile(
    name: 'sentieon_recalibrated.tsv', sort: true, storeDir: "${params.outdir}/Preprocessing/TSV"
)

tsv_sentieon_recal_sample
    .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") { idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        bam = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/RecalSentieon/${idSample}.recal.bam.bai"
        ["sentieon_recalibrated_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
}

// STEP 4.5: MERGING THE RECALIBRATED BAM FILES

process MergeBamRecal {
    label 'cpus_8'

    tag "${idPatient}-${idSample}"

    publishDir "${params.outdir}/Preprocessing/${idSample}/Recalibrated", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file(bam) from bam_recalibrated_to_merge

    output:
        set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai") into bam_recalibrated
        set idPatient, idSample, file("${idSample}.recal.bam") into bam_recalibrated_qc
        set idPatient, idSample into tsv_bam_recalibrated

    when: !(params.no_intervals)

    script:
    """
    samtools merge --threads ${task.cpus} ${idSample}.recal.bam ${bam}
    samtools index ${idSample}.recal.bam
    """
}

// STEP 4.5': INDEXING THE RECALIBRATED BAM FILES

process IndexBamRecal {
    label 'cpus_8'

    tag "${idPatient}-${idSample}"

    publishDir "${params.outdir}/Preprocessing/${idSample}/Recalibrated", mode: params.publish_dir_mode

    input:
        set idPatient, idSample, file("${idSample}.recal.bam") from bam_recalibrated_to_index

    output:
        set idPatient, idSample, file("${idSample}.recal.bam"), file("${idSample}.recal.bam.bai") into bam_recalibrated_indexed
        set idPatient, idSample, file("${idSample}.recal.bam") into bam_recalibrated_no_int_qc
        set idPatient, idSample into tsv_bam_recalibrated_no_int

    when: params.no_intervals

    script:
    """
    samtools index ${idSample}.recal.bam
    """
}

bam_recalibrated = bam_recalibrated.mix(bam_recalibrated_indexed)
bam_recalibrated_qc = bam_recalibrated_qc.mix(bam_recalibrated_no_int_qc)
tsv_bam_recalibrated = tsv_bam_recalibrated.mix(tsv_bam_recalibrated_no_int)

(bam_recalibrated_bamqc, bam_recalibrated_samtools_stats) = bam_recalibrated_qc.into(2)
(tsv_bam_recalibrated, tsv_bam_recalibrated_sample) = tsv_bam_recalibrated.into(2)

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

tsv_bam_recalibrated_sample
    .collectFile(storeDir: "${params.outdir}/Preprocessing/TSV") {
        idPatient, idSample ->
        status = statusMap[idPatient, idSample]
        gender = genderMap[idPatient]
        bam = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam"
        bai = "${params.outdir}/Preprocessing/${idSample}/Recalibrated/${idSample}.recal.bam.bai"
        ["recalibrated_${idSample}.tsv", "${idPatient}\t${gender}\t${status}\t${idSample}\t${bam}\t${bai}\n"]
}


























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
