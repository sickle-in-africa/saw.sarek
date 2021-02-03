include { nfcoreHeader } from "${params.modulesDir}/nfcore.nf"

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
def getInputTools(inputStep) {
    if (inputStep == 'controlfreec') {
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
            if (!params.genomes[params.genome].dbsnp && !params.genomes[params.genome].known_indels) {
              return "${params.outdir}/Preprocessing/TSV/duplicates_marked_no_table.tsv"
            }
            else {
              return "${params.outdir}/Preprocessing/TSV/recalibrated.tsv"
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
def getInputSampleListAsChannel(inputTsvPath, inputStep) {
  tsvFile = file(inputTsvPath)
  switch (inputStep) {
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
  if (!checkParameterList(inputAnnotationTools,annoList)) exit 1, 'Unknown tool(s) to annotate, see --help for more information'
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
  if ( !hasExtension(params.input, "tsv") ) exit 1, "your input sample tsv file has the wrong extension. Please check and try again."
  if ( hasExtension(params.input, "vcf") || hasExtension(params.input, "vcf.gz") ) exit 1, "you have specified a vcf input instead of tsv. That was perhaps meant for the annotation step."
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

def initializeDerivedParams(inputStep, inputToolsList) {

    def derivedParams = [:]

    // 'fasta' has to be the first one
    derivedParams['fasta'] = params.genome && !('annotate' in inputStep) ? params.genomes[params.genome].fasta ?: null : null
    derivedParams['ac_loci'] = params.genome && 'ascat' in inputToolsList ? params.genomes[params.genome].ac_loci ?: null : null
    derivedParams['ac_loci_gc'] = params.genome && 'ascat' in inputToolsList ? params.genomes[params.genome].ac_loci_gc ?: null : null
    derivedParams['bwa'] = params.genome && derivedParams['fasta'] && 'mapping' in inputStep ? params.genomes[params.genome].bwa ?: null : null
    derivedParams['chr_dir'] = params.genome && 'controlfreec' in inputToolsList ? params.genomes[params.genome].chr_dir ?: null : null
    derivedParams['chr_length'] = params.genome && 'controlfreec' in inputToolsList ? params.genomes[params.genome].chr_length ?: null : null
    derivedParams['dbsnp'] = params.genome && ('mapping' in inputStep || 'preparerecalibration' in inputStep || 'controlfreec' in inputToolsList || 'haplotypecaller' in inputToolsList || 'mutect2' in inputToolsList) ? params.genomes[params.genome].dbsnp ?: null : null
    derivedParams['dbsnp_index'] = params.genome && derivedParams['dbsnp'] ? params.genomes[params.genome].dbsnp_index ?: null : null
    derivedParams['dict'] = params.genome && derivedParams['fasta'] ? params.genomes[params.genome].dict ?: null : null
    derivedParams['fasta_fai'] = params.genome && derivedParams['fasta'] ? params.genomes[params.genome].fasta_fai ?: null : null
    derivedParams['germline_resource'] = params.genome && 'mutect2' in inputToolsList ? params.genomes[params.genome].germline_resource ?: null : null
    derivedParams['germline_resource_index'] = params.genome && derivedParams['germline_resource'] ? params.genomes[params.genome].germline_resource_index ?: null : null
    derivedParams['intervals'] = params.genome && !('annotate' in inputStep) ? params.genomes[params.genome].intervals ?: null : null
    derivedParams['known_indels'] = params.genome && ('mapping' in inputStep || 'preparerecalibration' in inputStep) ? params.genomes[params.genome].known_indels ?: null : null
    derivedParams['known_indels_index'] = params.genome && derivedParams['known_indels'] ? params.genomes[params.genome].known_indels_index ?: null : null
    derivedParams['mappability'] = params.genome && 'controlfreec' in inputToolsList ? params.genomes[params.genome].mappability ?: null : null
    derivedParams['snpeff_db'] = params.genome && ('snpeff' in inputToolsList || 'merge' in inputToolsList) ? params.genomes[params.genome].snpeff_db ?: null : null
    derivedParams['species'] = params.genome && ('vep' in inputToolsList || 'merge' in inputToolsList) ? params.genomes[params.genome].species ?: null : null
    derivedParams['vep_cache_version'] = params.genome && ('vep' in inputToolsList || 'merge' in inputToolsList) ? params.genomes[params.genome].vep_cache_version ?: null : null

    return derivedParams
}

def getSummaryMapFromParamsScopeAndArgs(derivedParams, step, custom_runName, skipQC, tools) {
  // requires params and workflow maps to be initialized

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