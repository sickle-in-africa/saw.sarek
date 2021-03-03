include {
    getInputTools;
    hasExtension;
    isChannelActive;
    getInactiveChannelFlag
} from "${params.modulesDir}/sarek.nf"

process FastQCFQ {
    label 'FastQC'
    label 'cpus_2'

    tag "${idPatient}-${idRun}"

    publishDir "${params.outdir}/Reports/${idSample}/FastQC/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz")

    output:
        file("*.{html,zip}")

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
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_${idRun}.bam")

    output:
        file("*.{html,zip}")

    when: !('fastqc' in skipQC)

    script:
    """
    fastqc -t 2 -q ${idSample}_${idRun}.bam
    """
}

process TrimReads {
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
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz")

    output:
        file("*.{html,zip,txt}")
        tuple val(idPatient), val(idSample), val(idRun), file("${outputPaths[0]}"), file("${outputPaths[1]}")

    when:  params.trim_fastq

    script:
        outputPaths = ["${idSample}_${idRun}_R1_val_1.fq.gz", "${idSample}_${idRun}_R2_val_2.fq.gz"]
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

// UMI - STEP 1 - ANNOTATE
// the process needs to convert fastq to unmapped bam
// and while doing the conversion, tag the bam field RX with the UMI sequence
process UMIFastqToBAM {
    publishDir "${params.outdir}/Reports/${idSample}/UMI/${idSample}_${idRun}", mode: params.publish_dir_mode

    input:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_${idRun}_R1.fastq.gz"), file("${idSample}_${idRun}_R2.fastq.gz")

    output:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_umi_converted.bam")

    when: params.umi

    // tmp folder for fgbio might be solved more elengantly?

    script:
        """
        mkdir tmp

        fgbio --tmp-dir=${PWD}/tmp \
        FastqToBam \
        -i "${idSample}_${idRun}_R1.fastq.gz" "${idSample}_${idRun}_R2.fastq.gz" \
        -o "${idSample}_umi_converted.bam" \
        --read-structures ${params.read_structure1} ${params.read_structure2} \
        --sample ${idSample} \
        --library ${idSample}
        """
}

// UMI - STEP 2 - MAP THE BAM FILE
// this is necessary because the UMI groups are created based on
// mapping position + same UMI tag
process UMIMapBamFile {

    input:
        tuple val(idPatient), val(idSample), val(idRun), file(convertedBam)
        file(bwaIndex)
        file(fasta)
        file(fastaFai) 
        
    output:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_umi_unsorted.bam")

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
        tuple val(idPatient), val(idSample), val(idRun), file(alignedBam)

    output:
        file("${idSample}_umi_histogram.txt")
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_umi-grouped.bam")

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
        tuple val(idPatient), val(idSample), val(idRun), file(groupedBamFile)

    output:
        tuple val(idPatient), val(idSample), val(idRun), file("${idSample}_umi-consensus.bam"), val("null")

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

def branchReadGroupsIntoBamOrFastqChannels(readGroups) {
    result = readGroups.branch {
        bam: hasExtension(it[3], "bam")
        pairReads: !hasExtension(it[3], "bam")
        }
    return [result.bam, result.pairReads]
}

def stripSecondInputFile(inputBam) {
    return inputBam.map {
        idPatient, idSample, idRun, inputFile1, inputFile2 ->
        [idPatient, idSample, idRun, inputFile1]
    }
}

def splitReadGroups(inputPairReads) {
    if (params.split_fastq){
        return inputPairReads
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
    else {
        return inputPairReads
    }
}

def branchReadGroupsIntoPreProcessingChannels(readGroups) {
    result = readGroups.branch {
        trimming: params.trim_fastq
        umiProcessing: params.umi
        noPreProcessing: (!params.trim_fastq && !params.umi)
    }
    return [result.trimming, result.umiProcessing, result.noPreProcessing]
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