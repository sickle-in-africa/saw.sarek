include {
    getInputStep;
    getInputTools;
    isChannelActive;
    getInactiveChannelFlag
} from "${params.modulesDir}/sarek.nf"

step = getInputStep()
tools = getInputTools(step)

initializeParamsScope(step, tools)


process GetBwaIndexes {
    // if the BWA indexes are available on igenomes, 
    // pull them, otherwise build them locally

    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/BWAIndex/${it}" : null }

    input:
        path(fasta)
        path(_bwa_)

    output:
        path("${fasta}.*"), includeInputs: true

    when: params.fasta

    script:
        aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
        if ( isChannelActive(_bwa_) == true )
            """
            : # do nothing. Indexes will get automatically pulled from iGenomes
            """
        else
            """
            ${aligner} index ${fasta}
            """
}

process GetGatkDictionary {
    // if the gatk dictionary is available on igenomes,
    // pull it, otherwise build it locally
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        path(fasta)
        path(_dict_)

    output:
        path("${fasta.baseName}.dict"), includeInputs: true

    when: params.fasta

    script:
        if ( isChannelActive(_dict_) == true )
            """
            : # do nothing. Dictionary will get automatically pulled from iGenomes
            """
        else
            """
            gatk --java-options "-Xmx${task.memory.toGiga()}g" \
                CreateSequenceDictionary \
                --REFERENCE ${fasta} \
                --OUTPUT ${fasta.baseName}.dict
            """

}

process GetSamtoolsFastaIndex {
    // if the samtools fasta index (fai) is available on iGenomes,
    // pull it, otherwise build it locally

    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        path(fasta)
        path(_fastaFai_)

    output:
        path("${fasta}.fai"), includeInputs: true

    when: params.fasta

    script:
        if ( isChannelActive(_fastaFai_) == true )
            """
            : do nothing. Fasta index will be pulled from iGenomes
            """
        else
            """
            samtools faidx ${fasta}
            """
}

process GetDbsnpIndex {
    // if the user has specified an input dbsnp file, 
    // build the index, otherwise, emit an inactive channel

    tag "${dbsnp}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        path(_dbsnp_)

    output:
        path("${outputPath}")

    script:
        outputPath = "${_dbsnp_}.tbi"
        if ( isChannelActive(_dbsnp_) == true )
            """
            tabix -p vcf ${_dbsnp_}
            """
        else
            outputPath = getInactiveChannelFlag()
            """
            touch ${outputPath}
            """
}

process GetKnownIndelsIndex {
    tag "${knownIndels}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        each path(_knownIndels_)

    output:
        path("${outputPath}")

    script:
        outputPath = "${_knownIndels_}.tbi"
        if ( isChannelActive(_knownIndels_) == true )
            """
            tabix -p vcf ${_knownIndels_}
            """
        else
            outputPath = getInactiveChannelFlag()
            """
            touch ${outputPath}
            """
}

process GetIntervalsList {
    tag "${fastaFai}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        path(fastaFai)
        path(_intervalsList_)

    output:
        path("${fastaFai.baseName}.bed"), includeInputs: true

    script:
        if ( isChannelActive(_intervalsList_) == true )
            """
            : # do nothing
            """
        else 
            """
            awk -v FS='\t' -v OFS='\t' '{ print \$1, \"0\", \$2 }' ${fastaFai} > ${fastaFai.baseName}.bed
            """ 
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