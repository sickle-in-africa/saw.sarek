include {
    getInputStep;
    getInputTools
} from "${params.modulesDir}/sarek.nf"

step = getInputStep()
tools = getInputTools(step)

initializeParamsScope(step, tools)


process GetBWAindexes {
    // if the BWA indexes are available on igenomes, 
    // pull them, otherwise build them locally

    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/BWAIndex/${it}" : null }

    input:
        file(fasta)
        file(_bwa_)

    output:
        path("${fasta}.*"), includeInputs: true

    when: params.fasta // && 'mapping' in step

    script:
    aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
    if ( !('NULL' =~ _bwa_) )
        """
        : # do nothing. They will get automatically pulled from iGenomes
        """
    else
        """
        ${aligner} index ${fasta}
        """
}

process BuildDict {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(fasta)

    output:
        file("${fasta.baseName}.dict")

    when: !(params.dict) && params.fasta && !('annotate' in step) && !('controlfreec' in step)

    script:
    """
    gatk --java-options "-Xmx${task.memory.toGiga()}g" \
        CreateSequenceDictionary \
        --REFERENCE ${fasta} \
        --OUTPUT ${fasta.baseName}.dict
    """
}

process BuildFastaFai {
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(fasta)

    output:
        file("${fasta}.fai")

    when: !(params.fasta_fai) && params.fasta && !('annotate' in step)

    script:
    """
    samtools faidx ${fasta}
    """
}

process BuildDbsnpIndex {
    tag "${dbsnp}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(dbsnp)

    output:
        file("${dbsnp}.tbi")

    when: !(params.dbsnp_index) && params.dbsnp && ('mapping' in step || 'preparerecalibration' in step || 'controlfreec' in tools || 'haplotypecaller' in tools || 'mutect2' in tools || 'tnscope' in tools)

    script:
    """
    tabix -p vcf ${dbsnp}
    """
}

process BuildGermlineResourceIndex {
    tag "${germlineResource}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(germlineResource)

    output:
        file("${germlineResource}.tbi")

    when: !(params.germline_resource_index) && params.germline_resource && 'mutect2' in tools

    script:
    """
    tabix -p vcf ${germlineResource}
    """
}

process BuildKnownIndelsIndex {
    tag "${knownIndels}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        each file(knownIndels)

    output:
        file("${knownIndels}.tbi")

    when: !(params.known_indels_index) && params.known_indels && ('mapping' in step || 'preparerecalibration' in step)

    script:
    """
    tabix -p vcf ${knownIndels}
    """
}

process BuildPonIndex {
    tag "${pon}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(pon)

    output:
        file("${pon}.tbi")

    when: !(params.pon_index) && params.pon && ('tnscope' in tools || 'mutect2' in tools)

    script:
    """
    tabix -p vcf ${pon}
    """
}

process BuildIntervals {
    tag "${fastaFai}"

    publishDir params.outdir, mode: params.publish_dir_mode,
    saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        file(fastaFai)

    output:
        file("${fastaFai.baseName}.bed")

    when: !(params.intervals) && !('annotate' in step) && !('controlfreec' in step) 

    script:
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